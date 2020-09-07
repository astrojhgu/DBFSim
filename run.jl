using Pkg
Pkg.activate(".")
using DBFSim
using DelimitedFiles
using Serialization
using PyPlot
using Statistics
using FFTW
using Serialization
using ThreadTools
using MPI

prefix=ARGS[1]

PayLoad=Union{Nothing, Tuple{Int, Vector{Float64}}}
outfile="$(prefix).dat"
msg_file="$(prefix)_msg.txt"

MPI.Init()
comm = MPI.COMM_WORLD
my_rank=MPI.Comm_rank(comm)
nranks=MPI.Comm_size(comm)

ant_loc=DBFSim.utils.load_ant_loc("antenna_positions.txt")


ds=DBFSim.sampler.DownSampler(256, 2048)
az, el, x, y, z=DBFSim.utils.sph_grid(1000, 500);

function write_msg(fname, s)
    open(fname, "a") do f
        write(f, s)
        write(f, "\n")
    end
end

function gen_sig(n)
    x=randn(n)
    y=fft(x)
    y./=abs.(y)
    real.(ifft(y))
end

if my_rank==0
    result=try 
        deserialize(outfile)
    catch SystemError
        Array{Union{Missing, Vector{Float64}}}(missing, size(az)...)
    end
    missing_idx=findall(reshape(map(ismissing,result), length(result)))
    
    job_list=zeros(Int, nranks-1)
    result_cnt=0

    

    for i in missing_idx
        #println("waiting...")
        write_msg(msg_file, "waiting...")
        (p, s)=MPI.recv(MPI.MPI_ANY_SOURCE, 1, comm)
        target=s.source
        buf=[i]
        MPI.Send(buf, s.source, 2, comm)
        job_list[target]=i

        writedlm("jobs.txt", job_list)
	
        if !isnothing(p)
	        global result_cnt;
	        result_cnt+=1
            write_msg(msg_file, "result fetched from $target")
            idx, payload=p
            result[idx]=payload
	    if result_cnt%(nranks-1)==0
	       serialize(outfile, result)
	    end
        end

        job_min=minimum(job_list)
        job_max=maximum(job_list)

        write_msg(msg_file, "$target $i $(az[i]) $(el[i]) $(job_min) $(job_max)")
    end
    serialize(outfile, result)
else
    MPI.send(nothing, 0, 1, comm)
    while true
        buf=[0]
        MPI.Recv!(buf, 0, 2, comm)
        az1=az[buf[1]]
        el1=el[buf[1]]
        println("$my_rank: $az1 $el1")

        output=mean(map(1:1) do i
            println("$(my_rank) $i")
            #signal=randn(4096*4096*1)
            signal=gen_sig(4096*4096*1)

            result=DBFSim.sampler.array_output(signal, ant_loc, ds, az1, el1)
            output=sum(result; dims=2)
            #signal1=DBFSim.sampler.sample(signal, ds, 0)
            #input=real.(DBFSim.Pfb.corr(signal1, signal1, 16384)[8193:end])
            output=real.(DBFSim.Pfb.corr(output, output, 16384)[8193:end])
            #output./input
            output
        end)
        println("$my_rank: finished")
        MPI.send((buf[1], output), 0, 1, comm)
        #println("$my_rank : $buf")
    end
end



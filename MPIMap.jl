module MPIMap
    using MPI:send,recv, Comm, Status, Barrier, Comm_rank, Comm_size, Send, Recv!, MPI_ANY_SOURCE
    export mpi_map

    function mpi_map(func::Function, data_in, comm::Comm; saved_result::Union{Nothing, AbstractArray{Union{Missing, U}}}=nothing, mgr_cb::Union{Nothing, Function}=nothing) where {T,U}
        my_rank=Comm_rank(comm)
        #println(my_rank)
        n_procs=Comm_size(comm)
        n_workers=n_procs-1
        data=collect(data_in)
        RT=Base.return_types(func, (eltype(data),))[1]
        
        if my_rank==0
            job_list=zeros(Int, n_workers)
            #result=similar(data)
            reply_cnt=0
            
            #result=Array{Union{Missing, RT}}(missing, size(data)...)
            result = if isnothing(saved_result)
                Array{Union{Missing, RT}}(missing, size(data)...)
            else
                #println(RT)
                #println(eltype(saved_result))
                @assert RT<:U
                @assert size(saved_result)==size(data)
                saved_result
            end
            
            missing_idx=findall(reshape(map(ismissing,result), length(result)))
            println("min_idx=",minimum(missing_idx))
            n_tasks=length(missing_idx)
            #println("waiting...")
            for i in missing_idx
                (p, s)=recv(MPI_ANY_SOURCE, 1, comm)::Tuple{Union{Nothing, Tuple{Int, RT}}, Status}
                
                reply_cnt+=1
                target=s.source
                job_list[target]=i
                job_min=minimum(job_list)
                job_max=maximum(job_list)

                send(i, target, 2, comm)

                if !isnothing(p)
                    #println("received from ",target)
                    r_idx, r=p
                    result[r_idx]=r
                    if !isnothing(mgr_cb)
                        mgr_cb(result;src=s,job_min=job_min, job_max=job_max)
                    end
                end
            end
            #println("shutting down...")
            while true
                (p, s)=recv(MPI_ANY_SOURCE, 1, comm)::Tuple{Union{Nothing, Tuple{Int, RT}}, Status}
                reply_cnt+=1
                #println("reply_cnt=", reply_cnt)
                target=s.source
                send(nothing, target, 2, comm)
                                
                if !isnothing(p)
                    r_idx, r=p
                    result[r_idx]=r
                    if !isnothing(mgr_cb)
                        mgr_cb(result)
                    end
                end
                if reply_cnt==n_tasks+n_workers
                    break
                end
            end
            Barrier(comm)
            
            #println("ntasks=", n_tasks)
            result
        else
            RT=Base.return_types(func, (eltype(data),))[1]
            result=Array{Union{Missing, RT}}(missing, size(data)...)
            send(nothing, 0, 1, comm)
            while true
                (idx, s)=recv(0, 2, comm)::Tuple{Union{Nothing, Int}, Status}
                if isnothing(idx)
                    break
                end
                output=func(data[idx])
                send((idx, output), 0, 1, comm)
            end
            Barrier(comm)
            result
        end
    end
end # module

module sampler
    using FFTW
    using .. utils:azel2xyz
    const light_speed=2.99792458e8
    const sample_interval=2.5e-9

    struct DownSampler
        ratio::Int
        filter_coeff_rev::Vector{Float64}
    end

    function DownSampler(ratio::Int, tap::Int, extra_cut=0)
        filter=ones(tap)
        filter[(tap÷2÷ratio+1-extra_cut):(end-tap÷2÷ratio+1+extra_cut)].=0
        fc=fftshift(real.(ifft(filter)))[end:-1:begin]
        DownSampler(ratio, fc)
    end

    function sample(signal, ds::DownSampler, skip::Int=0)
        first_point=1+skip
        tap=length(ds.filter_coeff_rev)
        map(first_point:ds.ratio:length(signal)-tap) do i
            sum(signal[i:i+tap-1].*ds.filter_coeff_rev)
        end
    end

    function array_output(signal::Vector, ant_loc::AbstractMatrix, ds::DownSampler, az, el)
        dir=azel2xyz(az, el)
        delays=map(eachrow(ant_loc)) do p
            Int(round(p'*dir/light_speed/sample_interval*ds.ratio))
        end
        delays.-=minimum(delays)
        println(delays)

        result=Vector(map(delays) do d
            sample(signal, ds, d)
        end)

        output_lengths=map(result) do r
            length(r)
        end

        min_len=minimum(output_lengths)
        result=map(result) do r
            r[begin:min_len]
        end

        hcat(result...)
    end

end

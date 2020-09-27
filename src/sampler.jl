module sampler
    export sample
    import Base:copy, deepcopy

    using FFTW
    using .. utils:azel2xyz
    import .. utils
    const light_speed=utils.light_speed
    const sample_interval=utils.sample_interval

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


    struct StatedDownSampler
        buffer::Vector{Float64}
        ds::DownSampler
        max_delay::Int
    end

    function copy(s::StatedDownSampler)::StatedDownSampler
        StatedDownSampler(copy(s.buffer), s.ds, s.max_delay)
    end

    function deepcopy(s::StatedDownSampler)::StatedDownSampler
        StatedDownSampler(deepcopy(s.buffer), s.ds, s.max_delay)
    end


    function StatedDownSampler(ds::DownSampler, max_delay::Int)
        #buffer_len=(2*max_delay+1)+length(ds.filter_coeff_rev)-1
        buffer=zeros(Float64, 0)
        StatedDownSampler(buffer, ds, max_delay)
    end

    function sample!(signal, dss::StatedDownSampler, delay::Int=0)
        @assert delay<=dss.max_delay
        extended_signal=[dss.buffer; signal]
        
        prefered_buf_len=2*dss.max_delay+1+length(dss.ds.filter_coeff_rev)-1
        actual_buf_len=length(extended_signal)-div(length(extended_signal)-prefered_buf_len, dss.ds.ratio)*dss.ds.ratio
        actual_buf_len=min(actual_buf_len, length(extended_signal))

        new_buffer=extended_signal[end-actual_buf_len+1:end]
        #println("actual_buf_len:", actual_buf_len,"  ", length(new_buffer))
        resize!(dss.buffer, actual_buf_len)
        dss.buffer[:]=new_buffer
        
        first_point=dss.max_delay+1
        last_point=length(extended_signal)-actual_buf_len+dss.max_delay
        
        #println("lastpoint=",last_point)
        tap=length(dss.ds.filter_coeff_rev)
        map(first_point:dss.ds.ratio:last_point) do i
            sum(extended_signal[i+delay:i+delay+tap-1].*dss.ds.filter_coeff_rev)
            #extended_signal[i+delay]
        end
    end

    function sample(signal, ds::DownSampler, skip::Int=0)
        first_point=1+skip
        tap=length(ds.filter_coeff_rev)
        map(first_point:ds.ratio:length(signal)-tap) do i
            sum(signal[i:i+tap-1].*ds.filter_coeff_rev)
        end
    end

    function stated_array_output(signal::Vector, ant_loc::AbstractMatrix, sds::AbstractVector{StatedDownSampler}, az, el; sample_interval=sample_interval)
        dir=azel2xyz(az, el)
        delays=map(zip(eachrow(ant_loc), sds)) do (p, sds1)
            df=-p'*dir/light_speed/sample_interval*sds1.ds.ratio
            d=Int(round(df))
            #println("delay: ", df)
            d
        end
        hcat(map(zip(delays, sds)) do (d, ds)
            sample!(signal, ds, d)
        end...)
    end

    function array_output(signal::Vector, ant_loc::AbstractMatrix, ds::DownSampler, az, el; sample_interval=sample_interval)
        dir=azel2xyz(az, el)
        delays=map(eachrow(ant_loc)) do p
            Int(round(p'*dir/light_speed/sample_interval*ds.ratio))
        end
        delays.-=minimum(delays)
        #println(delays)

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

module sampler
    using FFTW
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
end

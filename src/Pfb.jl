module Pfb

    using ThreadTools
    using SpecialFunctions
    using FFTW
    import Base:copy, deepcopy

    export filter_signal!, add_delay, coeff, FilterBank, Analyzer, Synthesizer, analyze!, synthesize!, Filter, corr

    default_coeffs=Dict(8=>4.853,
    10=>4.775,
    12=>5.257,
    14=>5.736,
    16=>5.856,
    18=>7.037,
    20=>6.499,
    22=>6.483,
    24=>7.410,
    26=>7.022,
    28=>7.097,
    30=>7.755,
    32=>7.452,
    48=>8.522,
    64=>9.396,
    96=>10.785,
    128=>11.5, 
    192=>11.5,
    256=>11.5)

    function rrerf(F, K, M)
        x=K*(2*M*F-0.5)
        sqrt(erfc(x)/2)
    end

    function coeff(N, L, K=-1.0)
        K=if haskey(default_coeffs, L)&&K<0.0  default_coeffs[L] else 8.0 end
        M=div(N, 2)
        F=(0:(L*M-1))/(L*M)
        A=rrerf.(F, K, M)
        N1=length(A)
        n=0:(div(N1,2)-2)
        A[N1.-n]=conj(A[n.+2])
        A[div(N1,2)+1]=0.0
        
        B=real.(ifft(A))
        B=fftshift(B)
        B./=sum(B)
        reshape(B, M, L)
    end

    struct Filter{T<:AbstractFloat, U}
        coeff_rev::Vector{T}
        initial_state::Vector{U}
    end

    function copy(f::Filter{T, U})::Filter{T, U} where {T<:AbstractFloat, U}
        Filter(copy(f.coeff_rev), copy(f.initial_state))
    end

    function deepcopy(f::Filter{T, U})::Filter{T, U} where {T<:AbstractFloat, U}
        Filter(deepcopy(f.coeff_rev), deepcopy(f.initial_state))
    end


    struct FilterBank{T<:AbstractFloat}
        filters_even::Vector{Filter{T, Complex{T}}}
        filters_odd::Vector{Filter{T, Complex{T}}}
    end

    function Filter(coeff::AbstractVector{T}, u::Type{U})::Filter{T} where {T<:AbstractFloat, U}
        initial_state=zeros(U, length(coeff)-1)
        Filter(Vector(coeff[end:-1:begin]), initial_state)
    end

    function filter_signal!(signal::AbstractVector, f::Filter{T, U})::Vector{U} where {T<:AbstractFloat, U}
        d=vcat(f.initial_state, signal)
        signal_length=length(signal)
        tap=length(f.coeff_rev)
        extended_signal=vcat(f.initial_state, signal)

        result=map(1:signal_length) do i
            sum(extended_signal[i:i+tap-1].*f.coeff_rev)
        end

        f.initial_state[:]=signal[signal_length-tap+2:end]
        result
    end


    function Analyzer(coeff::AbstractMatrix{T}) where {T<:AbstractFloat}
        n=size(coeff, 1)
        coeff=coeff[:, end:-1:begin]
        filters_even=Vector{Filter{T, Complex{T}}}(map(1:n) do i
            Filter(coeff[i, :], Complex{T})::Filter{T, Complex{T}}
        end)
        filters_odd=Vector{Filter{T, Complex{T}}}(map(1:n) do i
            Filter(coeff[i, :], Complex{T})::Filter{T, Complex{T}}
        end)
        FilterBank(filters_even, filters_odd)
    end

    function Synthesizer(coeff::AbstractMatrix{T}) where {T<:AbstractFloat}
        n=size(coeff, 1)
        filters_even=Vector{Filter{T, Complex{T}}}(map(1:n) do i
            Filter(coeff[i, :], Complex{T})::Filter{T, Complex{T}}
        end)

        filters_odd=Vector{Filter{T, Complex{T}}}(map(1:n) do i
            Filter(coeff[i, :], Complex{T})::Filter{T, Complex{T}}
        end)

        FilterBank(filters_even, filters_odd)
    end


    function analyze!(signal::AbstractVector{T}, pfb::FilterBank{T}) where {T<:AbstractFloat}
        n=length(pfb.filters_even)
        m=div(length(signal), n)
        x1=Complex{T}.(reshape(signal, n, m))
        x2=copy(x1)
        for i in 1:n
            x2[i, :].*=exp(1im*pi*(i-1)/n)
            x2[i, begin+1:2:end].*=-1.0
            x1[i, :]=filter_signal!(x1[i, :], pfb.filters_even[i])
            x2[i, :]=filter_signal!(x2[i, :], pfb.filters_odd[i])
        end
        

        x1=ifft(x1, 1)*n
        x2=ifft(x2, 1)*n
        result=Matrix{Complex{T}}(undef, n*2, m)
        result[begin:2:end, :]=x1
        result[begin+1:2:end, :] = x2
        result
    end


    function synthesize!(signal::AbstractMatrix{Complex{T}}, pfb::FilterBank{T}) where {T<:AbstractFloat}
        n=length(pfb.filters_even)
        m=size(signal, 2)
        y1=signal[begin:2:end, :]
        y2=signal[begin+1:2:end, :]
        Y1=fft(y1, 1)*n
        Y2=fft(y2, 1)*n
        for i in 1:n
            Y1[i, :]=filter_signal!(Y1[i, :], pfb.filters_even[i])
            Y2[i, :]=filter_signal!(Y2[i, :], pfb.filters_odd[i])

            Y2[i, :]*=exp(-1im*pi*(i-1)/n)
            Y2[i, begin+1:2:end].*=-1
        end

        reshape(Y1-Y2, n*m)
    end

    function add_delay(x, d)
        nch=size(x, 1)
        y=copy(x)
        freq=fftfreq(nch)
        pf=exp.(-2.0im*pi*freq*d)
        for i in 1:nch
            y[i,:].*=pf[i]
        end
        y
    end


    function corr(x, y, fold_len)
        n=length(x)÷fold_len
        x=reshape(x[begin:n*fold_len], fold_len, n)
        y=reshape(y[begin:n*fold_len], fold_len, n)
        fftshift(sum(fft(x, 1).*conj.(fft(y, 1)), dims=2))
    end

    

end # module

module PulsarSignal
    export PulsarParam, GaussianProfilePulsar, fetch_one_period, fetch_n_periods, calc_delay_s


    struct PulsarParam
        profile::Vector{Float64}
        dt::Float64
        period::Float64
    end

    function GaussianProfilePulsar(sigma::Float64, dt::Float64, period::Float64)::PulsarParam
        x=Int(round(sigma*5/dt))
        profile=exp.(-collect((-x:x)*dt).^2/(2*sigma^2)/2.)
        PulsarParam(profile, dt, period)
    end

    function fetch_one_period(p::PulsarParam)::Vector
        result=zeros(Int(round(p.period/p.dt)))
        result[1:length(p.profile)].=p.profile.*randn(length(p.profile))
        result
    end

    function fetch_n_periods(p::PulsarParam, n::Int)
        repeat(fetch_one_period(p), n)
    end

    function calc_delay_s(nu_MHz, dm_pc_per_cm3)
        1/2.410331e-4*dm_pc_per_cm3/nu_MHz^2
    end



end #module
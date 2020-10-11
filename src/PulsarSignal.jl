module PulsarSignal
    export Pulsar, GaussianProfilePulsar, fetch_one_period, fetch_n_periods, calc_delay_s, regulate, gen_1_period_signal, gen_1_period_signal_full

    using .. Utils:add_neg_freq_no_dc


    struct Pulsar{T<:AbstractFloat}
        nch_pos::Int # include both neg and pos freqs
        ch_min::Int
        dt::T # raw time domain resolution, pfb resolution=dt*nch
        dm::T
        period_dt_nch::Int
        ampl::Matrix{T}
    end

    function Pulsar(nch_pos::Int, dt::T, dm::T, period_dt_nch::Int, profile::Function, ch_min::Int, ::Type{T}=Float64)::Pulsar{T} where {T<:AbstractFloat}
        fs_MHz=1/dt/T(1e6)
        df_MHz=fs_MHz/(nch_pos*2)
        ampl=hcat(map(ch_min:nch_pos) do ch
            nu=ch*df_MHz
            delay=calc_delay_s(nu, dm)/(dt*nch_pos*2)/period_dt_nch
            profile.(regulate.((1:period_dt_nch).-delay, period_dt_nch)./period_dt_nch)
        end...)'
        ampl=[zeros(ch_min-1, period_dt_nch); ampl]
        Pulsar(nch_pos, ch_min,dt, dm, period_dt_nch, Matrix{T}(ampl))
    end

    function regulate(i, n::Int)
        1+mod(i-1, n)-div(n, 2)-1
    end

    function gen_1_period_signal(p::Pulsar{T})::Matrix{Complex{T}} where {T<:AbstractFloat}
        ampl=p.ampl
        sigma=sqrt(2)/2
        x=randn(size(ampl))*sigma+randn(size(ampl))*sigma*im
        x[end, :].=real.(x[end, :])*sqrt(2)
        #x[end, :].=0.0
        ampl.*x
    end

    function gen_1_period_signal_full(p::Pulsar{T})::Matrix{Complex{T}} where {T<:AbstractFloat}
        gen_1_period_signal(p) |> add_neg_freq_no_dc
    end

        
    function calc_delay_s(nu_MHz::T, dm_pc_per_cm3::T)::T where {T<:AbstractFloat}
        if nu_MHz!=0
            1/T(2.410331e-4)*dm_pc_per_cm3/nu_MHz^2
        else
            zero(T)
        end
    end

    function calc_t_dm_s(nu_MHz::T, dnu_MHz::T, dm::T)::T where {T<:AbstractFloat}
        T(8.3e3)*dm*dnu_MHz/nu_MHz^3
    end

end #module
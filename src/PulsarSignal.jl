module PulsarSignal
    export Pulsar, GaussianProfilePulsar, fetch_one_period, fetch_n_periods, calc_delay_s, regulate, gen_1_period_ampl, gen_1_period_signal


    struct Pulsar
        nch_pos::Int # include both neg and pos freqs
        dt::Float64 # raw time domain resolution, pfb resolution=dt*nch
        dm::Float64
        period_dt_nch::Int
    end

    function regulate(i, n::Int)
        1+mod(i-1, n)-div(n, 2)-1
    end

    function gen_1_period_ampl(p::Pulsar, profile::Function)
        fs_MHz=1/p.dt/1e6
        df_MHz=fs_MHz/(p.nch_pos*2)
        hcat(map(1:p.nch_pos) do ch
            nu=ch*df_MHz
            delay=calc_delay_s(nu, p.dm)/(p.dt*p.nch_pos*2)/p.period_dt_nch
            profile.(regulate.((1:p.period_dt_nch).-delay, p.period_dt_nch)./p.period_dt_nch)
        end...)'
    end

    function gen_1_period_signal(p::Pulsar, profile::Function)
        ampl=gen_1_period_ampl(p, profile)
        sigma=sqrt(2)/2
        x=randn(size(ampl))*sigma+randn(size(ampl))*sigma*im
        x[end, :].=real.(x[end, :])*sqrt(2)
        #x[end, :].=0.0
        ampl.*x
    end
        
    function calc_delay_s(nu_MHz, dm_pc_per_cm3)
        1/2.410331e-4*dm_pc_per_cm3/nu_MHz^2
    end


end #module
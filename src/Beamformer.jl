module Beamformer
    using .. Pfb:synthesize!, delay!, Delayer, FilterBank
    function beamformer_output!(signal::AbstractMatrix, beamformer::AbstractVector{Delayer{T}}, phase_factor::AbstractMatrix, pfb_syn::FilterBank{T}) where {T<:AbstractFloat}
        sumed=sum(map(zip(eachcol(signal), beamformer, eachcol(phase_factor))) do (s, b, pf)
            a=delay!(s, b, pf)
            println(size(s)," ", size(a))
            a
        end)
        real.(synthesize!(sumed, pfb_syn))
    end
end
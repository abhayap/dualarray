using DSP
using Random
using StaticArrays

"""Taken from here https://ccrma.stanford.edu/~jos/sasp/Example_Synthesis_1_F_Noise.html"""
function pinknoise(n::Integer)
    B = SVector{4}([0.049922035, -0.095993537, 0.050612699, -0.004408786])
    A = SVector{4}([1.0, -2.494956002, 2.017265875, -0.522189400])
    nT60 = 1430
    v = randn(n+nT60)
    x = filt(B, A, v)/0.08680587859687908
    @view x[nT60+1:end]
end

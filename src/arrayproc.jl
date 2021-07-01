using Convex
using DelimitedFiles
using FastTransforms
using LinearAlgebra
using Logging
using ProgressMeter
using SCS
using SphericalHarmonics
using Statistics

include("modalcoeff.jl")

const c = 343.0;

struct DualArray
    R
    ρ_r
    ρ_o
    θ_r
    ϕ_r
    θ_o
    ϕ_o
end

struct RigidArray
    R
    ρ
    θ
    ϕ
end

function cart2sph(x, y, z)
    θ = atan(√(x^2 + y^2), z)
    ϕ = atan(y, x)
    if ϕ < 0
        ϕ += 2π
    end
    r = √(x^2 + y^2 + z^2)
    θ, ϕ, r
end

function readtdesign(filename)
    td = readdlm(filename)
    td = reshape(td,3,:)
    td = permutedims(td)
    sph = cart2sph.(td[:,1],td[:,2],td[:,3])
    sph = permutedims(reinterpret(reshape, Float64, sph))
    θ = sph[:,1]
    ϕ = sph[:,2]
    r = sph[:,3]
    θ, ϕ, r
end

function acnorder(N::Integer)
    acn = 0:((N+1)^2-1)
    n = floor.(Int,.√acn)
    m = acn .-n.^2 .-n
    n, m
end

function Ysh(N::Integer, θ, ϕ)
    n, m = acnorder(N)
    Yres = computeYlm.(θ, ϕ, lmax = N, SHType = SphericalHarmonics.RealHarmonics())
    Ymat = reduce(hcat, Yres)
    # remove the Condon-Shortley phase
    condon = (-1).^m
    condon .* Ymat
end

function Y(N::Integer, θ, ϕ)
    n, m = acnorder(N)
    sphevaluate.(θ',ϕ',n,m)
end

function T(N::Integer, k::Number, R, ρ, Y)
    n, m = acnorder(N)
    w = bnrigid.(0:N, k*R*ρ, k*R)
    w = w[n .+ 1]
    permutedims(w .* Y)
end

function T(N::Integer, k::Number, R, ρ_r, ρ_o, Y_r, Y_o)
    T_r = T(N, k, R, ρ_r, Y_r)
    T_o = T(N, k, R, ρ_o, Y_o)
    vcat(T_r, T_o)
end

function E(T_Λ::Matrix, L::Integer, ampthresh_dB::Number)
    α = 10^(ampthresh_dB/20)
    β = 1/(2α)
    T_L = @view T_Λ[:,1:(L+1)^2]
    T_L' * inv((T_Λ * T_Λ') + (β^2 * I))
end

function E_L1(T_Λ::Matrix, L::Integer, ampthresh_dB::Number)
    α = 10^(ampthresh_dB/20)
    β = (1/(2α))^2
    N, Λ = size(T_Λ)
    Λ = floor(Int, √Λ) - 1
    A = Variable((L+1)^2, N)
    C = hcat( I((L+1)^2), zeros(((L+1)^2), (Λ+1)^2 - (L+1)^2) )
    problem = minimize(norm(vec(A * T_Λ - C), 2) + β * norm(vec(A), 1))
    solve!(problem, () -> SCS.Optimizer(verbose=false))
    A.value
end

function E(L::Integer, k::Number, R, ρ, Y)
    n, m = acnorder(L)
    w = bnrigid.(0:L, k*R*ρ, k*R)
    w = w[n .+ 1]
    (1 ./ w) .* pinv(permutedims(Y))
end

function AE(E_L, T_Λ, y_Λ)
    n = size(E_L, 1)
    y_hat = E_L * (T_Λ * y_Λ)
    y_L = @view y_Λ[1:n,:]
    T_L = @view T_Λ[:,1:n]
    y_id = E_L * (T_L * y_L)
    num = norm(y_hat - y_id)^2
    den = norm(y_id)^2
    num/den
end

function DI(E_L, T_Λ, y_Λ, y_V)
    n = size(E_L, 1)
    y_hat = E_L * (T_Λ * y_Λ)
    y_L = transpose(y_Λ[1:n,:])
    s_peak = abs(dot(y_L, y_hat))^2
    s_V = transpose(y_V) * y_hat
    df = s_peak / mean(abs.(s_V).^2)
    return df, transpose(abs.(s_V))
end

function WNG(E_L, T_Λ, y_Λ)
    n = size(E_L, 1)
    y_hat = E_L * (T_Λ * y_Λ)
    y_L = transpose(y_Λ[1:n,:])
    num = abs(dot(y_L, y_hat))^2
    den = abs(norm(y_L * E_L))^2
    num/den
end

function WNG_bw(E_L, T_Λ, y_Λ)
    n = size(E_L, 1)
    d = T_Λ * y_Λ
    W = transpose(E_L) * y_Λ[1:n,:]
    num = abs(dot(conj(W), d))^2
    den = abs(dot(conj(W), W))
    num/den
end

function sim_sma(sma::DualArray, N::Integer, θ, ϕ, ampthresh_dB, fs, nfft)
    f = (0:nfft/2) .* (fs/nfft)
    k = (2π/c) .* f
    kmax = k[end]
    Λ = floor(Int, min(85,2*kmax*sma.R*sma.ρ_o))
    @info "High SPH order is" Λ
    y_Λ = √(4π) .* Y(Λ, θ, ϕ)
    Y_r = √(4π) .* Y(Λ, sma.θ_r, sma.ϕ_r)
    Y_o = √(4π) .* Y(Λ, sma.θ_o, sma.ϕ_o)
    Hfilt = zeros(ComplexF64, (N+1)^2, size(θ,1), size(k,1))
    p = Progress(length(k); dt=0.5, desc="Computing filters: ")
    Threads.@threads for idx in 1:length(k)
        k_i = k[idx]
        T_k = (1/(4π)) .* T(Λ, k_i, sma.R, sma.ρ_r, sma.ρ_o, Y_r, Y_o)
        E_k = E(T_k, N, ampthresh_dB)
        Hfilt[:, :, idx] = E_k * (T_k * y_Λ)
        next!(p)
    end
    hfilt = copy(Hfilt)
    hfilt[:,:,end] = abs.(hfilt[:,:,end])
    hfilt = irfft(hfilt, nfft, 3)
    hfilt = fftshift(hfilt, 3)
    Hfilt, hfilt
end

function sim_sma(sma::RigidArray, N::Integer, θ, ϕ, ampthresh_dB, fs, nfft)
    f = (0:nfft/2) .* (fs/nfft)
    k = (2π/c) .* f
    kmax = k[end]
    Λ = floor(Int, min(85,2*kmax*sma.R))
    @info "High SPH order is" Λ
    y_Λ = √(4π) .* Y(Λ, θ, ϕ)
    Y_m = √(4π) .* Y(Λ, sma.θ, sma.ϕ)
    Hfilt = zeros(ComplexF64, (N+1)^2, size(θ,1), size(k,1))
    p = Progress(length(k); dt=0.5, desc="Computing filters: ")
    Threads.@threads for idx in 1:length(k)
        k_i = k[idx]
        T_k = (1/(4π)) .* T(Λ, k_i, sma.R, sma.ρ, Y_m)
        E_k = E(T_k, N, ampthresh_dB)
        Hfilt[:, :, idx] = E_k * (T_k * y_Λ)
        next(p)
    end
    hfilt = copy(Hfilt)
    hfilt[:,:,end] = abs.(hfilt[:,:,end])
    hfilt = irfft(hfilt, nfft, 3)
    hfilt = fftshift(hfilt, 3)
    Hfilt, hfilt
end

function WNG_DI(sma::DualArray, N::Integer, ampthresh_dB, fs, nfft)
    f = (0:nfft/2) .* (fs/nfft)
    k = (2π/c) .* f
    kmax = k[end]
    Λ = floor(Int, min(30,2*kmax*sma.R*sma.ρ_o))
    @info "High SPH order is" Λ
    y_Λ = √(4π) .* Y(Λ, π/2, 0)
    θ_v, ϕ_v, r_v = readtdesign("/Users/aparthy/dev/pie/spharrayproc/julia/res/des.3.240.21.txt")
    y_V = √(4π) .* Y(N, θ_v, ϕ_v)
    Y_r = √(4π) .* Y(Λ, sma.θ_r, sma.ϕ_r)
    Y_o = √(4π) .* Y(Λ, sma.θ_o, sma.ϕ_o)
    ϕ_s = (0:360).*π./180
    θ_s = repeat([π/2], size(ϕ_s,1))
    Y_s = √(4π) .* Y(N,θ_s,ϕ_s)
    y_k = √(4π) .* Y(Λ,π/2,2π/3)
    wng = zeros(size(k))
    di = zeros(size(k))
    di_hm = zeros(size(ϕ_s,1),size(k,1))
    ae = zeros(size(k))
    p = Progress(length(k); dt=0.5, desc="Computing WNG & DI: ")
    Threads.@threads for idx in 1:length(k)
        k_i = k[idx]
        T_k = (1/(4π)) .* T(Λ, k_i, sma.R, sma.ρ_r, sma.ρ_o, Y_r, Y_o)
        E_k = E_L1(T_k, N, ampthresh_dB)
        wng[idx] = WNG(E_k, T_k, y_Λ)
        di[idx], _ = DI(E_k, T_k, y_Λ, y_V)
        _, di_hm[:,idx] = DI(E_k, T_k, y_k, Y_s)
        ae[idx] = AE(E_k, T_k, y_Λ)
        next!(p)
    end
    return wng, di, di_hm, ae, f
end

function WNG_DI(sma::RigidArray, N::Integer, ampthresh_dB, fs, nfft)
    f = (0:nfft/2) .* (fs/nfft)
    k = (2π/c) .* f
    kmax = k[end]
    Λ = floor(Int, min(85,2*kmax*sma.R*sma.ρ))
    @info "High SPH order is" Λ
    y_Λ = √(4π) .* Y(Λ, π/2, 0)
    θ_v, ϕ_v, r_v = readtdesign("/Users/aparthy/dev/pie/spharrayproc/julia/res/des.3.240.21.txt")
    y_V = √(4π) .* Y(N, θ_v, ϕ_v)
    Y_m = √(4π) .* Y(Λ, sma.θ, sma.ϕ)
    ϕ_s = (0:360).*π./180
    θ_s = repeat([π/2], size(ϕ_s,1))
    Y_s = √(4π) .* Y(N,θ_s,ϕ_s)
    y_k = √(4π) .* Y(Λ,π/2,2π/3)
    wng = zeros(size(k))
    di = zeros(size(k))
    di_hm = zeros(size(ϕ_s,1),size(k,1))
    ae = zeros(size(k))
    p = Progress(length(k); dt=0.5, desc="Computing WNG & DI: ")
    Threads.@threads for idx in 1:length(k)
        k_i = k[idx]
        T_k = (1/(4π)) .* T(Λ, k_i, sma.R, sma.ρ, Y_m)
        E_k = E(T_k, N, ampthresh_dB)
        #E_k = E(N, k_i, sma.R, sma.ρ, Y_m[1:(N+1)^2,:])
        wng[idx] = WNG(E_k, T_k, y_Λ)
        di[idx], _ = DI(E_k, T_k, y_Λ, y_V)
        _, di_hm[:,idx] = DI(E_k, T_k, y_k, Y_s)
        ae[idx] = AE(E_k, T_k, y_Λ)
        next!(p)
    end
    return wng, di, di_hm, ae, f
end

function smafilters(sma::DualArray, N::Integer, ampthresh_dB, fs, nfft)
    f = (0:nfft/2) .* (fs/nfft)
    k = (2π/c) .* f
    kmax = k[end]
    Λ = floor(Int, min(85,2*kmax*sma.R*sma.ρ_o))
    @info "High SPH order is" Λ
    Y_r = √(4π) .* Y(Λ, sma.θ_r, sma.ϕ_r)
    Y_o = √(4π) .* Y(Λ, sma.θ_o, sma.ϕ_o)
    Hfilt = zeros(ComplexF64, (N+1)^2, length(sma.θ_r) + length(sma.θ_o), length(k))
    @showprogress 1 "Computing filter..." for (idx, k_i) in enumerate(k)
        T_k = (1/(4π)) .* T(Λ, k_i, sma.R, sma.ρ_r, sma.ρ_o, Y_r, Y_o)
        Hfilt[:,:,idx] = E(T_k, N, ampthresh_dB)
    end
    hfilt = copy(Hfilt)
    hfilt[:,:,end] = abs.(hfilt[:,:,end])
    hfilt = irfft(hfilt, nfft, 3)
    hfilt = fftshift(hfilt, 3)
    Hfilt, hfilt
end

function smafilters(sma::RigidArray, N::Integer, ampthresh_dB, fs, nfft)
    f = (0:nfft/2) .* (fs/nfft)
    k = (2π/c) .* f
    kmax = k[end]
    Λ = floor(Int, min(85,2*kmax*sma.R*sma.ρ))
    @info "High SPH order is" Λ
    Y_m = √(4π) .* Y(Λ, sma.θ, sma.ϕ)
    Hfilt = zeros(ComplexF64, (N+1)^2, length(sma.θ), length(k))
    @showprogress 1 "Computing filter..." for (idx, k_i) in enumerate(k)
        T_k = (1/(4π)) .* T(Λ, k_i, sma.R, sma.ρ, Y_m)
        Hfilt[:,:,idx] = E(T_k, N, ampthresh_dB)
    end
    hfilt = copy(Hfilt)
    hfilt[:,:,end] = abs.(hfilt[:,:,end])
    hfilt = irfft(hfilt, nfft, 3)
    hfilt = fftshift(hfilt, 3)
    Hfilt, hfilt
end

function gendual32x32(R=0.028, ρ_o=3.5)
    ρ_r = 1.0
    θ_r, ϕ_r, r_r = readtdesign("/Users/aparthy/dev/pie/spharrayproc/julia/res/des.3.32.7.txt")
    θ_o, ϕ_o, r_o = readtdesign("/Users/aparthy/dev/pie/spharrayproc/julia/res/des.3.32.7.txt")
    DualArray(R, ρ_r, ρ_o, θ_r, ϕ_r, θ_o, ϕ_o)
end

function gendual44x20(R=0.028, ρ_o=4.5)
    ρ_r = 1.0
    θ_r, ϕ_r, r_r = readtdesign("/Users/aparthy/dev/pie/spharrayproc/julia/res/des.3.44.8.txt")
    θ_o, ϕ_o, r_o = readtdesign("/Users/aparthy/dev/pie/spharrayproc/julia/res/des.3.20.5.txt")
    DualArray(R, ρ_r, ρ_o, θ_r, ϕ_r, θ_o, ϕ_o)
end

function genrigid32(R=0.042, ρ=1.0)
    θ, ϕ, r = readtdesign("/Users/aparthy/dev/pie/spharrayproc/julia/res/des.3.32.7.txt")
    RigidArray(R, ρ, θ, ϕ)
end

function genrigid64(R=0.042, ρ=1.0)
    θ, ϕ, r = readtdesign("/Users/aparthy/dev/pie/spharrayproc/julia/res/des.3.64.10.txt")
    RigidArray(R, ρ, θ, ϕ)
end

function genrigid4(R=0.0125, ρ=1.0)
    θ, ϕ, r = readtdesign("/Users/aparthy/dev/pie/spharrayproc/julia/res/des.3.4.2.txt")
    RigidArray(R, ρ, θ, ϕ)
end

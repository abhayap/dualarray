#%% exercise array processing
using DSP
using MAT
using Plots
using Random
using WAV

include("arrayproc.jl")
include("pinknoise.jl")


#%% compute and plot WNG for dual32x32
wng, di, di_hm, ae, f = WNG_DI(gendual32x32(), 4, 10, 40000, 1024)
# get rid of zero f
wng_3232 = wng[2:end]
di_3232 = di[2:end]
di_hm = di_hm[:,2:end]
ae = ae[2:end]
f = f[2:end]
plot(f, 10*log10.(hcat(wng_3232,di_3232,ae)), title="WNG DI AE [32x32]", label=["WNG" "DI" "AE"], legend=:bottomright, xaxis=("f (Hz)", :log), yaxis=("dB", [-40,20]), dpi=300)

#%% plot heatmap
heatmap(f, 0:360, di_hm, title="PWD Magnitude [32x32]", xaxis=("f (Hz)", :log), yaxis=("Azimuth (deg.)"), dpi=300)

#%% plot heatmap
heatmap(f/1e3, 0:360, di_hm, title="PWD Magnitude [32x32]", xaxis=("f (kHz)"), yaxis=("Azimuth (deg.)"), dpi=300)


#%% compute and plot WNG for dual44x20
wng, di, di_hm, ae, f = WNG_DI(gendual44x20(), 5, 10, 40000, 1024)
# get rid of zero f
wng_4420 = wng[2:end]
di_4420 = di[2:end]
di_hm = di_hm[:,2:end]
ae = ae[2:end]
f = f[2:end]
plot(f, 10*log10.(hcat(wng_4420,di_4420,ae)), title="WNG DI AE [44x20]", label=["WNG" "DI" "AE"], legend=:bottomright, xaxis=("f (Hz)", :log), yaxis=("dB", [-40,20]), dpi=300)

#%% plot without AE
plot(f, 10*log10.(hcat(wng,di)), title="WNG DI [44x20]", label=["WNG" "DI"], legend=:bottomright, xaxis=("f (Hz)", :log), yaxis=("dB", [-10,20]), dpi=300)

#%% plot order equivalent
di_db = 10*log10.(di)
ord_eq = 10 .^ (di_db ./ 20) .- 1
plot(f, ord_eq, label="Array Order", legend=true, linewidth=3, linecolor=:orange, xaxis=(:log), xticks=[], yaxis=([0,35]), yticks=[], background_color=:transparent, foreground_color=:black, dpi=300)

#%% plot heatmap
heatmap(f, 0:360, di_hm, title="PWD Magnitude [44x20]", xaxis=("f (Hz)", :log), yaxis=("Azimuth (deg.)"), dpi=300)

#%% plot heatmap
heatmap(f/1e3, 0:360, di_hm, title="PWD Magnitude [44x20]", xaxis=("f (kHz)"), yaxis=("Azimuth (deg.)"), dpi=300)


#%% compute and plot WNG for em64
wng, di, di_hm, ae, f = WNG_DI(genrigid64(), 6, 10, 40000, 1024)
# get rid of zero f
wng_em64 = wng[2:end]
di_em64 = di[2:end]
di_hm = di_hm[:,2:end]
ae = ae[2:end]
f = f[2:end]
plot(f, 10*log10.(hcat(wng_em64,di_em64,ae)), title="WNG DI AE [em64]", label=["WNG" "DI" "AE"], legend=:bottomright, xaxis=("f (Hz)", :log), yaxis=("dB", [-40,20]), dpi=300)

#%% plot heatmap
heatmap(f, 0:360, di_hm, title="PWD Magnitude [em64]", xaxis=("f (Hz)", :log), yaxis=("Azimuth (deg.)"), dpi=300)
#%% plot heatmap
heatmap(f/1e3, 0:360, di_hm, title="PWD Magnitude [em64]", xaxis=("f (kHz)"), yaxis=("Azimuth (deg.)"), dpi=300)

#%% compute and plot WNG for em32
wng, di, di_hm, ae, f = WNG_DI(genrigid32(), 4, 10, 40000, 1024)
# get rid of zero f
wng_em32 = wng[2:end]
di_em32 = di[2:end]
di_hm = di_hm[:,2:end]
ae = ae[2:end]
f = f[2:end]
plot(f, 10*log10.(hcat(wng_em32,di_em32,ae)), title="WNG DI AE [em32]", label=["WNG" "DI" "AE"], legend=:bottomright, xaxis=("f (Hz)", :log), yaxis=("dB", [-40,20]), dpi=300)

#%% plot heatmap
heatmap(f, 0:360, di_hm, title="PWD Magnitude [em32]", xaxis=("f (Hz)", :log), yaxis=("Azimuth (deg.)"), dpi=300)
#%% plot heatmap
heatmap(f/1e3, 0:360, di_hm, title="PWD Magnitude [em32]", xaxis=("f (kHz)"), yaxis=("Azimuth (deg.)"), dpi=300)

#%% compute and plot WNG for ambeo
wng, di, di_hm, ae, f = WNG_DI(genrigid4(), 1, 10, 40000, 1024)
# get rid of zero f
wng_amb = wng[2:end]
di_amb = di[2:end]
di_hm = di_hm[:,2:end]
ae = ae[2:end]
f = f[2:end]
plot(f, 10*log10.(hcat(wng_amb,di_amb,ae)), title="WNG DI AE [ambeo]", label=["WNG" "DI" "AE"], legend=:bottomright, xaxis=("f (Hz)", :log), yaxis=("dB", [-40,20]), dpi=300)

#%% plot heatmap
heatmap(f, 0:360, di_hm, title="PWD Magnitude [ambeo]", xaxis=("f (Hz)", :log), yaxis=("Azimuth (deg.)"), dpi=300)
#%% plot heatmap
heatmap(f/1e3, 0:360, di_hm, title="PWD Magnitude [ambeo]", xaxis=("f (kHz)"), yaxis=("Azimuth (deg.)"), dpi=300)


#%% plot WNG for all arrays
plot(f, 10*log10.(hcat(wng_amb,wng_em32,wng_em64,wng_4420)), title="WNG", label=["ambeo" "em32" "em64" "dual44x20"], legend=:bottomright, xaxis=("f (Hz)", :log), yaxis=("dB", [-10,20]), dpi=300)

#%% plot DI for all arrays
plot(f, 10*log10.(hcat(di_amb,di_em32,di_em64,di_4420)), title="Directivity Index", label=["ambeo" "em32" "em64" "dual44x20"], legend=:bottomright, xaxis=("f (Hz)", :log), yaxis=("dB", [0,20]), dpi=300)


#%% simulate HOA encoding of sound field
θ = [90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 30, 60, 90, 120, 150] * π/180
ϕ = [330, 0, 30, 60, 90, 120, 180, 240, 270, 300, 330, 300, 330, 0, 30, 60] * π/180
N = 5
ampthresh_dB = 10
fs = 48000
nfft = 1024
Hfilt, hfilt = sim_sma(gendual44x20(0.028, 3.5), N, θ, ϕ, ampthresh_dB, fs, nfft)
plot(permutedims(hfilt[:,1,:]), dpi=300)

#%% whitenoise = randn(floor(Int,5*fs))
pn = pinknoise(floor(Int,2*fs))
hoa_sig = zeros(floor(Int,0.2*fs), size(hfilt,1))
for j = 1:size(hfilt,2)
    currfilt = permutedims(hfilt[:,j,:])
    currsig = zeros(size(pn,1), size(currfilt,2))
    for (i, col) in enumerate(eachcol(currfilt))
        currsig[:,i] = fftfilt(col, pn)
    end
    hoa_sig = vcat(hoa_sig, zeros(floor(Int,0.2*fs),size(currfilt,2)))
    hoa_sig = vcat(hoa_sig, currsig)
end
hoa_sig = vcat(hoa_sig, zeros(floor(Int,0.2*fs),size(hfilt,1)))
hoa_sig = 0.8 .* hoa_sig ./ maximum(abs.(hoa_sig))
wavwrite(hoa_sig, "/Users/aparthy/dev/output/pinknoise_dual44x20_R28_r35_HOA5_N3D.wav", Fs=fs)

#%% simulate ideal HOA encoding of sound field
θ = [90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 30, 60, 90, 120, 150] * π/180
ϕ = [330, 0, 30, 60, 90, 120, 180, 240, 270, 300, 330, 300, 330, 0, 30, 60] * π/180
N = 6
fs = 48000
hfilt = Y(N, θ, ϕ)

#%% whitenoise = randn(floor(Int,5*fs))
pn = pinknoise(floor(Int,2*fs))
hoa_sig = zeros(floor(Int,0.2*fs), size(hfilt,1))
for j = 1:size(hfilt,2)
    currfilt = hfilt[:,j]
    currsig = pn .* permutedims(currfilt)
    hoa_sig = vcat(hoa_sig, zeros(floor(Int,0.2*fs),size(hfilt,1)))
    hoa_sig = vcat(hoa_sig, currsig)
end
hoa_sig = vcat(hoa_sig, zeros(floor(Int,0.2*fs),size(hfilt,1)))
hoa_sig = 0.8 .* hoa_sig ./ maximum(abs.(hoa_sig))
wavwrite(hoa_sig, "/Users/aparthy/dev/output/pinknoise_ideal_HOA6_N3D.wav", Fs=fs)


#%% test bf
k = (2π/343) * 1000
sma = gendual32x32()
Λ = 80
Y_r = √(4π) .* Y(Λ, sma.θ_r, sma.ϕ_r)
Y_o = √(4π) .* Y(Λ, sma.θ_o, sma.ϕ_o)
T_k = (1/(4π)) .* T(Λ, k, sma.R, sma.ρ_r, sma.ρ_o, Y_r, Y_o)
E_k = E(T_k, 3, 15)
p_k = T_k * Y(Λ, π/2, 0)
s = transpose(Y(3,π/2,0))*E_k*p_k
df, di_hm = DI(E_k, T_k, Y(3,π/2,0), transpose(Y(3,π/2,0)))


#%% WNG as defined in Brandstein and Ward
W = transpose(E_k)*Y(3,π/2,0);
d = p_k;
wng_ward = abs(dot(conj(W),d))^2 / abs(dot(conj(W),W))
10*log10(wng_ward)


#%% time Ylm
θ_test = rand(200);
ϕ_test = rand(200);
N = 3;
n,m = acnorder(N)
@time Y1 = computeYlm.(θ_test, ϕ_test, lmax = N, SHType = SphericalHarmonics.RealHarmonics());
@time Y3 = reduce(hcat, Y1);
@time Y2 = Y(N, θ_test, ϕ_test);

#%% Heatmap plot
xs = [string("x", i) for i = 1:10]
ys = [string("y", i) for i = 1:4]
z = float((1:4) * reshape(1:10, 1, :))
heatmap(xs, ys, z, aspect_ratio = 1)

#%% ideal beam
N = 4
ϕ = (0:360).*π./180
θ = repeat([π/2], size(ϕ,1))
Yl = Y(4,θ,ϕ)
Yk = Y(4,π/2,2π/3)
plot(ϕ, abs.(Yl'*Yk), proj = :polar, m = 2)

#%%
Ylk = repeat(Yl'*Yk, 1, 512)
heatmap(f/1e3, 0:360, Ylk, title="PWD Magnitude [ideal 4th]", xaxis=("f (kHz)"), yaxis=("Azimuth (deg.)"), dpi=300)

#%% compare filters
N = 6
ampthresh_dB = 10
fs = 48000
nfft = 3840
Hfilt, hfilt = smafilters(genrigid64(), N, ampthresh_dB, fs, nfft)


#%% compare filters dual44x20
N = 3
ampthresh_dB = 10
fs = 48000
nfft = 1024
Hfilt, hfilt = smafilters(gendual44x20(), N, ampthresh_dB, fs, nfft)
matfilts = matread(string(@__DIR__, "/test/matfilts_dual44x20.mat"))

#%% test pink noise
fs = 96000
n = 1*fs
x = pinknoise(n)
p = welch_pgram(x, 2^12; fs=fs)
plot(freq(p)./1000, 10*log10.(power(p)), xscale=:log10, xlims=(10, fs/2)./1000, ylims=(-60, 0), xlabel="Frequency (kHz)", ylabel="Power spectral density (dB/Hz)", label="Estimated PSD")
plot!(freq(p)./1000, 10*log10.(1 ./ (freq(p))), xscale=:log10, xlims=(10, fs/2)./1000, ylims=(-60, 0), xlabel="Frequency (kHz)", ylabel="Power spectral density (dB/Hz)", label="Expected roll-off")

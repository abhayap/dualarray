using SpecialFunctions

@enum Sphere open rigid cardioid

besseljs(n, x) = √(π/2x)*besselj(n+1/2, x)

besseljsd(n, x) = (n/x)*besseljs(n, x) - besseljs(n+1, x)

besselhs(n, x) = √(π/2x)*besselh(n+1/2, x)

besselhsd(n, x) = (n/x)*besselhs(n, x) - besselhs(n+1, x)

function bnopen(n, kr)
    return 4π*im^n * besseljs(n,kr)
end

function bnrigid(n, kr, ka)
    if kr ≈ 0.0 atol=1e-14
        if n == 0
            return 4π
        else
            return 0
        end
    else
        return 4π*im^n * ( besseljs(n,kr) - ( besseljsd(n,ka) / conj(besselhsd(n,ka)) ) * conj(besselhs(n,kr)) )
    end
end

function bncardioid(n, kr)
    return 4π*im^n * (besseljs(n,kr) - im*besseljsd(n,kr));
end

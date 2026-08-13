

# elasticae

using Elliptic
 
# Analytic curvature κ(s) from equation (22) — Pfefferle et al. (2018)
# κ0 - maximum curvature
# p - elliptic modulus, must satify p <= w (controls amplitude)
# w - shape parameter, must satisfy p<<w<<1 (controls frequency)
function kappa(κ0, p, w, s)
    t = (κ0 / 2w) * s
    return sqrt(κ0^2 * (1 - (p^2/w^2) * (Elliptic.Jacobi.sn(t, p^2))^2))
end
 
# Analytic torsion τ(s) from equation (17) — Pfefferle et al. (2018)
# c - Lagrange multiplier for the torsion constraint
# λ2 - Lagrange multiplier for the closure constraint
function tau(c, λ2, κ0, p, w, s)
    return (c / kappa(κ0, p, w, s)^2 + λ2) / 2
end
 
# Closure condition A(p) from equation (50) — solves R² = A(p) for p
function A(p)
    return 2 * Elliptic.E(p^2) / Elliptic.K(p^2) - 1
end


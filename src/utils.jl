
using DifferentialEquations, Plots, Elliptic, Roots, Optim, LinearAlgebra, QuadGK, ForwardDiff

# Frenet-Serret ODE system
function frenet_serret(du, u, p, s)
    κ, τ = p[1](s), p[2](s)

    x = u[1:3] # Position
    t = u[4:6] # Tangent
    n = u[7:9] # Normal
    b = u[10:12] # Binormal

    # dx/ds = t
    du[1:3] = t
    # dt/ds = κ*n
    du[4:6] = κ * n
    # dn/ds = -κ*t + τ*b
    du[7:9] = -κ * t + τ * b
    # db/ds = -τ*n
    du[10:12] = -τ * n
end

# curvature function kappa 
function kappa(κ0, p, w, s)
    t = (κ0/2w)*s
    return sqrt(κ0^2*(1-(p^2/w^2) *(Elliptic.Jacobi.sn(t,p^2))^2))
end

# torsion τ = (c/κ^2 + λ2 )/2
function tau(c, λ2, κ0, p, w, s)
    return (c/(kappa(κ0, p, w, s))^2 + λ2 )/2
end

# closure condition function, solves R^2 = A(p) for p 
function A(p)
    return 2 * Elliptic.E(p) / Elliptic.K(p) - 1
end



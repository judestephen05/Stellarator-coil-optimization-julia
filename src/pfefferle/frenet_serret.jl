
using DifferentialEquations, Plots, Elliptic
include("elasticae.jl")
 
# Frenet-Serret ODE system 
function frenet_serret(du, u, p, s)
    κ, τ = p[1](s), p[2](s)
    t = u[4:6]; n = u[7:9]; b = u[10:12]
    du[1:3] = t
    du[4:6] = κ * n
    du[7:9] = -κ * t + τ * b
    du[10:12] = -τ * n
end

# initial frame, curve starts at origin pointing in x direction
x0 = [0.0, 0.0, 0.0]
t0 = [1.0, 0.0, 0.0]
n0 = [0.0, 1.0, 0.0]
b0 = [0.0, 0.0, 1.0]
u0 = vcat(x0, t0, n0, b0)
 
κ0 = 1.0 # maximum curvature
l  = 0.8590380923408762 # (p) elliptic modulus
w  = 0.8776938213822927 # shape parameter
c  = sign(0.18) * (κ0^3 / w^2) * sqrt((w^2 - l^2) * (1 - w^2)) # Lagrange multiplier for the torsion constraint
λ2 = 0.33 / w # Lagrange multiplier for the clousre constraint 

# test parameters, p=0 gives constant curvature, produces near-circular curve
κ_func = s -> kappa(1, 0.8590380923408762, 0.8776938213822927, s)
τ_func = s -> tau(c, λ2, 1.0, 0.8590380923408762, 0.8776938213822927, s)
p = (κ_func, τ_func)

L = 3 * 4 * 0.8776938213822927 * Elliptic.K(0.8590380923408762^2) / 1.0
sspan = (0.0, L)
#sspan = (0.0, 3 * 4 * 0.36711 * Elliptic.K(.225101^2) / 1)
 
prob = ODEProblem(frenet_serret, u0, sspan, p)
sol = solve(prob, Tsit5(), saveat=L/300, abstol=1e-10, reltol=1e-10)
 
x_vals = [sol.u[i][1] for i in 1:length(sol.u)]
y_vals = [sol.u[i][2] for i in 1:length(sol.u)]
z_vals = [sol.u[i][3] for i in 1:length(sol.u)]
 
display(plot3d(x_vals, y_vals, z_vals))
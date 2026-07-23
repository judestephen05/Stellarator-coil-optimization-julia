
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
 
# test parameters, p=0 gives constant curvature, produces near-circular curve
κ_func = s -> kappa(1.3, 0, 0.99, s)
τ_func = s -> tau(0, 0, 1.0, 0, 0.99, s)
p = (κ_func, τ_func)
sspan = (0.0, 2π)
 
prob = ODEProblem(frenet_serret, u0, sspan, p)
sol = solve(prob, saveat=0.01)
 
x_vals = [sol.u[i][1] for i in 1:length(sol.u)]
y_vals = [sol.u[i][2] for i in 1:length(sol.u)]
z_vals = [sol.u[i][3] for i in 1:length(sol.u)]
 
display(plot3d(x_vals, y_vals, z_vals))

using Optim, QuadGK, Plots, ADTypes
include("fourier_curves.jl")
 
# Minimizes bending energy E = ∫κ²ds subject to:
#   - fixed integrated torsion Tw
#   - fixed total arc length L
# Uses ForwardDiff for exact gradients (much faster than numerical differentiation)
 
function objective(p)
    penalty = 1000.0
    Tw_target = 0.5
    L_target = 2π
 
    E,  _ = quadgk(t -> κ_sq(t, p) * ds(t, p), 0, 2π, rtol=1e-2, maxevals=50)
    Tw, _ = quadgk(t -> τ_geom(t, p) * ds(t, p), 0, 2π, rtol=1e-2, maxevals=50)
    L,  _ = quadgk(t -> ds(t, p), 0, 2π, rtol=1e-2, maxevals=50)
    Tw /= 2π
 
    return E + penalty*(Tw - Tw_target)^2 + penalty*(L - L_target)^2
end
 
result = optimize(objective, params, BFGS(), autodiff=AutoForwardDiff())
opt_p = Optim.minimizer(result)
 
t_vals = range(0, 2π, length=200)
x_vals = [curve_x(t, opt_p) for t in t_vals]
y_vals = [curve_y(t, opt_p) for t in t_vals]
z_vals = [curve_z(t, opt_p) for t in t_vals]
 
display(plot3d(x_vals, y_vals, z_vals, title="Optimized curve — minimum bending energy"))
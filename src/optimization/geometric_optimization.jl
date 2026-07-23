
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
 
result = optimize(objective, params, BFGS(), Optim.Options(store_trace=true, extended_trace=true), autodiff=AutoForwardDiff())
history = [tr.value for tr in Optim.trace(result)]
opt_p = Optim.minimizer(result)

println("Iteration    Cost           Cost reduction    Optimality")
let prev_cost = NaN
    for tr in Optim.trace(result)
        cost = tr.value
        reduction = isnan(prev_cost) ? NaN : prev_cost - cost
        optimality = tr.g_norm
        println(rpad(tr.iteration, 12), rpad(round(cost, digits=6), 16),
                rpad(isnan(reduction) ? "—" : round(reduction, digits=6), 18),
                round(optimality, digits=6))
        prev_cost = cost
    end
end

#=
println("Iter   E (bending)   Tw (torsion)   L (length)   Total Cost")
for (i, tr) in enumerate(Optim.trace(result))
    p_i = tr.metadata["x"]  # parameters at this iteration
    E_i,  _ = quadgk(t -> κ_sq(t, p_i)   * ds(t, p_i), 0, 2π, rtol=1e-2, maxevals=50)
    Tw_i, _ = quadgk(t -> τ_geom(t, p_i) * ds(t, p_i), 0, 2π, rtol=1e-2, maxevals=50)
    Tw_i /= 2π
    L_i,  _ = quadgk(t -> ds(t, p_i),                   0, 2π, rtol=1e-2, maxevals=50)
    println(rpad(i,6), rpad(round(E_i,digits=4),14), rpad(round(Tw_i,digits=4),15),
            rpad(round(L_i,digits=4),12), round(tr.value,digits=4))
end
=#

t_vals = range(0, 2π, length=200)
x_vals = [curve_x(t, opt_p) for t in t_vals]
y_vals = [curve_y(t, opt_p) for t in t_vals]
z_vals = [curve_z(t, opt_p) for t in t_vals]
 
display(plot3d(x_vals, y_vals, z_vals, title="Optimized curve — minimum bending energy"))

# plot history
display(plot(history, xlabel="iteration", ylabel="objective value", title="Convergence history", yscale=:log10))
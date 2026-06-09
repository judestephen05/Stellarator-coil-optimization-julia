
using DifferentialEquations, Plots, Elliptic, Roots
include("elasticae.jl")
include("frenet_serret.jl")
include("unit_disk.jl")
 
# Finds closed curve families on the ΔΦ = -1/3 contour.
# For each valid (X,Y) point integrates Frenet-Serret over 3 fundamental lengths.
# Based on Section VI of Pfefferle et al. (2018).
 
# collect (X, Y, p, w) for points near the ΔΦ = -1/3 contour
contour_points = []
for (i, X) in enumerate(x_grid)
    for (j, Y) in enumerate(y_grid)
        if X^2 + Y^2 <= 1
            if abs(ΔΦ[i,j] - (-1/3)) < 0.001
                push!(contour_points, (X, Y, p_vals[i,j], w_vals[i,j]))
            end
        end
    end
end
 
fig = plot3d()
for pt in contour_points
    X, Y, p, w = pt
    λ2 = X / w
    κ0 = 1.0
 
    # c from equation (24)
    c = (κ0^3 / w^2) * sqrt((w^2 - p^2) * (1 - w^2))
 
    κ_funcnew = s -> kappa(κ0, p, w, s)
    τ_funcnew = s -> tau(c, λ2, κ0, p, w, s)
    ode_p = (κ_funcnew, τ_funcnew)
 
    # 3 fundamental lengths from eq (46): l = 4wK(p²)/κ0
    sspan = (0.0, 3 * 4 * w * Elliptic.K(p^2))
 
    probnew = ODEProblem(frenet_serret, u0, sspan, ode_p)
    solnew = solve(probnew)
 
    xv = [solnew.u[k][1] for k in 1:length(solnew.u)]
    yv = [solnew.u[k][2] for k in 1:length(solnew.u)]
    zv = [solnew.u[k][3] for k in 1:length(solnew.u)]
 
    plot3d!(fig, xv, yv, zv)
 
    # print diagnostics for the curve closest to the QIV edge point (0.943, -0.333)
    if sqrt((X - 0.943)^2 + (Y + 0.333)^2) < 0.1
        println("sspan = ", sspan)
        println("p = ", p, "  w = ", w)
        println("start: ", solnew.u[1][1:3])
        println("end:   ", solnew.u[end][1:3])
    end
end
 
display(length(contour_points))
display(fig)

using DifferentialEquations, Plots, Elliptic, Roots, Statistics
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
 
closure_pos_err = Float64[]
closure_tan_err = Float64[]
c_vals          = Float64[]
λ2_vals         = Float64[]
κ0_vals = Float64[]


fig = plot3d()
for pt in contour_points
    X, Y, p, w = pt
    λ2 = X / w
    κ0 = 1.0
 
    # c from equation (24)
    c = sign(Y) * (κ0^3 / w^2) * sqrt((w^2 - p^2) * (1 - w^2))
 
    κ_funcnew = s -> kappa(κ0, p, w, s)
    τ_funcnew = s -> tau(c, λ2, κ0, p, w, s)
    ode_p = (κ_funcnew, τ_funcnew)

    push!(c_vals,   c)
    push!(λ2_vals,  λ2)
    push!(κ0_vals,  κ0)
    
    # 3 fundamental lengths from eq (46): l = 4wK(p²)/κ0
    sspan = (0.0, 3 * 4 * w * Elliptic.K(p^2))
 
    probnew = ODEProblem(frenet_serret, u0, sspan, ode_p)
    solnew = solve(probnew, saveat=sspan[2]/300)
 
    xv = [solnew.u[k][1] for k in 1:length(solnew.u)]
    yv = [solnew.u[k][2] for k in 1:length(solnew.u)]
    zv = [solnew.u[k][3] for k in 1:length(solnew.u)]
 
    plot3d!(fig, xv, yv, zv)

    pos_err = norm(solnew.u[end][1:3] - solnew.u[1][1:3])
    tan_err = norm(solnew.u[end][4:6] - solnew.u[1][4:6])
    push!(closure_pos_err, pos_err)
    push!(closure_tan_err, tan_err)
 
    #= print diagnostics for the curve closest to the QIV edge point (0.943, -0.333)
    if sqrt((X - 0.943)^2 + (Y + 0.333)^2) < 0.1
        println("sspan = ", sspan)
        println("p = ", p, "  w = ", w)
        println("start: ", solnew.u[1][1:3])
        println("end:   ", solnew.u[end][1:3])
    end=#
end
 

println("median position error: ", median(closure_pos_err))
println("median tangent error:  ", median(closure_tan_err))
println("worst position error:  ", maximum(closure_pos_err))

display(length(contour_points))
display(fig)


for (idx, pt) in enumerate(contour_points)
    X, Y, p, w = pt
    # recompute pos_err for this point same as in your loop
    println(idx, "  X=", X, " Y=", Y, " κ0=", κ0_vals[idx]," p=", p, " w=", w, " c=", c_vals[idx], " λ2=", λ2_vals[idx], "  err=", closure_pos_err[idx])
end
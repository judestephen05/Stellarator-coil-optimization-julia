# Closure conditions (Section VI, Pfefferle et al. 2018):
# 1. Zero vertical displacement: A(p) = 2E(p)/K(p) - 1 = R² - Eqn (50)
#    determines p² as a function of (X,Y), restricts parameters to unit disk
# 2. Rational toroidal displacement: ΔΦ/2πn = -1/n (m=1) - Eqn (52)
#    curves that close after exactly 1 toroidal revolution in n field periods
#    restricts parameters further to the ΔΦ = -1/n contour within the unit disk


using Plots, Elliptic, Roots, DifferentialEquations
 
# Computes U, V and ΔΦ over the unit disk X²+Y²≤1.
# Reproduces Figures 1, 2, 3 from Pfefferle
 
x_grid = -1.0:0.01:1.0
y_grid = -1.0:0.01:1.0
 
R = zeros(length(x_grid), length(y_grid)) # Radial coordinate R = √(X²+Y²) - Eqn (29)
p_vals = zeros(length(x_grid), length(y_grid)) # Parametization solutions - 0 <= p <= w <= 1, solved from closure condition A(p) = R²
w_vals = zeros(length(x_grid), length(y_grid)) # Parametization solutions w = √(Y²+p²) - Eqn (30), satisfies p ≤ w ≤ 1
J = zeros(length(x_grid), length(y_grid)) # Normalized Noether charge J̃ - Eqn (32)
U = zeros(length(x_grid), length(y_grid)) # Transformed coordinate X - Eqn (33)
V = zeros(length(x_grid), length(y_grid)) # Conjugate of U - Eqn (34) also U² + V² = X² + Y² = R²
M = zeros(length(x_grid), length(y_grid)) # Eqn (40)
ΔΦ = zeros(length(x_grid), length(y_grid)) # Toroidal displacement per segment Δϕ/2πn - Eqn (52) 

# Closure condition Eqn (50) - solves R² = A(p) for p
A = p -> 2 * Elliptic.E(p^2) / Elliptic.K(p^2) - 1
n = 4
contour_points = []

# loop over every point in the unit disk 
for (i, x) in enumerate(x_grid)
    for (j, y) in enumerate(y_grid)
        if x^2 + y^2 <= 1 
            R[i, j] = sqrt(x^2 + y^2) # Eqn (29)
            p_vals[i, j] = find_zero(p -> A(p) - R[i,j]^2, (0, 1)) # Eqn (49) and solve for 0
            w_vals[i, j] = sqrt(y^2 + p_vals[i,j]^2) # Eqn (30)
            wp = sqrt(1 - w_vals[i,j]^2)
            pp = sqrt(1 - p_vals[i,j]^2) 
 
            J[i, j] = sqrt((1 - R[i,j]^2)^2 + 4*(y*wp - x*w_vals[i,j])^2) / (4*w_vals[i,j]^2)
            U[i, j] = (x*(1 - R[i,j]^2) + 2*(y*wp - x*w_vals[i,j])*w_vals[i,j]) / (4*J[i,j]*w_vals[i,j]^2)
            V[i, j] = (y*(1 - R[i,j]^2) - 2*(y*wp - x*w_vals[i,j])*wp) / (4*J[i,j]*w_vals[i,j]^2)
 
            M[i, j] = p_vals[i,j]^2 / (p_vals[i,j]^2 + V[i,j]^2)
            N = U[i,j] * ((1 - R[i,j]^2) / (2*(p_vals[i,j]^2 + V[i,j]^2)) - 1) - (4*w_vals[i,j]^2 * J[i,j] * x) / (2*(p_vals[i,j]^2 + V[i,j]^2))
 
            η = acos(sqrt(V[i,j]^2 / (1 - p_vals[i,j]^2))) # angle for incomplete elliptic integrals - Eqn (53)
            
 
            K_val = Elliptic.K(p_vals[i,j]^2)
            E_complete = Elliptic.E(p_vals[i,j]^2) # complete elliptic integral of second kind
            E_inc = Elliptic.E(η, pp^2) # incomplete elliptic integral of second kind
            F_inc = Elliptic.F(η, pp^2) # incomplete elliptic integral of first kind
 
            ΔΦ[i, j] = (1/2) * ((U[i,j] * K_val)/(pi/2) + (sign(N)/(pi/2)) * (K_val*E_inc - (K_val - E_complete)*F_inc))

            if abs(ΔΦ[i,j] - (-1/n)) < 0.001 # Keep values within 0.001 of -1/n (ΔΦ = -1/n where n = 3)
                push!(contour_points, (x, y, p_vals[i,j], w_vals[i,j]))
            end

        end
    end
end
 
display(contour(x_grid, y_grid, U', title="U contours (Fig 1)"))
display(contour(x_grid, y_grid, V', title="V contours (Fig 2)"))
display(contour(x_grid, y_grid, ΔΦ', title="ΔΦ/2πn contours (Fig 3)"))



# curvature κ(s) from Eqn (22)
# κ0 - maximum curvature
# p - elliptic modulus, must satify p <= w (controls amplitude)
# w - shape parameter, must satisfy p<<w<<1 (controls frequency)
function kappa(κ0, p, w, s)
    t = (κ0 / 2w) * s
    return sqrt(κ0^2 * (1 - (p^2/w^2) * (Elliptic.Jacobi.sn(t, p^2))^2))
end
 
# torsion τ(s) from Eqn (17) 
# c - Lagrange multiplier for the torsion constraint
# λ2 - Lagrange multiplier for the closure constraint
function tau(c, λ2, κ0, p, w, s)
    return (c / kappa(κ0, p, w, s)^2 + λ2) / 2
end



# Frenet-Serret ODE system 
# u = [x, y, z, Tx, Ty, Tz, Nx, Ny, Nz, Bx, By, Bz]
# κ(s) and τ(s) are functions passed in (creating the shape of the coil)
function frenet_serret(du, u, p, s)
    κ, τ = p[1](s), p[2](s)

    T = u[4:6]
    N = u[7:9]
    B = u[10:12]

    du[1:3] = T         # dx/ds = T
    du[4:6] = κ*N       # dT/ds = κN
    du[7:9] = -κ*T + τ*B      # dN/ds = -κT + τB
    du[10:12] = -τ*N          # dB/ds = -τN
end

# initial frame, curve starts at origin pointing in x direction
u0 = [0.0,0.0,0.0, 1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0]


c_vals = Float64[]
λ2_vals = Float64[]
κ0_vals = Float64[]

fig = plot3d()
# Finding the lagrange parameters and solving the ODE for each contour point
for i in contour_points
    x, y, p, w = i

    κ0 = 100 # coefficient of curvature 
    λ2 = x*κ0/w # Eqn (28) and Eqn (27)
    c = sign(y) * (κ0^3/ w^2) * sqrt((w^2 - p^2) * (1 - w^2)) # Eqn (24)

    # TEST PARAMETERS
    κ_func = s -> kappa(κ0, p, w, s) 
    τ_func = s -> tau(c, λ2, κ0, p, w, s) 

    L = n*4*w* Elliptic.K(p^2)/κ0 # Eqn (46)
    ode_p = (κ_func, τ_func)
    sspan = (0.0, L)

    prob = ODEProblem(frenet_serret, u0, sspan, ode_p)
    sol = solve(prob, saveat=L/300, abstol=1e-10, reltol=1e-10)

    # Putting the soltions to the ODE into a matrix 
    X = [sol.u[i][1] for i in 1:length(sol.u)]
    Y = [sol.u[i][2] for i in 1:length(sol.u)]
    Z = [sol.u[i][3] for i in 1:length(sol.u)]
 
    plot3d!(fig, X, Y, Z)
    display(plot3d(X, Y, Z))


    push!(c_vals,   c)
    push!(λ2_vals,  λ2)
    push!(κ0_vals,  κ0)



    pos_err = norm(sol.u[end][1:3] - sol.u[1][1:3])
    println("x=", round(x,digits=2), " y=", round(y,digits=2), 
            "  L=", round(L,digits=4),
            "  err=", round(pos_err, digits=6))
end


display(fig)



# diagnostics printout 
#for (idx, pt) in enumerate(contour_points)
#    X, Y, p, w = pt
    # recompute pos_err for this point same as in your loop
#    println(idx, "  X=", X, " Y=", Y, " κ0=", κ0_vals[idx]," p=", p, " w=", w, " c=", c_vals[idx], " λ2=", λ2_vals[idx])
#end


# compute c and lambda2 (lambda2 = x/w, c = sign(y) * (κ0^3/w^2)*sqrt((w^2-p^2)(1-w^2))) kappa0 scale of curvature 
# set the arc length  L = 3 * 4wK(p^2)/κ0
# integrate frenet serret 


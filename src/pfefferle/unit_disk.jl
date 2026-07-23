

using Plots, Elliptic, Roots
include("elasticae.jl")
 
# Computes U, V and ΔΦ over the unit disk X²+Y²≤1.
# Reproduces Figures 1, 2, 3 from Pfefferle
 
x_grid = -1.0:0.01:1.0
y_grid = -1.0:0.01:1.0
 
R = zeros(length(x_grid), length(y_grid))
p_vals = zeros(length(x_grid), length(y_grid))
w_vals = zeros(length(x_grid), length(y_grid))
J = zeros(length(x_grid), length(y_grid))
U = zeros(length(x_grid), length(y_grid))
V = zeros(length(x_grid), length(y_grid))
M = zeros(length(x_grid), length(y_grid))
ΔΦ = zeros(length(x_grid), length(y_grid))
 
for (i, x) in enumerate(x_grid)
    for (j, y) in enumerate(y_grid)
        if x^2 + y^2 <= 1
            R[i, j] = sqrt(x^2 + y^2)
            p_vals[i, j] = find_zero(p -> A(p) - R[i,j]^2, (0, 1))
            w_vals[i, j] = sqrt(y^2 + p_vals[i,j]^2)
            w0 = sqrt(1 - w_vals[i,j]^2)
 
            J[i, j] = sqrt((1 - R[i,j]^2)^2 + 4*(y*w0 - x*w_vals[i,j])^2) / (4*w_vals[i,j]^2)
            U[i, j] = (x*(1 - R[i,j]^2) + 2*(y*w0 - x*w_vals[i,j])*w_vals[i,j]) / (4*J[i,j]*w_vals[i,j]^2)
            V[i, j] = (y*(1 - R[i,j]^2) - 2*(y*w0 - x*w_vals[i,j])*w0) / (4*J[i,j]*w_vals[i,j]^2)
 
            M[i, j] = p_vals[i,j]^2 / (p_vals[i,j]^2 + V[i,j]^2)
            N = U[i,j] * ((1 - R[i,j]^2) / (2*(p_vals[i,j]^2 + V[i,j]^2)) - 1) - (4*w_vals[i,j]^2 * J[i,j] * x) / (2*(p_vals[i,j]^2 + V[i,j]^2))
 
            η = acos(sqrt(V[i,j]^2 / (1 - p_vals[i,j]^2)))
            p0 = sqrt(1 - p_vals[i,j]^2)
 
            K_val = Elliptic.K(p_vals[i,j]^2)
            E_complete = Elliptic.E(p_vals[i,j]^2)
            E_inc = Elliptic.E(η, p0^2)
            F_inc = Elliptic.F(η, p0^2)
 
            ΔΦ[i, j] = (1/2) * ((U[i,j] * K_val)/(pi/2) + (sign(N)/(pi/2)) * (K_val*E_inc - (K_val - E_complete)*F_inc))
        end
    end
end
 
display(contour(x_grid, y_grid, U', title="U contours (Fig 1)"))
display(contour(x_grid, y_grid, V', title="V contours (Fig 2)"))
display(contour(x_grid, y_grid, ΔΦ', title="ΔΦ/2πn contours (Fig 3)"))
using LinearAlgebra
include("fourier_curves.jl")

# Physical constant
const μ0 = 4*pi + 1e-7 #permeability of free space

# Core biot-savart law for arbitrary current sources
# re: (n_eval, 3) - evaluation points in cartesian
# rs: (n_src, 3) - source points
# JdV: (n_src, 3) - current element I*dl at each source point
function biot_savart_general(re, rs, JdV)
    re   = reshape(re,   :, 3)
    rs   = reshape(rs,   :, 3)
    JdV  = reshape(JdV,  :, 3)
    
    
    n_eval = size(re, 1)
    B = zeros(n_eval, 3)
    for i in 1:n_eval
        for j in 1:size(rs, 1)
            dr = rs[j,:] - re[i,:] # r_src - r_eval
            num = cross(dr, JdV[j,:]) # dr x J dV
            den = norm(dr)^3
            B[i,:] += cross(dr, JdV[j,:]) * ifelse(den > eps(), 1/den, zero(den)) # ifelse(condition, a (True), b (False))
        end
    end
    return B * μ0 / (4*pi)
end


# Filamentary coil wrapper (similar to DESC's biot_savart_quad)
function biot_savart_quad(eval_pts, coil_pts, tangents, current)
    return biot_savart_general(eval_pts, coil_pts, current*tangents)
end

# Discritizes a Fourier coil (similar to DESC _Coil._compute_A_or_B)
# eval_pts: (n_eval, 3) - evaluation points
# params: coil Fourier coefficients
# current: coil current in Amps
# N_src: numer of quadrature points along the coil
function compute_B_coil(eval_pts, params, current, N_src=100)
    # discretize the coil
    t_vals = range(0, 2*pi, length=N_src+1)[1:end-1]
    dt = 2*pi/N_src

    # compute source points and tangents
    # mirrors DESC self compute
    coil_pts = zeros(N_src, 3)
    tangents = zeros(N_src, 3)
    for (k, t) in enumerate(t_vals)
        # coil position x(t)
        coil_pts[k,:] = [curve_x(t, params), curve_y(t, params), curve_z(t, params)]

        #tangent vector x_s(t) * ds
        tangents[k,:] = [dcurve_x(t, params), dcurve_y(t, params), dcurve_z(t, params)] * dt
    end
    return biot_savart_quad(eval_pts, coil_pts, tangents, current)
end





# ============================================================
# TEST — circular loop of radius R, analytic solution on axis:
# Bz = μ0 * I * R² / (2 * (R² + z²)^(3/2))
# ============================================================
function test_biot_savart()
    println("Testing Biot-Savart on circular loop...")
 
    R = 1.0     # coil radius (m)
    I = 1.0     # current (A)
    z = 0.5     # height above coil center (m)
 
    # unit circle in xy plane using params layout from fourier_curves.jl
    # x(t) = R*cos(t), y(t) = R*sin(t), z(t) = 0
    test_params = zeros(6*NF)
    test_params[1]  = R    # Xc[1] → x = R*cos(t)
    test_params[13] = R    # Ys[1] → y = R*sin(t)
 
    # evaluate B on axis at (0, 0, z)
    eval_pts = [0.0 0.0 z]
    B = compute_B_coil(eval_pts, test_params, I, 1000)
 
    # analytic solution
    Bz_analytic = μ0 * I * R^2 / (2 * (R^2 + z^2)^(3/2))
 
    println("  Computed  Bz = ", B[1, 3])
    println("  Analytic  Bz = ", Bz_analytic)
    println("  Error        = ", abs(B[1, 3] - Bz_analytic) / abs(Bz_analytic) * 100, "%")
end
 
test_biot_savart()
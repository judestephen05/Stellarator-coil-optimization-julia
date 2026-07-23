using LinearAlgebra
include("fourier_curves.jl")
include("biot_savart.jl")


# Compute B*n on the plasma surface (QuadraticFlux.compute())
# coil_params: Fourier coeifficients of the coil
# current: coil current in Amps
# plasma_pts: (n_plasma, 3) - points on plasma surface in cartesian
# plasma_nrho: (n_plasma, 3) - outward normal vectors on plasma surface
function quadratic_flux(coil_params, current, plasma_pts, plasma_nrho, N_src=100)
    # compute B at every plasma surface point from this coil
    B = compute_B_coil(plasma_pts, coil_params, current, N_src)

    # dot with surface normal at each point
    Bn = [dot(B[i, :], plasma_nrho[i, :]) for i in 1:size(plasma_pts, 1)]

    # return sum of squares
    return sum(Bn.^2)/2
end





# ============================================================
# TEST — circular coil of radius R at z=0, flat plasma surface
# at height z=h with normal n=[0,0,1] everywhere.
# On this surface B·n = Bz, which we know analytically at any
# point (r, 0, h) from the full off-axis Biot-Savart formula.
#
# On axis (r=0) the analytic Bz is:
#   Bz = μ0 * I * R² / (2 * (R² + h²)^(3/2))
# ============================================================
function test_quadratic_flux()
    println("Testing quadratic flux on flat surface above circular coil...")
 
    R = 1.0   # coil radius (m)
    I = 1.0   # current (A)
    h = 0.5   # height of plasma surface above coil (m)
 
    # circular coil in xy plane
    test_params = zeros(6*NF)
    test_params[1]  = R    # Xc[1] → x = R*cos(t)
    test_params[13] = R    # Ys[1] → y = R*sin(t)
 
    # plasma surface: a ring of points at radius r=0.2, height z=h
    # normal is [0,0,1] everywhere since surface is flat
    n_pts = 16
    phi_vals = range(0, 2π, length=n_pts+1)[1:end-1]
    r_plasma = 0.2
 
    plasma_pts  = zeros(n_pts, 3)
    plasma_nrho = zeros(n_pts, 3)
    for (k, phi) in enumerate(phi_vals)
        plasma_pts[k, :]  = [r_plasma*cos(phi), r_plasma*sin(phi), h]
        plasma_nrho[k, :] = [0.0, 0.0, 1.0]
    end
 
    # compute quadratic flux
    chi2 = quadratic_flux(test_params, I, plasma_pts, plasma_nrho, 1000)
 
    # verify by computing Bz at each plasma point directly and squaring
    B = compute_B_coil(plasma_pts, test_params, I, 1000)
    Bn_vals = [dot(B[i, :], plasma_nrho[i, :]) for i in 1:n_pts]
    chi2_direct = sum(Bn_vals.^2) / 2
 
    println("  χ²_B from quadratic_flux() = ", chi2)
    println("  χ²_B from direct B·n sum   = ", chi2_direct)
    println("  Match: ", isapprox(chi2, chi2_direct, rtol=1e-10))
 
    # sanity check: on-axis Bz matches analytic
    eval_axis = [0.0 0.0 h]
    B_axis = compute_B_coil(eval_axis, test_params, I, 1000)
    Bz_computed = B_axis[1, 3]
    Bz_analytic = μ0 * I * R^2 / (2 * (R^2 + h^2)^(3/2))
    println("  On-axis Bz computed = ", Bz_computed)
    println("  On-axis Bz analytic = ", Bz_analytic)
    println("  On-axis error       = ", abs(Bz_computed - Bz_analytic)/abs(Bz_analytic) * 100, "%")
end
 
test_quadratic_flux()




using LinearAlgebra

# initial guess — unit circle in xy plane
NF = 4
params = zeros(6*NF)
params[1]  = 1.0    # Xc[1] → x(t) = cos(t)
params[13] = 1.0    # Ys[1] → y(t) = sin(t)
 
curve_x(t, p) = sum(p[n]*cos(n*t)    + p[n+4]*sin(n*t)  for n in 1:NF)
curve_y(t, p) = sum(p[n+8]*cos(n*t)  + p[n+12]*sin(n*t) for n in 1:NF)
curve_z(t, p) = sum(p[n+16]*cos(n*t) + p[n+20]*sin(n*t) for n in 1:NF)
 
dcurve_x(t, p) = sum(n*p[n+4]*cos(n*t)   - n*p[n]*sin(n*t)    for n in 1:NF)
dcurve_y(t, p) = sum(n*p[n+12]*cos(n*t)  - n*p[n+8]*sin(n*t)  for n in 1:NF)
dcurve_z(t, p) = sum(n*p[n+20]*cos(n*t)  - n*p[n+16]*sin(n*t) for n in 1:NF)
 
d2curve_x(t, p) = sum(-n^2*p[n]*cos(n*t)    - n^2*p[n+4]*sin(n*t)  for n in 1:NF)
d2curve_y(t, p) = sum(-n^2*p[n+8]*cos(n*t)  - n^2*p[n+12]*sin(n*t) for n in 1:NF)
d2curve_z(t, p) = sum(-n^2*p[n+16]*cos(n*t) - n^2*p[n+20]*sin(n*t) for n in 1:NF)
 
d3curve_x(t, p) = sum(-n^3*p[n+4]*cos(n*t)  + n^3*p[n]*sin(n*t)    for n in 1:NF)
d3curve_y(t, p) = sum(-n^3*p[n+12]*cos(n*t) + n^3*p[n+8]*sin(n*t)  for n in 1:NF)
d3curve_z(t, p) = sum(-n^3*p[n+20]*cos(n*t) + n^3*p[n+16]*sin(n*t) for n in 1:NF)
 
# arc length element ds/dt = |x'|
ds(t, p) = norm([dcurve_x(t,p), dcurve_y(t,p), dcurve_z(t,p)])
 
# κ² = |x' × x''|² / |x'|⁶
function κ_sq(t, p)
    xp  = [dcurve_x(t,p),  dcurve_y(t,p),  dcurve_z(t,p)]
    xpp = [d2curve_x(t,p), d2curve_y(t,p), d2curve_z(t,p)]
    return norm(cross(xp, xpp))^2 / norm(xp)^6
end
 
# τ = (x' × x'') · x''' / |x' × x''|²
function τ_geom(t, p)
    xp   = [dcurve_x(t,p),  dcurve_y(t,p),  dcurve_z(t,p)]
    xpp  = [d2curve_x(t,p), d2curve_y(t,p), d2curve_z(t,p)]
    xppp = [d3curve_x(t,p), d3curve_y(t,p), d3curve_z(t,p)]
    return dot(cross(xp, xpp), xppp) / norm(cross(xp, xpp))^2
end
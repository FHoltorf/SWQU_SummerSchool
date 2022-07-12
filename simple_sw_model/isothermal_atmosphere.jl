cd(@__DIR__)
using Pkg
Pkg.activate(".")

using DifferentialEquations, CairoMakie, LinearAlgebra

# write out differential equation:
# q = dT/dz
# dP/dz = - n*M*g
# dn/dz = 1/T * (n*q - n*M*g/k)
# dq/dz = -Q_EUV(Z,T,z_break,sharpness)
#


Q_EUV(z,t,z_break,sharpness) = 20*1/2*(1 + tanh((z-z_break)*sharpness))

# Temperature Profile
N = 100
z_range = range(0, 1, length = N)
Δz = step(z_range)

# Laplacian with dirichlet and neumann bc
dl = 1/Δz * [ones(N-1)...]
du = [0, 1/Δz * ones(N-2)...]
d = [1, -2/Δz * ones(N-2)..., -1/Δz]
b = [200.0, map(x -> -Q_EUV.(x, 0, 0.3, 15.0), z_range[2:end-1])...,0.0]
Δ = Tridiagonal(dl, d, du)

T_profile = Δ\b
lines(T_profile)

function T(z,T_profile,Δz) 
    i = min(floor(Int64, z/Δz) + 1, length(T_profile)-1)
    return T_profile[i] + (T_profile[i+1] - T_profile[i])*(z-(i-1)*Δz)/Δz
end

function dTdz(z,T_profile, Δz)
    i = min(floor(Int64, z/Δz) + 1, length(T_profile)-1)
    return (T_profile[i+1] - T_profile[i])/Δz
end

M = 28.97/6.0221408e23*1e-3 # g/mol *  
g = 9.81
k = 1.380649e-23
M̃ = M*g/k
# number density
#dn/dz = 1/T * (n*q - n*M*g/k)

function dndz(n, p, z)
    T, dTdz, M̃, s = p
    return 1/T(z) * (n*dTdz(z) - n*M̃*s)
end

n0 = 1e19
p = (z->T(z,T_profile, Δz), z->dTdz(z,T_profile, Δz), M̃, 700e3)
prob = ODEProblem(dndz, n0, (0.0, 1.0), p)
sol = solve(prob, dt = Δz)

fig = Figure()
ax = Axis(fig[1,1], xscale=log10)
lines!(ax, map(sol, z_range), 600*z_range .+ 200, label = "computed T profile")
fig

p = (z->800*(1-exp(-z/0.05)) + 200, z-> 800*0.05*exp(-z/0.05), M̃, 700e3)
prob = ODEProblem(dndz, n0, (0.0, 1.0), p)
sol = solve(prob, dt =0.01)
lines!(ax, map(sol, z_range), 600*z_range .+ 200, label = "exponential T profile")

p = (z->900.0, z->0.0, M̃, 700e3)
prob = ODEProblem(dndz, n0, (0.0, 1.0), p)
sol = solve(prob, dt =0.01)

lines!(ax, map(sol, z_range), 600*z_range .+ 200, label = "constant T profile")
axislegend(ax, position = :lb)
fig

# solving nonlinear PDE
dl = 1/Δz * [ones(N-1)...]
du = [0, 1/Δz * ones(N-2)...]
d = [1, -2/Δz * ones(N-2)..., -1/Δz]
thermal_conductivity(T) = 3.6e-4*T^0.75
rhs(T) = [200.0, (map(x -> -0.1*Q_EUV(x, 0, 0.3, 15.0), z_range[2:end-1]) ./ thermal_conductivity.(T[2:end-1])) ...,0.0]
Δ = Tridiagonal(dl, d, du)

# fixed point iteration
error = Inf
iter = 0.0
Temp = 800.0 * ones(size(z_range))
while error > 1e-6 && iter < 100
    Temp_new = Δ\rhs(Temp)
    error = norm(Temp_new - Temp)
    Temp .= Temp_new
    iter += 1
end

p = (z->T(z,Temp, Δz), z->dTdz(z,Temp, Δz), M̃, 700e3)
prob = ODEProblem(dndz, n0, (0.0, 1.0), p)
sol = solve(prob, dt = Δz)

lines!(ax, map(sol, z_range), 600*z_range .+ 200, label = "accurate T profile")
axislegend(ax, position = :lb)
fig

# time-dependent problem (temperature)
# ∂T/∂t = Q_EUV/(ρ*cp) + 1/(ρ*cp)*(∂λ/∂z⋅∂T/∂z + λ ∂²T/∂z²)
# ∂n_i/∂z = -1/T * (n∂T/∂z + n*M̃) + s_iᵀ*r(n)
# 

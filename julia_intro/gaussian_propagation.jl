using LinearAlgebra, ForwardDiff, CairoMakie

struct LinearGaussianApproximation
    μ
    Σ
end

const LGA = LinearGaussianApproximation
import Base.+
import Base.-
import Base.*

function *(a::LGA, b::Number)
    return LGA(a.μ*b, b^2*a.Σ)
end
function *(a::Number, b::LGA)
    return b*a
end
function +(a::LGA, b::LGA) 
    return LGA(a.μ + b.μ, a.Σ + b.Σ)
end
function +(a::LGA, b::Number) 
    return LGA(a.μ + b, a.Σ)
end
function +(a::Number, b::LGA) 
    return b+a
end
function -(a::LGA, b::LGA) 
    return LGA(a.μ-b.μ, a.Σ + b.Σ)
end
function -(a::LGA, b::Number) 
    return LGA(a.μ-b, a.Σ)
end
function -(a::Number, b::LGA) 
    return -1*(b - a)
end

k = [1.0, 3.0, 0.4]
S = [-1 0 0;
      1 -1 1
      0 2 -2]

rxns(x,k) = [k[1]*x[1],k[2]*x[2],k[3]*x[3]*x[3]]

f(x,S,k) = S*rxns(x,k)

function f(x::LGA, S, k)
    J = ForwardDiff.jacobian(y->f(y,S,k), x.μ)
    μ = f(x.μ,S,k)
    return LGA(μ, J*x.Σ*J')
end

function reactor_recursion(x, S, k, dt)
    return x + dt*f(x,S,k)
end

function lotka_volterra(x,c)
    dx = [c[1]*x[1] - c[2]*x[1]*x[2];
          -c[3]*x[2] + c[4]*x[1]*x[2]]
    return dx
end
function lotka_volterra(x::LGA, c)
    J = ForwardDiff.jacobian(y->lotka_volterra(y,c), x.μ)
    μ = lotka_volterra(x.μ,c)
    return LGA(μ, J*x.Σ*J')
end

function lotka_volterra_recursion(x,c,dt)
    return x + dt*lotka_volterra(x,c)
end

N = 12000
dt = 0.001
μ_0 = [1.0, 0.25]
V = randn(2,2)
Σ_0 = 0.1*V'*V
xs = [LGA(μ_0,Σ_0)]
c = [1,2,1,2]
for i in 1:N
    push!(xs, lotka_volterra_recursion(xs[end], c, dt))
end

fig = Figure();
ax = Axis(fig[1,1])
lines!(ax,[x.μ[1] for x in xs])
lines!(ax,[x.μ[2] for x in xs])
display(fig)


fig = Figure();
ax = Axis(fig[1,1])
limits!(ax, 0, 2, 0, 2)
probing = -5:0.05:5
probing_grid = [[p1,p2] for p1 in probing, p2 in probing]
mean = Observable(xs[1].μ)
std = Observable(xs[1].Σ)
gaussian = Observable([p'*std.val*p for p in probing_grid])
shifted_x = Observable(mean.val[1] .+ probing)
shifted_y = Observable(mean.val[2] .+ probing)

contour!(ax, shifted_x, shifted_y, gaussian, levels = 0.002:0.002:0.02)

record(fig, string(@__DIR__, "/gaussian_rxn_network.mp4"), 1:10:length(xs), framerate=24) do i
    mean[] = xs[i].μ
    std[] = xs[i].Σ
    gaussian[] = [p'*std.val*p for p in probing_grid]
    shifted_x[] = mean.val[1] .+ probing
    shifted_y[] = mean.val[2] .+ probing
end

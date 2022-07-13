using DifferentialEquations, Measurements
import Base: +, -, *, /


struct GaussApprox
    mean
    var
end

const GA = GaussApprox

function +(a::GA, b::GA)
    return GA(a.mean+b.mean, a.var+b.var)
end

function +(a::Number, b::GA) 
    return GA(a + b.mean, b.var)
end

+(a::GA, b::Number) = b+a

function *(a::Number, b::GA)
    return GA(a*b.mean, a^2*b.var)
end

*(a::GA, b::Number) = b*a

function *(a::GA, b::GA)
    return GA(a.mean*b.mean, b.mean^2*a.var + a.mean^2*b.var)
end

function -(a::GA, b::GA)
    return GA(a.mean-b.mean, a.var + b.var)
end
-(a::GA) = GA(-a.mean, a.var)
-(a::GA, b::Number) = a + (-b)
-(a::Number, b::GA) = a + (-b)

/(a::GA, b::Number) = GA(a.mean/b, a.var/b^2)

function lotka_volterra(x,c)
    dx = [c[1]*x[1] - c[2]*x[1]*x[2];
          -c[3]*x[2] + c[4]*x[1]*x[2]]
    return dx
end

function lotka_volterra_recursion(x,c,dt)
    return x + dt*lotka_volterra(x,c)
end

N = 12000
dt = 0.001
x0 = [GA(1.0, 0.1), GA(0.25, 0.05)]
xs = [x0]
c = [GA(1.0, 0.1),GA(2, 0.3),1,2]
for i in 1:N
    push!(xs, lotka_volterra_recursion(xs[end], c, dt))
end

fig = Figure()
ax = Axis(fig[1,1])
lines!(ax, [x[1].var for x in xs])

function simple_pendulum(x,(g,L),t)
    dx = [x[2], -(g/L)*x[1]]
    return dx
end 

function simple_pendulum_recursion(x,(g,L),dt)
    return x + dt*simple_pendulum(x,(g,L),0)
end

N = 12000
dt = 0.001
x0 = [GA(π/60, 0.01), GA(0.0, 0.01)]
xs = [x0]
g = GA(9.81, 0.02)
L = 1.0
for i in 1:N
    push!(xs, simple_pendulum_recursion(xs[end], (g,L), dt))
end

fig = Figure()
ax = Axis(fig[1,1])
lines!(ax, [x[1].mean for x in xs])
fig


N = 12000
dt = 0.001
x0 = [π/60 ± 0.01, 0.0 ± 0]
xs = [x0]
g = 9.81 ± 0.02
L = 1.0
for i in 1:N
    push!(xs, simple_pendulum_recursion(xs[end], (g,L), dt))
end

fig = Figure()
ax = Axis(fig[1,1])
lines!(ax, [x[1].mean for x in xs])
fig

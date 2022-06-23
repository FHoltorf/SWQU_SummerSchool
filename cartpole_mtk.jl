using DifferentialEquations, ModelingToolkit, GLMakie

const g = 9.81
@variables t, θ(t), x(t), v(t), ω(t)
@parameters m, M, L, F
d = Differential(t)

rhs = 1/(cos(θ)^2 * M  - L*(m+M)) * [M*L*cos(θ) -L;
                                        -(m+M) cos(θ)] * [-g*sin(θ);
                                                          F + M * L * ω^2 * sin(θ)]

dynamics = [d(x) ~ v;
            d(θ) ~ ω;
            d(v) ~ rhs[1];
            d(ω) ~ rhs[2]]

cart_pole = ODESystem(dynamics, name = :cart_pole)

u0 = Dict(x => 0.0, θ => π-0.2, v => 0.0, ω => 0.0)
ps = Dict(m => 0.1, M => 1.0, L => 1.0, F => 0.0)
cart_pole_problem = ODEProblem(cart_pole, u0, (0.0, 10.0), ps)
cart_pole_solution = solve(cart_pole_problem)

# visualization
fig = Figure();
ax = Axis(fig[1,1]);
xlims!(ax, (-1.5*ps[L],1.5*ps[L]))
ylims!(ax, (-1.5*ps[L],1.5*ps[L]))
pendulum = Observable(Point2f(u0[x] + sin(u0[θ])*ps[L], -cos(u0[θ])*ps[L]))
cart = Observable(Rect(u0[x]-0.2, -0.1, 0.4, 0.2))
rod = Observable([pendulum.val, Point2f(u0[x],0.0)])
scatter!(ax, pendulum, color = :black, markersize = 20)
poly!(ax, cart, color = :red)
lines!(ax, rod, color = :black, linewidth = 5)
display(fig)

record(fig, string(@__DIR__,"/uncontrolled_cart.mp4"), 0:0.01:10.0, framerate = 24) do t
    u = cart_pole_solution(t)
    pendulum[] = Point2f(u[1] + ps[L]*sin(u[2]), -ps[L]*cos(u[2]))
    cart[] = Rect(u[1]-0.2, -0.1, 0.4, 0.2)
    rod[] = [pendulum.val, Point2f(u[1],0.0)]
end





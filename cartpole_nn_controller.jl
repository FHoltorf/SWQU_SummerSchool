cd(@__DIR__)
using Pkg; Pkg.activate("."); Pkg.instantiate()

using Flux, DiffEqFlux, DiffEqSensitivity
using DifferentialEquations
using Statistics
using Zygote, ForwardDiff
using LinearAlgebra
using ControlSystems
using Random
using GLMakie, ColorSchemes
using FileIO

const g = 9.81
const m = 0.1
const M = 1.0
const L = 1.0

T = 5.0
n_T = 100.0
Δt = T/(n_T-1) 

n_inits = 10
nn = FastChain(FastDense(4, 10, tanh, initW = Flux.glorot_uniform, initb = Flux.glorot_uniform),
                FastDense(10, 10, tanh, initW = Flux.glorot_uniform, initb = Flux.glorot_uniform),
                FastDense(10, 1, x-> 10*tanh(x), initW = Flux.glorot_uniform, initb = Flux.glorot_uniform))

p_nn = initial_params(nn)

θ0_init = range(2*π*1/n_inits,2π*(1-1/n_inits), length = n_inits)
u0_init = [[0.0, θ, 0.0, 0.0] for θ in θ0_init]
ref_point = [0.0, π, 0.0, 0.0] 
Q = Diagonal([1.0, 10.0, 0.001, 0.001])
R = reshape([0.1],1,1)

function cart_pole_dynamics(u,F)
    return [u[3];
            u[4];
            1/(cos(u[2])^2 * L * M  - L*(m+M)) * (-g*sin(u[2])*M*L*cos(u[2])-L*(F + M * L * u[4]^2 * sin(u[2])));
            1/(cos(u[2])^2 * L * M  - L*(m+M)) * (g*sin(u[2])*(m+M) + cos(u[2])*(F + M * L * u[4]^2 * sin(u[2])))]
end

function tracking_error(x,ref_point)
    return [x[1] - ref_point[1], cos(x[2]) - cos(ref_point[2]), x[3] - ref_point[3], x[4] - ref_point[4]]
end

J = ForwardDiff.jacobian(x->cart_pole_dynamics(x[1:4],x[5]), vcat(ref_point, 0))
A = J[1:4, 1:4]
B = J[:,5]

# LQR 
K_lqr = lqr(Continuous, A, B, Q, R)

function soft_step(x, lb, ub, alpha=10)
    return 1/4*(1+tanh((ub-x)*alpha))*(1+tanh((x-lb)*alpha))
end

function cart_pole!(du, u, (K, ref_point, p_nn), t)
    α = soft_step(cos(ref_point[2]) - cos(u[2]), -0.25, 0.25, 20) #exp(-5*(cos(ref_point[2]) - cos(u[2]))^2)
    #error = tracking_error(u, ref_point)
    F = - ((K * (u-ref_point))[1] * α + (1-α) * (nn(u-ref_point, p_nn) - nn(ref_point-u, p_nn))[1])*soft_step(u[1], -10.0, 10.0)
    du[1] = u[3]
    du[2] = u[4]
    du[3] = 1/(cos(u[2])^2 * L * M  - L*(m+M)) * (-g*sin(u[2])*M*L*cos(u[2])-L*(F + M * L * u[4]^2 * sin(u[2])))
    du[4] = 1/(cos(u[2])^2 * L * M  - L*(m+M)) * (g*sin(u[2])*(m+M) + cos(u[2])*(F + M * L * u[4]^2 * sin(u[2])))
    return nothing
end

# spin condition
condition(u,t,integrator) = (abs(u[2]) > 3*2*π)
function affect!(integrator)
    correction = floor(abs(integrator.u[2])/(2π))
    integrator.u[2] -= sign(integrator.u[2])*correction*2π 
end
cb = DiscreteCallback(condition,affect!,save_positions=(false,false))

function loss(p, u0_init, T, Δt, ref_point, cb;
            K = zeros(1,4),
            obj_weights = ones(4), 
            sensealg = ForwardDiffSensitivity())
    cart_pole_prob = ODEProblem((du,u,p,t) -> cart_pole!(du, u, (K,ref_point,p), t), u0_init[1], (0.0, T), p, callback=cb)
    scaled_weights = sqrt.(Δt * obj_weights)
    cum_error = 0.0
    for (i,u0) in enumerate(u0_init)
        _prob = remake(cart_pole_prob, u0=u0)
        sol_object = solve(_prob, sensealg=sensealg, saveat=Δt)
        #if sol_object.retcode == :Success
        sol = Array(sol_object)
        cum_error += sum(abs2, scaled_weights[2] .* (cos.(sol[2,:]) .- cos(ref_point[2])))
        cum_error += sum(abs2, scaled_weights[[1,3,4]] .* (sol[[1,3,4],:] .- ref_point[[1,3,4]]))
        #else # for linesearch
        #    println("initial condition $i was instable")
        #    tracking_error += sum(abs2, scaled_weights[2] .* (cos.(sol[2,:]) .- cos(ref_point[2])))
        #    tracking_error += sum(abs2, scaled_weights[[1,3,4]] .* (sol[[1,3,4],:] .- ref_point[[1,3,4]]))
        #    tracking_error += 1e4
        #end
    end
    return cum_error
end

sensealg = ForwardDiffSensitivity()
#=
# differentiation tests -> why is reverse diff so bad
dp = @time Zygote.gradient(p -> loss(p, u0_init, T, Δt, ref_point, sensealg=sensealg), p_nn)
l, back = @time Zygote.pullback(p -> loss(p, u0_init, T, Δt, ref_point, sensealg=sensealg), p_nn)
dp = @time back(l)[1]
=#
function linesearch!(p, dp, f, f0, α = 1.0; β = 0.5, iter_max = 100)
    f1 = Inf
    s = α
    iter = 0
    while f1 >= f0 && iter < iter_max
        f1 = f(p - s*dp)
        s = s*β
        iter += 1
    end
    p .-= (s/β)*dp
end
#optimizer = ADAM(stepsize) # richtiger schmutz
l_opt = Inf
p_opt = deepcopy.(p_nn)

losses = []
obj(p) = loss(p, u0_init, T, Δt, ref_point, cb, K = K_lqr, sensealg=sensealg, obj_weights = 1/10*diag(Q))
obj(p_opt)
stepsize = 1.0
steps = 200
ls_iter_max = 10
p_nn .= deepcopy(p_opt)
for iter in 1:steps
    l, back = @time Zygote.pullback(obj, p_nn)
    if l < l_opt
        l_opt = l
        p_opt .= deepcopy(p_nn)
    end
    dp = @time back(one(l))[1]
    println("iter: $iter - loss = $l - ||dp|| = $(norm(dp))")
    push!(losses, l)
    linesearch!(p_nn, dp/norm(dp), obj, l, stepsize, iter_max = ls_iter_max)
    println("")
end

opacity = 0.2
T_test = 10.0
cart_pole_solution = Dict()
u0_init_test = [[0.0, 2π*rand(), 0.0, 0.0] for i in 1:10]
for u0 in u0_init_test
    prob = ODEProblem(cart_pole!, u0, (0.0, T_test), [K_lqr, ref_point, p_opt])
    cart_pole_solution[u0] = solve(prob)
end

fig = Figure();
ax = Axis(fig[1,1], aspect = 5.0/1.5);
xlims!(ax, (-5.0*L,5.0*L))
ylims!(ax, (-1.5*L,1.5*L))
pendulum = Dict(u0 => Observable(Point2f(u0[1] + sin(u0[2])*L, -cos(u0[2])*L)) for u0 in u0_init_test)
cart = Dict(u0 => Observable(Rect(u0[1]-0.2, -0.1, 0.4, 0.2)) for u0 in u0_init_test)
rod = Dict(u0 => Observable([pendulum[u0].val, Point2f(u0[1],0.0)]) for u0 in u0_init_test)
for u0 in u0_init_test
    scatter!(ax, pendulum[u0], color = (:black,opacity), markersize = 20)
    poly!(ax, cart[u0], color = (:red,opacity))
    lines!(ax, rod[u0], color = (:black,opacity), linewidth = 5)
end
display(fig)

record(fig, string(@__DIR__,"/nn_controlled_cart.mp4"), 0:0.01:T_test, framerate = 48) do t
    for u0 in u0_init_test
        u = cart_pole_solution[u0](t)
        pendulum[u0][] = Point2f(u[1] + L*sin(u[2]), -L*cos(u[2]))
        cart[u0][] = Rect(u[1]-0.2, -0.1, 0.4, 0.2)
        rod[u0][] = [pendulum[u0].val, Point2f(u[1],0.0)]
    end
end


fig = Figure()
N = 100
opacity = 0.2
T_test = 5.0
cart_pole_solution = Dict()
u0_init_test = [[0.0, θ, 0.0, 0.0] for θ in range(2π/N, 2π*(1-1/N), length = N)]
for u0 in u0_init_test
    prob = ODEProblem(cart_pole!, u0, (0.0, T_test), [K_lqr, ref_point, p_opt])
    cart_pole_solution[u0] = solve(prob)
end

ax = Axis(fig[1,1], aspect = 5.0/1.5);
xlims!(ax, (-5.0*L,5.0*L))
ylims!(ax, (-1.5*L,1.5*L))
pendulum = Dict(u0 => Observable(Point2f(u0[1] + sin(u0[2])*L, -cos(u0[2])*L)) for u0 in u0_init_test)
trajectory = Dict(u0 => Observable([Point2f(u0[1] + sin(u0[2])*L, -cos(u0[2])*L)]) for u0 in u0_init_test)
for (i,u0) in enumerate(u0_init_test)
    scatter!(ax, pendulum[u0], color = get(colorschemes[:hsv], i/N), markersize = 10)
    lines!(ax, trajectory[u0], color = get(colorschemes[:hsv], i/N),  linewidth = 2)
end
display(fig)

record(fig, string(@__DIR__,"/nn_controlled_cart_trajectories.mp4"), 0:0.01:T_test, framerate = 48) do t
    for u0 in u0_init_test
        u = cart_pole_solution[u0](t)
        pendulum[u0][] = Point2f(u[1] + L*sin(u[2]), -L*cos(u[2]))
        push!(trajectory[u0].val, Point2f(u[1] + L*sin(u[2]), -L*cos(u[2])))
        trajectory[u0][] = trajectory[u0][]
    end
end

save("nn_controller_params.jld", "params", p_nn)

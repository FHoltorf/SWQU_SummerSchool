cd(@__DIR__)
using Pkg; Pkg.activate("."); Pkg.instantiate()

# load packages
using Flux, DiffEqFlux, DiffEqSensitivity
using DifferentialEquations
using Statistics
using Zygote, ForwardDiff
using LinearAlgebra
using Random

const g = 9.81
const m = 0.1
const M = 1.0
const L = 1.0

#################################################
# read parameters from command line
lr = 0.0015f0 #parse(Float64,ARGS[1]) #0.05
epochs = 400 #parse(Int,ARGS[2]) #100
slurmidx = 1 #parse(Int,ARGS[3]) #1

numtraj = 10 # number of trajectories in parallel simulations for training
numtrajplot = numtraj # .. for plotting

# time range for the solver
dt = 0.01f0 #0.001f0
tinterval = 0.02f0
tstart = 0.0f0
Nintervals = 500 # total number of intervals, total time = t_interval*Nintervals
tspan = (tstart,tinterval*Nintervals)

# Hamiltonian parameters
Δ = 20.0f0
Ωmax = 10.0f0 # control parameter (maximum amplitude)

# loss hyperparameters
C1 = Float32(1.0)  # evolution state fidelity
C2 = Float32(0.0) # action amplitudes
C3 = Float32(0.0) # evolution state fidelity for last few steps!

struct Parameters{flType,intType,tType}
	lr::flType
	epochs::intType
	numtraj::intType
	numtrajplot::intType
	dt::flType
	tinterval::flType
	tspan::tType
	Nintervals::intType
	Δ::flType
	Ωmax::flType
	C1::flType
	C2::flType
	C3::flType
end

myparameters = Parameters{typeof(dt),typeof(numtraj), typeof(tspan)}(
  lr, epochs, numtraj, numtrajplot, dt, tinterval, tspan, Nintervals,
  Δ, Ωmax, C1, C2, C3)


################################################
# Define Neural Network

# state-aware
nn = FastChain(
	   FastDense(4, 10, relu, initW = Flux.glorot_uniform, initb = Flux.glorot_uniform),
	   FastDense(10, 10, relu, initW = Flux.glorot_uniform, initb = Flux.glorot_uniform),
#	   FastDense(10, 10, relu, initW = Flux.glorot_uniform, initb = Flux.glorot_uniform),
	   FastDense(10, 1, softsign, initW = Flux.glorot_uniform, initb = Flux.glorot_uniform))
p_nn = initial_params(nn)
# nn = Chain(
# 	Dense(4, 4, relu, initW = Flux.glorot_uniform, initb = Flux.glorot_uniform),
# 	#Dense(64, 64, relu, initW = Flux.glorot_uniform, initb = Flux.glorot_uniform),
# 	#Dense(64, 32, relu, initW = Flux.glorot_uniform, initb = Flux.glorot_uniform),
# 	Dense(4, 1, softsign, initW = Flux.glorot_uniform, initb = Flux.glorot_uniform))
# p_nn, re = Flux.destructure(nn)

###############################################
# initial state anywhere on the Bloch sphere
function prepare_initial(n_par)
  # shape 4 x n_par
  # input number of parallel realizations and dt for type inference
  # random position on the Bloch sphere
  θ = 2π*rand(n_par)  # uniform sampling for cos(theta) between -1 and 1
  # real and imaginary parts ceR, cdR, ceI, cdI
  u0 = [zeros(n_par), θ, zeros(n_par), zeros(n_par)]
  return vcat(transpose.(u0)...) # build matrix
end

# target state
# ψtar = |up>

u0 = prepare_initial(myparameters.numtraj)

# var(u0[1,:].^2+u0[3,:].^2)

###############################################
# Define cart pole

function qubit!(du,u,p,t)
  Ω = nn(u, p)[1]
  @inbounds begin
    du[1] = u[3]
    du[2] = u[4]
    du[3] = 1/(cos(u[2])^2 * M  - L*(m+M)) * (-g*sin(u[2])*M*L*cos(u[2])-L*( Ω + M * L * u[4]^2 * sin(u[2])))
    du[4] = 1/(cos(u[2])^2 * M  - L*(m+M)) * (g*sin(u[2])*(m+M) + cos(u[2])*(Ω + M * L * u[4]^2 * sin(u[2])))
  end
  return nothing
end

  

# define SDE problem
prob = ODEProblem(qubit!, vec(u0[:,1]), myparameters.tspan, p_nn)

sol = @time solve(prob, Tsit5(),
                  dt=myparameters.dt,
                  saveat=myparameters.tinterval,
                  dtmax=myparameters.dt*2,
                  adaptive=true, abstol=1e-6, reltol=1e-6)

#########################################
# compute loss
function obj(u)
  x = @view u[1,:,:]
  θ = @view u[2,:,:]
  v = @view u[3,:,:]
  ω = @view u[4,:,:]
  return mean(x.^2) + mean((θ .- π).^2) + mean(v.^2) + mean(ω.^2) 
end

function loss(p, u0, myparameters::Parameters; sensealg = ForwardDiffSensitivity())
  pars = p

  function prob_func(prob, i, repeat)
    # prepare initial state and applied control pulse
	u0tmp = deepcopy(vec(u0[:,i]))
    remake(prob,
	    p = pars,
	  	u0 = u0tmp
		)
  end

  ensembleprob = EnsembleProblem(prob,
   prob_func = prob_func,
   safetycopy = true
   )

  _sol = solve(ensembleprob, Tsit5(), EnsembleThreads(),
  	sensealg=sensealg,
	saveat=myparameters.tinterval,
	dt=myparameters.dt,
	adaptive=true, abstol=1e-6, reltol=1e-6,
	trajectories=myparameters.numtraj)#, batch_size=myparameters.numtraj)

  A = convert(Array,_sol)
  loss = obj(A)

  return loss
end

function test_loss(p, u0, myparameters::Parameters)
    pars = p

  function prob_func(prob, i, repeat)
    # prepare initial state and applied control pulse
	u0tmp = deepcopy(vec(u0[:,i]))
    remake(prob,
	    p = pars,
	  	u0 = u0tmp
		)
  end

  ensembleprob = EnsembleProblem(prob,
   prob_func = prob_func,
   safetycopy = true
   )

  _sol = solve(ensembleprob, Tsit5(), EnsembleThreads(),
	saveat=myparameters.tinterval,
	dt=myparameters.dt,
	adaptive=true, abstol=1e-6, reltol=1e-6,
	trajectories=myparameters.numtraj, batch_size=myparameters.numtraj)

  A = convert(Array,_sol)
  loss = obj(A)

  return loss
end

opt = ADAM(myparameters.lr)
losses = []
for epoch in 1:myparameters.epochs
  println("epoch: $epoch / $(myparameters.epochs)")
  local u0 = prepare_initial(myparameters.numtraj)
  _dy, back = @time Zygote.pullback(p -> loss(p, u0, myparameters)#, sensealg=InterpolatingAdjoint())
                                          , p_nn)
  gs = @time back(one(_dy))[1]
  push!(losses, _dy)
  if epoch % myparameters.epochs == 0
    local u0 = prepare_initial(myparameters.numtrajplot)
	  println("Loss (epoch: $epoch): $(test_loss(myparameters.p_nn, u0, myparameters))")
  end
  Flux.Optimise.update!(opt, p_nn, gs)
  println("")
end

u0 = vec(prepare_initial(myparameters.numtrajplot)[:,1])
prob = ODEProblem(qubit!, u0, myparameters.tspan, p_nn)
cart_pole_solution = solve(prob)

using GLMakie
fig = Figure();
ax = Axis(fig[1,1]);
xlims!(ax, (-1.5*L,1.5*L))
ylims!(ax, (-1.5*L,1.5*L))
pendulum = Observable(Point2f(u0[1] + sin(u0[2])*L, -cos(u0[2])*L))
cart = Observable(Rect(u0[1]-0.2, -0.1, 0.4, 0.2))
rod = Observable([pendulum.val, Point2f(u0[1],0.0)])
scatter!(ax, pendulum, color = :black, markersize = 20)
poly!(ax, cart, color = :red)
lines!(ax, rod, color = :black, linewidth = 5)
display(fig)

record(fig, string(@__DIR__,"/uncontrolled_cart.mp4"), 0:0.01:10.0, framerate = 24) do t
    u = cart_pole_solution(t)
    pendulum[] = Point2f(u[1] + L*sin(u[2]), -L*cos(u[2]))
    cart[] = Rect(u[1]-0.2, -0.1, 0.4, 0.2)
    rod[] = [pendulum.val, Point2f(u[1],0.0)]
end

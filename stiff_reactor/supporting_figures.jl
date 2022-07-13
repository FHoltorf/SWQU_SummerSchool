using DifferentialEquations, CairoMakie, 
      LaTeXStrings, LinearAlgebra

# stiffness image
const S = [-1 0 0;
           1 -1 1;
           0 2 -2;]

kinetics(x,p) = p .* [x[1],x[2],x[3]^2]
function reactor(x, p, t)
    return S*kinetics(x,p)
end

prob = ODEProblem(reactor, [1.0,0,0], (0.0,5.0), [100.0, 0.25, 1.0])
sol = solve(prob, saveat=0.001)

fig = Figure(fontsize=24);
ax = Axis(fig[1,1], ylabel = "concentration", xlabel = "time")
ax_mag = Axis(fig[1,2], ylabel = "concentration", xlabel = "time")
colors = [:red, :dodgerblue, :black]
labels = ["A", "B", "C"]
for i in 1:3
    lines!(ax, 0:0.001:5.0, [u[i] for u in sol.u], linewidth = 2, color = colors[i], label = labels[i])
    lines!(ax_mag, 0:0.001:0.1, [u[i] for u in sol.u[1:101]], linewidth = 2, color = colors[i], label = labels[i])
end
axislegend(ax_mag, position = :rc)
save(string(@__DIR__, "/assets/stiffness.png"), fig)
save(string(@__DIR__, "/assets/stiffness.pdf"), fig)

# region of stability: Euler
fig = Figure(fontsize=24)
ax = Axis(fig[1,1], ylabel = latexstring("Im(h \\lambda)"), xlabel = latexstring("Re(h \\lambda)"))
ylims!(ax,(-3,3))
xlims!(ax,(-3,3))
discrete_pts = range(-2,0,length=1000)
band!(ax, discrete_pts, zeros(length(discrete_pts)), [sqrt(1-(1+d)^2) for d in discrete_pts],
          color = (:dodgerblue, 0.5))
band!(ax, discrete_pts, zeros(length(discrete_pts)), [-sqrt(1-(1+d)^2) for d in discrete_pts],
          color = (:dodgerblue, 0.5), label = "Region of stability")
hlines!(ax, 0, color = :black)
vlines!(ax, 0, color = :black)
axislegend(ax, position=:rt)
fig

# region of stability: Explicit Euler
fig = Figure(fontsize=24)
ax = Axis(fig[1,1],
          ylabel = latexstring("Im(h \\lambda)"), 
          xlabel = latexstring("Re(h \\lambda)"),
          aspect = 1.0)
function explicit_Euler_RoS!(ax)
    ylims!(ax,(-3,3))
    xlims!(ax,(-3,3))
    discrete_pts = range(-2,0,length=1000)
    band!(ax, discrete_pts, zeros(length(discrete_pts)), [sqrt(1-(1+d)^2) for d in discrete_pts],
            color = (:dodgerblue, 0.5))
    band!(ax, discrete_pts, zeros(length(discrete_pts)), [-sqrt(1-(1+d)^2) for d in discrete_pts],
            color = (:dodgerblue, 0.5), label = "Region of stability")
    hlines!(ax, 0, color = :black)
    vlines!(ax, 0, color = :black)
end

explicit_Euler_RoS!(ax)
save(string(@__DIR__,"/assets/explicit_euler_ros.pdf"), fig)
save(string(@__DIR__,"/assets/explicit_euler_ros.png"), fig)

# region of stability: Implicit Euler
fig = Figure(fontsize=24)
ax = Axis(fig[1,1],
          ylabel = latexstring("Im(h \\lambda)"), 
          xlabel = latexstring("Re(h \\lambda)"),
          aspect = 1.0)

function implicit_Euler_RoS!(ax)
    ylims!(ax,(-3,3))
    xlims!(ax,(-3,3))
    discrete_pts = range(-3,0,length=2)
    band!(ax, discrete_pts,-3*ones(length(discrete_pts)), 3*ones(length(discrete_pts)),
            color = (:dodgerblue, 0.5))
    discrete_pts = range(2,3,length=2)
    band!(ax, discrete_pts,-3*ones(length(discrete_pts)), 3*ones(length(discrete_pts)),
            color = (:dodgerblue, 0.5))
    discrete_pts = range(0,2,length=1000)
    band!(ax, discrete_pts, [sqrt(1-(d-1)^2) for d in discrete_pts], 3*ones(length(discrete_pts)),
            color = (:dodgerblue, 0.5))
    band!(ax, discrete_pts, -3*ones(length(discrete_pts)), [-sqrt(1-(d-1)^2) for d in discrete_pts],
            color = (:dodgerblue, 0.5), label = "Region of stability")
    hlines!(ax, 0, color = :black)
    vlines!(ax, 0, color = :black)
end

implicit_Euler_RoS!(ax)
save(string(@__DIR__,"/assets/implicit_euler_ros.pdf"), fig)
save(string(@__DIR__,"/assets/implicit_euler_ros.png"), fig)

# correct region of stability
fig = Figure(fontsize=24)
ax = Axis(fig[1,1],
          ylabel = latexstring("Im(h \\lambda)"), 
          xlabel = latexstring("Re(h \\lambda)"),
          aspect = 1.0)
function correct_RoS!(ax)
    ylims!(ax,(-3,3))
    xlims!(ax,(-3,3))
    discrete_pts = range(-3,0,length=2)
    band!(ax, discrete_pts, -3*ones(length(discrete_pts)), 3*ones(length(discrete_pts)),
              color = (:dodgerblue, 0.5))
    hlines!(ax, 0, color = :black)
    vlines!(ax, 0, color = :black)
end
correct_RoS!(ax)
save(string(@__DIR__,"/assets/correct_ros.pdf"), fig)
save(string(@__DIR__,"/assets/correct_ros.png"), fig)

fig = Figure(fontsize=24)
names = ["correct", "explicit Euler", "implicit Euler"]
axs = [Axis(fig[1,i],
          ylabel = latexstring("Im(h \\lambda)"), 
          xlabel = latexstring("Re(h \\lambda)"),
          aspect = 1.0,
          title = names[i]) for i in 1:3]
correct_RoS!(axs[1])
explicit_Euler_RoS!(axs[2])
implicit_Euler_RoS!(axs[3])
save(string(@__DIR__,"/assets/ros.pdf"), fig)
save(string(@__DIR__,"/assets/ros.png"), fig)
fig


# flow field
x = y = -2:0.005:2
f(z) = 1 / (z * (z^2 - z - 1 - 3im))
fvals = [f(u + 1im * v) for u in x, v in y]
fvalues = abs.(fvals)
fargs = angle.(fvals)
polya(x, y) = Point2f(real(f(x + 1im * y)), -imag(f(x + 1im * y)))

fig = Figure(resolution = (900, 400))
ax = Axis(fig[1, 1], aspect = 1)
streamplot!(ax, polya, -2 .. 2, -2 .. 2, colormap = [:black, :black],
        gridsize = (40, 40), arrow_size = 6, linewidth = 1)
limits!(ax, -2, 2, -2, 2)

colsize!(fig.layout, 1, Aspect(1, 1.0))
display(fig)
save(string(@__DIR__,"/assets/polya.png"), fig)
save(string(@__DIR__,"/assets/polya.pdf"), fig)
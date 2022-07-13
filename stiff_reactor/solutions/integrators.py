from scipy.optimize import fsolve
import numpy as np
import matplotlib.pyplot as plt
import time

def reaction_rates(x,k):
    return k*x

def reactor(x,t,k,S):
    return S @ reaction_rates(x,k)

def explicit_euler_step(x0,t,f,dt):
    return x0 + dt * f(x0, t)

def implicit_euler_step(x0,t,f,dt):
    return fsolve(lambda x : x0+dt*f(x,t+dt) - x, x0)

def crank_nicolson_step(x0,t,f,dt):
    e_step = dt/2*f(x0,t)
    x_ee = x0 + e_step
    update = lambda x : x_ee + dt/2*f(x,t+dt) - x
    x = fsolve(update, x0)
    return x

def heun_step(x0,t,f,dt):
    e_step = f(x0,t)
    x_ee = x0 + dt*e_step
    return x0 + dt/2*(e_step + f(x_ee,t+dt))

def integrate(f, x0, tspan, dt, int):
    t, tf = tspan
    x = x0
    trajectory = [x0]
    ts = [t]
    while t < tf:
        dt_eff = min(dt, tf-t)
        x = int(x,t,f,dt_eff)
        t = min(t+dt_eff, tf)
        trajectory.append(x)
        ts.append(t)
    return trajectory, ts

# stoichiometry matrix 
S = np.array([[-1, 0, 0],
              [1, -1, 1],
              [0, 2, -2]])

k = np.array([100.0, 0.25, 1.0])

x0 = np.array([1.0, 0.0, 0.0])

dt_range = 10**np.linspace(-4,1,6)
integrators = [implicit_euler_step, explicit_euler_step, crank_nicolson_step, heun_step] 
performance = dict()
for dt in dt_range:
    for integrator in integrators:
        t0 = time.time()
        sol = integrate(lambda x,t : reactor(x,t,k,S), x0, (0.0,10.0), dt, integrator)
        tf = time.time()
        performance[dt, integrator, "CPU_time"] = tf-t0
        performance[dt, integrator, "sol"] = sol 
   
# compare stability and computational performance
def plot_traces(ax, t_range, sol):
    colors = ["red", "dodgerblue", "black"]
    labels = ["A", "B", "C"]
    for i in range(3):
        ax.plot(t_range, [s[i] for s in sol], color = colors[i], linewidth = 2, label = labels[i])
    ax.set_xlabel("t")
    ax.set_ylabel("C(t)")
    plt.tight_layout()
        
integrator2name = {crank_nicolson_step : "Crank-Nicolson",
                   explicit_euler_step : "Explicit Euler",
                   implicit_euler_step : "Implicit Euler",
                   heun_step : "Heun"}

for integrator in integrators:
    fig, axs = plt.subplots(3,2)
    fig.suptitle(integrator2name[integrator])
    for i in range(3):
        for j in range(2):
            dt = dt_range[i*2+j]
            print(dt)
            sol, ts = performance[dt, integrator, "sol"]
            axs[i,j].set_title("dt = "+str(dt))
            plot_traces(axs[i,j], ts, sol)   
    fig.show()
    fig.savefig(integrator2name[integrator]+"_traces.pdf")

def compute_error(sol, ref_sol):
    return sum(np.linalg.norm(sol[i]-ref_sol[i]) for i in range(len(ref_sol)))/len(sol)

benchmark = np.exp(-k[0]*0.01)*x0[0]
fig, axs = plt.subplots(2)
for i in range(2):
    axs[i].set_yscale('log')
    axs[i].set_xscale('log')
    
colors  = {integrators[0] : "red",
           integrators[1] : "dodgerblue",
           integrators[2] : "green",
           integrators[3] : "black"}

axs[1].set_xlabel('h')
axs[0].set_ylabel('error')
axs[1].set_ylabel('CPU time [s]')
for integrator in integrators:
    times = [performance[dt, integrator, "CPU_time"] for dt in dt_range]
    accuracy = [np.linalg.norm(performance[dt, integrator, "sol"][0][round(0.01/dt)][0] - benchmark) for dt in dt_range]
    axs[0].plot(dt_range, accuracy, 'o-', color = colors[integrator], linewidth = 2, label = integrator2name[integrator])
    axs[1].plot(dt_range, times, 'o-', color = colors[integrator], linewidth = 2)#, label = integrator2name[integrator])
axs[0].legend(loc = 'lower right')
plt.tight_layout()
fig.savefig("integrator_performance.pdf")
fig.show()

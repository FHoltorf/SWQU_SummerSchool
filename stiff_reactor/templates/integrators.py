from scipy.optimize import fsolve # use fsolve for implicit time stepping
import time # use time.time() to track integration time 

def explicit_euler_step(x0,t,f,dt):
    # implment explicit euler time stepping
    return None 

def implicit_euler_step(x0,t,f,dt):
    # implement implicit euler time stepping
    return None
    
def crank_nicolson_step(x0,t,f,dt):
    # implement crank-nicolson time stepping
    return None 
    
def heun_step(x0,t,f,dt):
    # implement heun time stepping
    return None
    
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

# Simulate reactor:
# A -> B <-> 2D



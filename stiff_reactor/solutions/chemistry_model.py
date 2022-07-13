import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# 1,2,7,10
class IonosphereChemistry:
    def __init__(self):
            self.spec2idx = {'e' : 0,
                             'O' : 1,
                             'O+' : 2,
                             'O2' : 3,
                             'O2+' : 4,
                             'N' : 5,
                             'N+' : 6,
                             'N2' : 7,
                             'N2+' : 8,
                             'NO' : 9,
                             'NO+' : 10,
                             'He' : 11,
                             'He+' : 12}
            self.idx2spec = {self.spec2idx[key] : key for key in self.spec2idx.keys()}
            stoich_mat = np.zeros((13, 17))
            #electrons
            stoich_mat[0,[6,11,12]] = -1
            #O
            stoich_mat[1,[0,1,2,5,12,14]] = 1
            stoich_mat[1,6] = 2
            stoich_mat[1,[7,8]] = -1
            #O+
            stoich_mat[2,8] = 1
            stoich_mat[2,[0,1,2]] = -1
            #O2
            stoich_mat[3,4] = 1
            stoich_mat[3,[1,9,14,15]] = -1
            #O2+
            stoich_mat[4,[1,9,15]] = 1
            stoich_mat[4,[3,4,5,6]] = -1
            #N
            stoich_mat[5,5] = -1
            stoich_mat[5,[7,11,12,13,15,16]] = [1,2,1,1,1,1]
            #N+
            stoich_mat[6,13] = 1
            stoich_mat[6,[14,15,16]] = -1
            #N2
            stoich_mat[7,[8,9,10]] = 1
            stoich_mat[7,[0,3,13]] = -1
            #N2+
            stoich_mat[8,[7,8,9,10,11]] = -1
            #NO
            stoich_mat[9,3] = 1
            stoich_mat[9,[2,4,10,16]] = -1
            #NO+
            stoich_mat[10,[0,2,3,4,5,7,10,14,16]] = 1
            stoich_mat[10,12] = -1
            #He 
            stoich_mat[11,13] = 1
            #He+
            stoich_mat[12,13] = -1

            self.S = stoich_mat

    def k1(self,T):
        if T >= 300 and T <= 1700:
            k1 = 1.533e-12 - 5.92e-13*(T/300.0) + 8.6e-14*(T/300.0)**2
        else: 
            k1 = 2.73e-12 - 1.155e-12*(T/300.0) + 1.483e-13*(T/300.0)**2 
        return k1

    def k2(self,T):
        if 300 <= T: 
            k2 = 2.82e-11-7.74e-12*(T/300.0)+1.073e-12*(T/300.0)**2-5.17e-14*(T/300.0)**3 + 9.65e-16*(T/300.0)**4
        return k2

    def k3(self,T):
        if 320 < T and T < 1500:
            k3 = 8.36e-13 - 2.02e-13*T/300.0 + 6.95e-14*(T/300.0)**2 
        else: 
            k3 = 5.33e-13 - 1.64e-14*(T/300.0) + 4.72e-14*(T/300.0)**2 - 7.05e-16*(T/300.0)**3 
        return k3

    def rate_constants(self,T):
        ks = np.array([self.k1(T),
              self.k2(T), 
              self.k3(T),
              5e-16,
              4.4e-10,
              1.2e-10,
              1.6e-7*(300/T)**0.55,
              1.4e-10*(300/T)**0.44 if T <= 1500 else 5.2e-11*(T/300)**0.2,
              1e-11*(300/T)**0.23 if T <= 1500 else 3.6e-12*(T/300)**0.41,
              5e-11*300/T,
              3.3e-10,
              1.8e-7*(300/T)**0.39,
              4.2e-7*(300/T)**0.85,
              1e-9,
              2e-10,
              4e-10,
              2e-12])
        return ks 

    def rates(self,c,T):
        ks = self.rate_constants(T)
        aux = ks * [c[2]*c[7], 
              c[2]*c[3], 
              c[2]*c[9], 
              c[4]*c[7],
              c[4]*c[9], 
              c[4]*c[5], 
              c[4]*c[0], 
              c[8]*c[1],
              c[8]*c[1],
              c[8]*c[3],
              c[8]*c[9],
              c[8]*c[0],
              c[10]*c[0],
              c[12]*c[7],
              c[6]*c[3],
              c[6]*c[3],
              c[6]*c[9]]
        return aux

def reactor_model(t,c,chem_model,D,c_in,T):
    return D*(c_in - c) + chem_model.S @ chem_model.rates(c,T)

chem_model = IonosphereChemistry()

D = 0.0 # [-]
T = 2000 # [K] 
scale = 1e8
ions = [0, 2, 4, 6, 8, 10, 12]
c0 = np.zeros(13)
c0[chem_model.spec2idx['O2']] = 0.21*scale #0.21*p/(8.314*T)
c0[chem_model.spec2idx['N2']] = 0.79*scale #0.79*p/(8.314*T)
c0[ions[1:]] = scale*np.ones(6)
c0[0] = sum(c0[ions]) 

def my_reactor_model(t,c):
    return chem_model.S @ chem_model.rates(c,T)

tspan_short = (0.0, 2.0)
tspan_long = (0.0, 50.0)

sol = solve_ivp(my_reactor_model, tspan_long, c0)
fig, ax = plt.subplots()
for i in ions:
    ax.plot(sol.t, sol.y[i], label = chem_model.idx2spec[i])
ax.legend(loc = 'upper right')
fig.show()
fig.savefig("ion_concentrations_long.pdf")

fig, ax = plt.subplots()
for i in set(range(13)) - set(ions):
    ax.plot(sol.t, sol.y[i], label = chem_model.idx2spec[i])
ax.legend(loc = 'upper right')
fig.show()
fig.savefig("species_concentrations_long.pdf")


sol = solve_ivp(my_reactor_model, tspan_short, c0)
fig, ax = plt.subplots()
for i in ions:
    ax.plot(sol.t, sol.y[i], label = chem_model.idx2spec[i])
ax.legend(loc = 'upper right')
fig.show()
fig.savefig("ion_concentrations_short.pdf")

fig, ax = plt.subplots()
for i in set(range(13)) - set(ions):
    ax.plot(sol.t, sol.y[i], label = chem_model.idx2spec[i])
ax.legend(loc = 'upper right')
fig.show()
fig.savefig("species_concentrations_short.pdf")

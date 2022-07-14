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
                             'N2' : 5,
                             'N2+' : 6,
                             'NO+' : 7,
                             'N' : 8}
            # NO, N, N+, He, He+
            self.idx2spec = {self.spec2idx[key] : key for key in self.spec2idx.keys()}
            
            # reaction network
            # chemical reactions
            # (1) O+ + N2 -> NO+ + N
            # (2) O+ + O2 -> O + O2+
            # (3) O2+ + e -> 2O
            # (4) N2+ + O -> O+ + N2
            # (5) N2+ + O2 -> O2+ + N2
            # (6) O2+ + N -> NO+ + O
            # (7) NO+ + e -> N + O
            # photo-ionization
            # (8/9) O <-> O+ + e
            # (10/11) O2 <-> O2+ + e
            # (12/13) N2 <-> N2+ + e
            
            stoich_mat = np.zeros((9,13))
            # rxn 1: O+ + N2 -> NO+ + N
            stoich_mat[[2,5],0] = - 1
            stoich_mat[[7,8],0] = 1
            # rxn 2: O+ + O2 -> O + O2+
            stoich_mat[[2,3],1] = - 1
            stoich_mat[[1,4],1] = 1
            # rxn 3: O2+ + e -> 2O
            stoich_mat[[0,4],2] = -1
            stoich_mat[1,2] = 2
            # rxn 4: N2+ + O -> O+ + N2
            stoich_mat[[6,1],3] = -1
            stoich_mat[[2,5],3] = 1
            # rxn 5: N2+ + O2 -> O2+ + N2
            stoich_mat[[6,3],4] = -1
            stoich_mat[[4,5],4] = 1
            # rxn 6: O2+ + N -> NO+ + O
            stoich_mat[[4,8],5] = -1
            stoich_mat[[7,1],5] = 1
            # rxn 7: NO+ + e -> N + O
            stoich_mat[[7,0],6] = -1
            stoich_mat[[8,1],6] = 1
            # rxn 8/9: O <-> O+ + e
            stoich_mat[[1],7] = -1
            stoich_mat[[0,2],7] = 1
            stoich_mat[[1],8] = 1
            stoich_mat[[0,2],8] = -1
            # rxn 10/11: O2 <-> O2+ + e
            stoich_mat[[3],9] = -1
            stoich_mat[[0,4],9] = 1
            stoich_mat[[3],10] = 1
            stoich_mat[[0,4],10] = -1
            # rxn 12/13: N2 <-> N2+ + e
            stoich_mat[[5],11] = -1
            stoich_mat[[0,6],11] = 1
            stoich_mat[[5],12] = 1
            stoich_mat[[0,6],12] = -1

            self.S = stoich_mat

    def k1(self,T):
        if T >= 300 and T <= 1700:
            k1 = 1.533e-12 - 5.92e-13*(T/300.0) + 8.6e-14*(T/300.0)**2
        else: 
            k1 = 2.73e-12 - 1.155e-12*(T/300.0) + 1.483e-13*(T/300.0)**2 
        return k1

    def rate_constants(self,T):
        # reaction network
        # chemical reactions
        # (1) O+ + N2 -> NO+ + O 
        # (2) O2+ + e -> 2O
        # (3) N2+ + O -> O+ + N2
        # (4) N2+ + O2 -> O2+ + N2
        # (5) O2+ + N -> NO+ + O
        # (6) NO+ + e -> N + O
        # photo-ionization
        # (7/8) O <-> O+ + e
        # (9/10) O2 <-> O2+ + e
        # (11/12) N2 <-> N2+ + e
        ks = np.array([
                self.k1(T),
                2.82e-11,
                1.6e-7*(300/T)**0.55,
                1e-11*(300/T)**0.23 if T <= 1500 else 3.6e-12*(T/300)**0.41,
                5e-11*300/T,
                1.2e-10,
                4.2e-7*(300/T)**0.85,
                1e-4, 
                1e-8,
                1e-4,
                1e-8,
                1e-4,
                1e-8])
        return ks 

    def rates(self,c,T):
        ks = self.rate_constants(T)
        return ks * [c[2]*c[5], c[2]*c[3], c[0]*c[4], c[6]*c[1], c[6]*c[3], c[4]*c[8], 
                     c[7]*c[0], c[1], c[2]*c[0], c[3], c[4]*c[0], c[5], c[6]*c[0]]

def reactor_model(t,c,chem_model,T):
    return chem_model.S @ chem_model.rates(c,T)

chem_model = IonosphereChemistry()

c0_100 = np.zeros(9)
c0_100[chem_model.spec2idx['O']] = 4.26e11
c0_100[chem_model.spec2idx['O2']] = 2.21e12
c0_100[chem_model.spec2idx['N2']] = 9.22e12
Te_100 = 300

c0_300 = np.zeros(9)
c0_300[chem_model.spec2idx['O']] = 3.21e8
c0_300[chem_model.spec2idx['O2']] = 1.03e6
c0_300[chem_model.spec2idx['N2']] = 2.74e7
Te_300 = 1200

c0 = {100 : c0_100, 300 : c0_300}
Te = {100 : Te_100, 300 : Te_300}

alts = [100, 300]
tspan = {100 : (0, 10000.0), 300 : (0, 40000.0)}
t_eval = {alt : None for alt in alts} 
sol = {alt : solve_ivp(lambda t, c : reactor_model(t,c,chem_model,Te[alt]), tspan[alt], c0[alt], t_eval = t_eval[alt], method = 'LSODA') for alt in alts}

ions = [0,2,4,6,7]
#print(sum(sol[100].y[ions[1:],-1]) - sol[100].y[0,-1])
#print(sol[100].y[:,-1])
#print(sum(sol[100].y[[1,2],-1]) + 2*sum(sol[100].y[[3,4],-1]) - c0[100][1] - 2*c0[100][3])
#print(chem_model.S)
for alt in alts:
    fig, ax = plt.subplots()
    for i in ions:
        ax.plot(sol[alt].t, sol[alt].y[i], label = chem_model.idx2spec[i])
    ax.legend(loc = 'upper right')
    fig.savefig("ion_concentrations_"+str(alt)+".pdf")

    fig, ax = plt.subplots()
    for i in set(range(8)) - set(ions):
        ax.plot(sol[alt].t, sol[alt].y[i]/c0[alt][i], label = chem_model.idx2spec[i])
    ax.legend(loc = 'upper right')
    fig.savefig("species_concentrations_"+str(alt)+".pdf")



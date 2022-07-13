import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp 

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

    # implement class methods for computation of reaction rates
    
def reactor_model(args):
    # implement (closed) isothermal and isochoric CSTR reactor model
    return None

# integrate model for given initial condition and temperature

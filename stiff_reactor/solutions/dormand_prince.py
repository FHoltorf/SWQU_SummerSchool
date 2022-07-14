
class DormandPrince:
    def __init__(self):
        self.a = [[1/5],
                  [3/40, 9/40], 
                  [44/45, -56/15, 32/9], 
                  [19372/6561, -25360/2187, 64448/6561, -212/729],  
                  [9017/3168, -355/33, 46732/5247, 49/176, -5103/18656], 
                  [35/384, 0, 500/1113, 125/192, -2187/6784, 11/84]]
        self.b = [35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0]
        self.c = [0, 1/5, 3/10, 4/5, 8/9, 1, 1]

def explicit_RK_stepper(x,t,f,h,butcher_tableau):
    n = len(butcher_tableau.c)
    ks = [f(x,t)]
    x_new = x + h*butcher_tableau.b[0]*ks[0]
    for i in range(n-1):
        y = x + h*sum(butcher_tableau.a[i][j]*ks[j] for j in range(i+1))
        ks.append(f(y, t+h*butcher_tableau.c[i+1]))
        x_new += h*butcher_tableau.b[i+1]*ks[-1]
    return x_new

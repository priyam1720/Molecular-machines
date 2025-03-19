import numpy as np
import matplotlib.pyplot as plt

def add(u,v):
    i = 0
    Ans = []
    while i < len(u):
        a = u[i]+v[i]
        Ans.append(a)
        i = i + 1
    return Ans

def add1(S):
    s0 = 0
    for x in S:
        s0 = s0 + x
    return s0

def Normalize(S):
    s0 = add1(S)
    L = [0]
    s = 0
    for x in S:
        s = s + x/s0
        L.append(s)
    return L

def gen():
    r1 = np.random.uniform()
    r2 = np.random.uniform()
    return r1,r2

Temp = 283.15
def rate(E):
    kb = 1.38*(10**(-23))
    h = 6.625*(10**(-34))
    A = (kb*Temp)/h
    k = np.exp((-1)*(E*1000)/(1.987*(Temp)))
    return A*k

Rxns = [[-1,1,0,0,-1,1],[0,0,1,-1,-1,1],[0,0,-1,1,0,0],[1,-1,0,0,0,0],[0,-1,1,0,0,0],[0,1,-1,0,0,0]]
G_form = 11.2
delta_G_form = 0.1
G_hyd = 8.6
delta_G_hyd = 0.5
G_mech2 = 4.3
G_mech3 = 15.6
k1 = rate(G_form)
k2 = rate(G_form + delta_G_form)
k3 = rate(G_hyd)
k4 = rate(G_hyd + delta_G_hyd)
k5 = rate(G_mech2)
k6 = rate(G_mech2)
cat = 2400
Volume=600000
def set_rxn(X):
    a,b,c,d,f,w = X
    r1 = k1*a*f/Volume
    r2 = k2*d*f/Volume
    r3 = k3*c*cat/Volume
    r4 = k4*b*cat/Volume
    r5 = k5*b
    r6 = k6*c
    r0 = r1 + r2 + r3 + r4 + r5 + r6
    P = [r1/r0,(r1+r2)/r0,(r1+r2+r3)/r0,(r1+r2+r3+r4)/r0,(r1+r2+r3+r4+r5)/r0]
    return P,r0

def step(X,t):
    P,a0 = set_rxn(X)
    a1,a2,a3,a4,a5 = P
    r1,r2 = gen()
    if r1 < a1:
        rxn = Rxns[0]
    elif r1 < a2:
        rxn = Rxns[1]
    elif r1 < a3:
        rxn = Rxns[2]
    elif r1 < a4:
        rxn = Rxns[3]
    elif r1 < a5:
        rxn = Rxns[4]
    elif r1 < 1:
        rxn = Rxns[5]
    tau = (1/a0)*(np.log(1/r2))
    Ans = add(X,rxn)
    t1 = t + tau
    return Ans,t1,Rxns.index(rxn)

sample = 1
k = 0
S = []
Time = []
No_of_rxns = []
Type_of_rxn = []
# Open a file for writing
output_file = open('simulation_output.txt', 'w')
step_counter = 0  # Counter to track steps
while k < sample:
    X = [24,0,0,24,4800,0]
    t = 0
    T = [t]
    K = [X]
    i = 0
    I = [i]
    Type = []
    while X[4] > 0:
        X,t,n = step(X,t)
        K.append(X)
        T.append(t)
        i = i + 1
        I.append(i)
        Type.append(n)
        step_counter += 1

        # Check if 10000 steps have been completed
        if step_counter % 10000 == 0:
            output_file.write(f"After {step_counter} steps: {X}\n")

    while X[1] > 0 or X[2] > 0:
        X,t,n = step(X,t)
        K.append(X)
        T.append(t)
        i = i + 1
        I.append(i)
        Type.append(n)
        step_counter += 1

        if step_counter % 10000 == 0:
            output_file.write(f"After {step_counter} steps: {X}\n")

    S.append(K)
    Time.append(T)
    No_of_rxns.append(I)
    Type_of_rxn.append(Type)
    k = k + 1
print(X)


As = []
Bs = []
Cs = []
Ds = []
Fs = []
Ws = []
for K in S:
    A = []
    B = []
    C = []
    D = []
    F = []
    W = []
    for X in K:
        a,b,c,d,f,w = X
        A.append(a)
        B.append(b)
        C.append(c)
        D.append(d)
        F.append(f)
        W.append(w)
    As.append(A)
    Bs.append(B)
    Cs.append(C)
    Ds.append(D)
    Fs.append(F)
    Ws.append(W)

i = 0
EE = []
while i < sample:
    a = As[i][-1]
    d = Ds[i][-1]
    x = (d-a)/(a+d)
    EE.append(x*100)
    plt.plot(Time[i],As[i], label='1', color= 'red')
    plt.plot(Time[i],Bs[i], label='2', color= 'blue')
    plt.plot(Time[i],Cs[i], label = '3', color= 'purple')
    plt.plot(Time[i],Ds[i], label='4', color= 'green')
    i = i + 1

m = np.mean(EE)
s = np.std(EE)
CI = (s/(np.sqrt(sample)))
print(f'''
mean enantiomeric excess = {m} +- {CI} %

CI at 90% = {m - 1.64*CI} to {m + 1.64*CI}
CI at 95% = {m - 1.96*CI} to {m + 1.96*CI}

This is at
delta_delta_G_form = {delta_G_form} kcal/mol and
delta_delta_G_hyd = {delta_G_hyd} kcal/mol''')

plt.legend()
plt.show()

i = 0
while i < sample:
    plt.plot(Time[i],Fs[i], label='Fuel', color= 'red')
    plt.plot(Time[i],Ws[i], label='Waste', color= 'blue')
    i = i + 1
plt.legend()
plt.show()


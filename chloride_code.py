import numpy as np
import matplotlib.pyplot as plt



def add(u, v):
    i = 0
    Ans = []
    while i < len(u):
        a = u[i] + v[i]
        if i == 2:  # Apply constraint only to x[2]
            a = max(0, a)  # Ensure x[2] is non-negative
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

Rxns = [[0,-1,1,1],[1,-1,0,1],[1,0,-1,0],[-1,0,1,0]]
G_hyd = 8.8
delta_G_hyd = -0.3
G_mech1 = 13.6
k1 = rate(G_hyd)
k2 = rate(G_hyd + delta_G_hyd)
km1 = rate(G_mech1)
cat = 600
Volume=600000
def set_rxn1(X):
    a,x,d,w = X
    r1 = (k1*x*cat)/(2*Volume)
    r2 = (k2*x*cat)/(2*Volume)
    r3 = km1*d
    r4 = km1*a
    r0 = r1 + r2 + r3 + r4
    P = []
    sumr = 0
    for r in [r1,r2,r3]:
        sumr = sumr + r
        P.append(sumr/r0)
    return P,r0


def step1(X,t):
    P,a0 = set_rxn1(X)
    a1,a2,a3 = P
    r1,r2 = gen()
    m = 0
    if r1 < a1:
        rxn = Rxns[0]
    elif r1 < a2:
        rxn = Rxns[1]
    elif r1 < a3:
        rxn = Rxns[2]
        m = 1
    elif r1 < 1:
        rxn = Rxns[3]
        m = -1
    tau = (1/a0)*(np.log(1/r2))
    Ans = add(X,rxn)
    t1 = t + tau
    return Ans,t1,Rxns.index(rxn),m
    Ans[1] = max(0, Ans[1])

sample = 100
k = 0
S = []
Time = []
No_of_rxns = []
Type_of_rxn = []
Mech_rxns = []
# Open a file for writing
output_file = open('simulation_output.txt', 'w')
step_counter = 0  # Counter to track steps
while k < sample:
    X = [0,600,0,0]
    t = 0
    T = [t]
    K = [X]
    i = 0
    I = [i]
    Type = []
    M = []
    while X[3] < 600 :
        X,t,n,m = step1(X,t)
        K.append(X)
        T.append(t)
        i = i + 1
        I.append(i)
        Type.append(n)
        M.append(m)
        step_counter += 1

        # Check if 10000 steps have been completed
        if step_counter % 10 == 0:
            output_file.write(f"After {step_counter} steps: {X}\n")


    S.append(K)
    Time.append(T)
    No_of_rxns.append(I)
    Type_of_rxn.append(Type)
    Mech_rxns.append(M)
    k = k + 1
print(X)

As = []
Xs = []
Ds = []
Ws = []
for K in S:
    A = []
    X = []
    D = []
    W = []
    for X in K:
        a,x,d,w = X
        A.append(a)
        X.append(x)
        D.append(d)
        W.append(w)
    As.append(A)
    Xs.append(X)
    Ds.append(D)
    Ws.append(W)

i = 0
EE = []
while i < sample:
    a = As[i][-1]
    d = Ds[i][-1]
    if (a + d) == 0:
        x = 0  # or np.nan, depending on your use case
    else:
        x = (d - a) / (a + d)
    EE.append(x * 100)
    plt.plot(Time[i],As[i], label='1', color= 'red')
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
delta_delta_G_hyd = {delta_G_hyd} kcal/mol''')

plt.legend()
plt.show()

Rotations = []
Total_mech_rxns = []
for M in Mech_rxns:
    r = 0
    l = len(M)
    for m in M:
        r = r + m
        if m == 0:
            l = l - 1
    Rotations.append(r)
    Total_mech_rxns.append(l)

print(f"total no. of mechanical rxns = {np.mean(Total_mech_rxns)} +- {np.std(Total_mech_rxns)}")
print(f"D_count = {np.mean(Rotations)} +- {np.std(Rotations)}")
print(f"D_count/total no. of rxns = {(np.mean(Rotations))/(np.mean(Total_mech_rxns))}")

i = 0
while i < sample:
    plt.plot(Time[i],Ws[i], label='Waste', color= 'blue')
    i = i + 1
plt.legend()
plt.show()

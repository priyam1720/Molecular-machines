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

Rxns = [[-1,1,0,-1,1],[0,1,-1,-1,1],[0,-1,1,0,0],[1,-1,0,0,0],[1,0,-1,0,0],[-1,0,1,0,0]]
G_form = 11.2
delta_G_form = 0.1
G_hyd = 8.6
delta_G_hyd = 0.5
G_mech1 = 15.6
Volume = 600000
k1 = rate(G_form)
k2 = rate(G_form + delta_G_form)
k3 = rate(G_hyd)
k4 = rate(G_hyd + delta_G_hyd)
km1 = rate(G_mech1)
cat = 600
def set_rxn1(X):
    a,x,d,f,w = X
    r1 = (k1*a*f)/Volume
    r2 = (k2*d*f)/Volume
    r3 = (k3*x*cat)/(2*Volume)
    r4 = (k4*x*cat)/(2*Volume)
    r5 = km1*d
    r6 = km1*a
    r0 = r1 + r2 + r3 + r4 + r5 + r6
    P = []
    sumr = 0
    for r in [r1,r2,r3,r4,r5]:
        sumr = sumr + r
        P.append(sumr/r0)
    return P,r0


def step1(X,t):
    P,a0 = set_rxn1(X)
    a1,a2,a3,a4,a5 = P
    r1,r2 = gen()
    m = 0
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
        m = 1
    elif r1 < 1:
        rxn = Rxns[5]
        m = -1
    tau = (1/a0)*(np.log(1/r2))
    Ans = add(X,rxn)
    t1 = t + tau
    return Ans,t1,Rxns.index(rxn),m
    Ans[4] = max(0, Ans[4])

S = []
Time = []
No_of_rxns = []
Type_of_rxn = []
Mech_rxns = []

# Initialize lists for additional data
X3_values = []
rxn4_counts = []
rxn5_counts = []
mech_rxn_diff = []
Fuel_at_start = 7500
sample = 1
k = 0
while k < sample:
    X = [300, 0, 300, 7500, 0]
    t = 0
    T = [t]
    K = [X]
    i = 0
    I = [i]
    Type = []
    M = []

    # Initialize counts for each step
    rxn4_count = 0
    rxn5_count = 0
    D = 0

    # Collect data throughout the simulation
    X3_sample = [X[3]]
    rxn4_sample = [rxn4_count]
    rxn5_sample = [rxn5_count]
    mech_diff = [rxn4_count-rxn5_count]

    while X[3] > 0:
        X, t, n, m = step1(X, t)
        K.append(X)
        T.append(t)
        i += 1
        I.append(i)
        Type.append(n)
        M.append(m)
        # Collect the number of reactions
        if n == 4:
            rxn4_count += 1
            d = 1
        elif n == 5:
            rxn5_count += 1
            d = -1
        else:
            d = 0
        D = D + d
        # Append the current values to lists
        X3_sample.append(X[3])
        rxn4_sample.append(rxn4_count)
        rxn5_sample.append(rxn5_count)
        mech_diff.append(d)

    # Append the sample data
    X3_values.append(X3_sample)
    rxn4_counts.append(rxn4_sample)
    rxn5_counts.append(rxn5_sample)
    mech_rxn_diff.append(mech_diff)

    S.append(K)
    Time.append(T)
    No_of_rxns.append(I)
    Type_of_rxn.append(Type)
    Mech_rxns.append(M)
    k += 1

# Calculate the difference between the number of Rxn[4] and Rxn[5]

#4 - r5 for r4, r5 in zip(rxn4_counts, rxn5_counts)
print(X)
D_Rxn4 = []
D_Rxn5 = []
D_mech_diff = []
for i in range(0,Fuel_at_start):
    d_rxn4 = []
    d_rxn5 = []
    d_mech_diff = []
    k = 0
    while k < len(X3_sample):
        if X3_sample[k] == i:
            d_rxn4.append(rxn4_sample[k])
            d_rxn5.append(rxn5_sample[k])
            d_mech_diff.append(mech_diff[k])
        k = k + 1
    D_Rxn4.append(d_rxn4[-1] - d_rxn4[0])
    D_Rxn5.append(d_rxn5[-1] - d_rxn5[0])
    D_mech_diff.append(d_mech_diff[-1] - d_mech_diff[0])

As = []
Bs = []
Ds = []
Fs = []
Ws = []
for K in S:
    A = []
    B = []
    D = []
    F = []
    W = []
    for X in K:
        a,b,d,f,w = X
        A.append(a)
        B.append(b)
        D.append(d)
        F.append(f)
        W.append(w)
    As.append(A)
    Bs.append(B)
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
    plt.plot(Time[i],Bs[i], label='3', color= 'violet')
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

Rotations = []
Total_mech_rxns = []
for M in Mech_rxns:
    r = 0
    l = len(M)
    for m in M:
        r += m
        if m == 0:
            l -= 1
    Rotations.append(r)
    Total_mech_rxns.append(l)

print(f"total no. of mechanical rxns = {np.mean(Total_mech_rxns)} +- {np.std(Total_mech_rxns)}")
print(f"D_count = {np.mean(Rotations)} +- {np.std(Rotations)}")
print(f"D_count/total no. of rxns = {np.mean(Rotations) / np.mean(Total_mech_rxns)}")

# Calculate the difference between the number of Rxn[4] and Rxn[5]
#rxn_diff = [r4 - r5 for r4, r5 in zip(rxn4_counts, rxn5_counts)]

# Print debug information
plt.figure(figsize=(12, 6))
plt.scatter(X3_values, rxn4_counts, s=10, color='red', label='Number of Rxn[4]', alpha=0.7)
plt.xlabel('X[3]')
plt.ylabel('Number of Rxn[4]')
plt.title('X[3] vs Number of Rxn[4]')
plt.legend()
plt.grid(True)
plt.show()

# Plot X[3] vs number of Rxn[5]
plt.figure(figsize=(12, 6))
plt.scatter(X3_values, rxn5_counts, s=10, color='blue', label='Number of Rxn[5]', alpha=0.7)
plt.xlabel('X[3]')
plt.ylabel('Number of Rxn[5]')
plt.title('X[3] vs Number of Rxn[5]')
plt.legend()
plt.grid(True)
plt.show()

# Plot X[3] vs (number of Rxn[4] - number of Rxn[5])
plt.figure(figsize=(12, 6))
plt.scatter(range(0,Fuel_at_start), D_Rxn4, s=10, color='green', label='Number of Rxn[4] - Number of Rxn[5]', alpha=0.7)
plt.xlabel('X[3]')
plt.ylabel('Number of Rxn[4] - Number of Rxn[5]')
plt.title('X[3] vs Number of Rxn[4] - Number of Rxn[5]')
plt.legend()
plt.grid(True)
plt.show()

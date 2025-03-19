import numpy as np
import matplotlib.pyplot as plt

# Helper functions
def add(u, v):
    return [ui + vi for ui, vi in zip(u, v)]

def gen():
    return np.random.uniform(), np.random.uniform()

# Constants
Temp = 283.15
kb = 1.38 * (10 ** (-23))
h = 6.625 * (10 ** (-34))
A = (kb * Temp) / h

def rate(E):
    return A * np.exp((-1) * (E * 1000) / (1.987 * Temp))

# Reaction definitions
Rxns = [[-1, 1, -1, 1], [-1, 1, -1, 1], [1, -1, 0, 0], [1, -1, 0, 0]]
Rxns2 = [[1, -1, 0, 0], [1, -1, 0, 0]]

# Reaction parameters
G_form = 11.2
delta_G_form = 0.1
G_hyd = 8.6
delta_G_hyd = 0.5
k1 = rate(G_form)
k2 = rate(G_form + delta_G_form)
k3 = rate(G_hyd)
k4 = rate(G_hyd + delta_G_hyd)
cat = 600
Volume = 600000

def set_rxn1(X):
    a, x, f, w = X
    r1 = (k1 * a * f) / (2 * Volume)
    r2 = (k2 * a * f) / (2 * Volume)
    r3 = (k3 * x * cat) / (2 * Volume)
    r4 = (k4 * x * cat) / (2 * Volume)
    r0 = r1 + r2 + r3 + r4
    P = [sum([r1, r2, r3, r4][:i+1]) / r0 for i in range(4)]
    return P, r0

def set_rxn2(X):
    a, x, f, w = X
    r1 = (k3 * x * cat) / (2 * Volume)
    r2 = (k4 * x * cat) / (2 * Volume)
    r0 = r1 + r2
    P = [sum([r1, r2][:i+1]) / r0 for i in range(2)]
    return P, r0

def step1(X, t, reaction_counts):
    P, a0 = set_rxn1(X)
    r1, r2 = gen()
    if r1 < P[0]:
        rxn = Rxns[0]
        reaction_counts["Rxn1"] += 1
    elif r1 < P[1]:
        rxn = Rxns[1]
        reaction_counts["Rxn2"] += 1
    elif r1 < P[2]:
        rxn = Rxns[2]
        reaction_counts["Rxn3"] += 1
    else:
        rxn = Rxns[3]
        reaction_counts["Rxn4"] += 1
    tau = (1 / a0) * np.log(1 / r2)
    return add(X, rxn), t + tau, reaction_counts

def step2(X, t, reaction_counts):
    P, a0 = set_rxn2(X)
    r1, r2 = gen()
    if r1 < P[0]:
        rxn = Rxns2[0]
        reaction_counts["Rxn3"] += 1
    else:
        rxn = Rxns2[1]
        reaction_counts["Rxn4"] += 1
    tau = (1 / a0) * np.log(1 / r2)
    return add(X, rxn), t + tau, reaction_counts

# Simulation parameters
sample = 10
max_steps = 15000
reaction_sums = {"Rxn1": 0, "Rxn2": 0, "Rxn3": 0, "Rxn4": 0}
X2_values = []
reaction_data = {"Rxn1": [], "Rxn2": [], "Rxn3": [], "Rxn4": []}

for run in range(sample):
    reaction_counts = {"Rxn1": 0, "Rxn2": 0, "Rxn3": 0, "Rxn4": 0}
    X = [600, 0, 7500, 0]
    t = 0
    step_counter = 0

    while X[2] > 0 and step_counter < max_steps:
        X, t, reaction_counts = step1(X, t, reaction_counts)
        X2_values.append(X[2])
        for key in reaction_counts:
            reaction_data[key].append(reaction_counts[key])
        step_counter += 1

    while X[1] > 0 and step_counter < max_steps:
        X, t, reaction_counts = step2(X, t, reaction_counts)
        X2_values.append(X[2])
        for key in reaction_counts:
            reaction_data[key].append(reaction_counts[key])
        step_counter += 1

    for key in reaction_counts:
        reaction_sums[key] += reaction_counts[key]

# Compute averages
average_rxns = {key: reaction_sums[key] / sample for key in reaction_sums}
print("Average reaction counts over", sample, "simulations:")
for key, value in average_rxns.items():
    print(f"{key}: {value}")

# Plot results
plt.figure(figsize=(10, 6))
for key in reaction_data:
    plt.plot(X2_values, reaction_data[key], label=key)
plt.xlabel("X[2] values")
plt.ylabel("Reaction Counts")
plt.title("Reaction Counts vs X[2] Over Simulations")
plt.legend()
plt.show()

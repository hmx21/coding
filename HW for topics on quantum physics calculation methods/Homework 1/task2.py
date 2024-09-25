import numpy as np
from scipy import sparse
import matplotlib.pyplot as plt

# Parameters
N = 101
F = 0.1
k_0 = np.pi / 2
alpha = 0.15
N_0 = 51

# Hamiltonian matrix
H = np.zeros((N, N))
for i in range(N):
    if i > 0:
        H[i, i - 1] = -1
    if i < N - 1:
        H[i, i + 1] = -1
    H[i, i] = F * (i + 1)

# Diagonalize the Hamiltonian
eigenvalues, eigenvectors = sparse.linalg.eigsh(H, k=N, which="SM")

# Initial wave packet
psi0 = [
    np.exp(-0.5 * (alpha * (i - N_0)) ** 2) * np.exp(1j * k_0 * i) for i in range(N)
]

# Normalize the wave packet
psi0 /= np.linalg.norm(psi0)


# Time evolution
def time_evolution(psi0, eigenvalues, eigenvectors, t):
    # Calculate the coefficients c of each eigenstate
    c = np.dot(eigenvectors.T, psi0)
    psi_t = np.dot(eigenvectors, c * np.exp(-1j * eigenvalues * t))
    return psi_t


# Output the results

# Q1:Show the lowest 10 eigenvalues
print(f"Lowest 10 eigenvalues:{eigenvalues[:10]}")

# Q2:Show |\psi(j, t)|^2 for t=42 and j=10,20,30,40,50
t = 42
j_values = [10, 20, 30, 40, 50]
print("The probability density at t=42 and j=10,20,30,40,50:")
for j in j_values:
    psi_t = time_evolution(psi0, eigenvalues, eigenvectors, t)
    print(f"|\psi({j}, {t})|^2 = {np.abs(psi_t[j]) ** 2}")

# Q3:Plot |\psi(j, t)|^2 as functions of j and t

# Time points and probability density at each site and time
times = np.linspace(0, 100, 500)
prob_density = np.zeros((N, len(times)))

# Calculate the probability density
for idx, t in enumerate(times):
    psi_t = time_evolution(psi0, eigenvalues, eigenvectors, t)
    prob_density[:, idx] = np.abs(psi_t) ** 2

# Plot the probability density as a function of time and site index
plt.imshow(prob_density, extent=[0, times[-1], 0, N], aspect="auto", origin="lower")
plt.colorbar(label="Probability Density")
plt.xlabel("time t")
plt.ylabel("position j")
plt.savefig("2021010167_何铭鑫_hw1_task2_Q3.png")
plt.show()

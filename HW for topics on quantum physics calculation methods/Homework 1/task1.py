from BitOperation import *
from scipy import sparse
import numpy as np

# Parameters
t = 1.0
U = 6.0
V = 0.3
mu = 3.0
n_sites = 6
n_spin = 2
dim = 2 ** (n_sites * n_spin)

# Additional term list(both nearest and second nearest terms)
nearest_list = [(0, 2), (1, 3), (2, 4), (3, 5), (1, 2), (3, 4)]
second_nearest_list = [(0, 1), (2, 3), (1, 4), (4, 5)]

# Store the hamiltonian matrix(non-zero elements) for each subspace, row index and column index and value
HI = []
HJ = []
HV = []


# Seperate the Hilbert space into different subspaces according to the number of spin-up and spin-down electrons
# The number of spin-up/down electrons is from 0 to n_sites
# Function to count the number of spin-up and spin-down electrons
def count_spins(state, n_sites):
    spin_up = 0
    spin_down = 0
    for i in range(n_sites):
        if ReadBit(state, i) == 1:
            spin_up += 1
        if ReadBit(state, i + n_sites) == 1:
            spin_down += 1
    return spin_up, spin_down


# Dictionary to store subspaces
# The key is the tag of subspace, labeled by the number of spin-up and spin-down electrons
# The value is a list of the states in the subspace
subspaces = {}

# Iterate over all possible states
for state in range(dim):
    spin_up, spin_down = count_spins(state, n_sites)
    if (spin_up, spin_down) not in subspaces:
        subspaces[(spin_up, spin_down)] = []
    subspaces[(spin_up, spin_down)].append(state)

# Re-label the states within each subspace, and store the mapping in a dictionary
# The key is the tag of subspace
# The value is a dictionary, where the key is the original state and the value is the new index
relabeled_subspaces = {}
for key, states in subspaces.items():
    relabeled_subspaces[key] = {state: idx for idx, state in enumerate(states)}

# The Hamiltonian matrix is stored in a dictionary
# The key is the tag of subspace
# The value is the Hamiltonian matrix in the basis of the states in the subspace
hamiltonian_matrices = {}

# Construct the Hamiltonian matrix for each subspace
for config in subspaces:
    for state in subspaces[config]:
        # The diagonal elements: U, V and mu terms
        HI.append(relabeled_subspaces[config][state])
        HJ.append(relabeled_subspaces[config][state])
        U_term = (
            sum(
                [
                    ReadBit(state, i) * ReadBit(state, i + n_sites)
                    for i in range(n_sites)
                ]
            )
            * U
        )
        V_term = (
            sum(
                [
                    (ReadBit(state, i) + ReadBit(state, i + n_sites))
                    * (ReadBit(state, j) + ReadBit(state, j + n_sites))
                    for i, j in second_nearest_list
                ]
            )
            * V
        )
        mu_term = PopCntBit(state) * mu
        HV.append(U_term - V_term - mu_term)
        # The off-diagonal elements: t term
        for pair in nearest_list:
            i, j = pair
            # "minus_sign" is used to determine the sign of the hopping term due to the anticommutation relation of fermionic operators
            minus_sign_up = 0
            minus_sign_down = 0
            # Spin up hopping term
            if ReadBit(state, i) != ReadBit(state, j):
                if ReadBit(state, i) == 1:
                    minus_sign_up += PopCntBit(PickBits(state, i, 0))
                    new_state = FlipBit(state, i)
                    minus_sign_up += PopCntBit(PickBits(new_state, j, 0))
                    new_state = FlipBit(new_state, j)
                else:
                    minus_sign_up += PopCntBit(PickBits(state, j, 0))
                    new_state = FlipBit(state, j)
                    minus_sign_up += PopCntBit(PickBits(new_state, i, 0))
                    new_state = FlipBit(new_state, i)
                HI.append(relabeled_subspaces[config][new_state])
                HJ.append(relabeled_subspaces[config][state])
                HV.append(-t * (-1) ** minus_sign_up)
            # Spin down hopping term
            if ReadBit(state, i + n_sites) != ReadBit(state, j + n_sites):
                if ReadBit(state, i + n_sites) == 1:
                    minus_sign_down += PopCntBit(PickBits(state, i + n_sites, 0))
                    new_state = FlipBit(state, i + n_sites)
                    minus_sign_down += PopCntBit(PickBits(new_state, j + n_sites, 0))
                    new_state = FlipBit(new_state, j + n_sites)
                else:
                    minus_sign_down += PopCntBit(PickBits(state, j + n_sites, 0))
                    new_state = FlipBit(state, j + n_sites)
                    minus_sign_down += PopCntBit(
                        PickBits(new_state, i + n_sites - 1, 0)
                    )
                    new_state = FlipBit(new_state, i + n_sites)
                HI.append(relabeled_subspaces[config][new_state])
                HJ.append(relabeled_subspaces[config][state])
                HV.append(-t * (-1) ** minus_sign_down)
    # Generate the Hamiltonian matrix for the subspace
    hamiltonian_matrices[config] = sparse.coo_matrix(
        (HV, (HI, HJ)), shape=(len(subspaces[config]), len(subspaces[config]))
    ).tocsc()
    # Clear the lists for the next subspace
    HI.clear()
    HJ.clear()
    HV.clear()

# Find the eigenvalues of the subspace (2, 4)
subspace_key = (2, 4)
eigenvalues_of_subspace, _ = sparse.linalg.eigsh(
    hamiltonian_matrices[subspace_key].toarray(),
    k=len(subspaces[subspace_key]),
    which="SM",
)

# Store the ground state energy of the whole system
# The ground state energy is stored as a list, as there may be degenerate ground states
# Each element in the list is a tuple, where the first element is the energy and the second element is the eigenvector, the third element is the tag of the subspace
ground_states = []

# Find all the eigenvalues of the Hamiltonian and sort them in increasing order, as well as find the ground state
eigenvalues = []
for config in subspaces:
    if config in hamiltonian_matrices:
        eigvals, eigstats = sparse.linalg.eigsh(
            hamiltonian_matrices[config].toarray(), k=len(subspaces[config]), which="SM"
        )
        eigenvalues.extend(eigvals.tolist())
        # Store the ground state energy and corresponding eigenstate, as well as the tag of the subspace
        if len(ground_states) == 0:
            ground_states.append((eigvals[0], eigstats[:, 0], config))
        elif eigvals[0] < min([v[0] for v in ground_states]):
            ground_states.clear()
            ground_states.append((eigvals[0], eigstats[:, 0], config))
        elif eigvals[0] == min([v[0] for v in ground_states]):
            ground_states.append((eigvals[0], eigstats[:, 0], config))
eigenvalues.sort()

# Store the expectation values of the number of spin-up and spin-down electrons for the ground state at each site
Density_spin_up = [0] * n_sites
Density_spin_down = [0] * n_sites

# Tag for the ground state
gs_tag = 0

# Output the results

# Q1:Print the eigenvalues of the subspace (2, 4)
print(
    f"6 smallest eigenvalues of the subspace {subspace_key}: {eigenvalues_of_subspace[:6]}"
)

# Q2:Print the first 20 eigenvalues
print(f"20 smallest eigenvalues: {eigenvalues[:20]}")

# Check if there are degenerate ground states
if len(ground_states) > 1:
    print("There are degenerate ground states.")
else:
    print("There is only one ground state.")

# Q3:Print the expectation values of the number of spin-up and spin-down electrons for the ground state at each site
for gs in ground_states:
    gs_tag += 1
    for site in range(n_sites):
        Density_spin_up[site] = 0
        Density_spin_down[site] = 0
        # Using the tag of subspace, we can find the initial state of each relabeld index
        # Here gs[2] is the tag of the subspace, and relabeled_subspaces[gs[2]] is the dictionary that stores the mapping from the initial state to the relabeled index
        for state, relabel in relabeled_subspaces[gs[2]].items():
            Density_spin_up[site] += np.abs(gs[1][relabel]) ** 2 * ReadBit(state, site)
            Density_spin_down[site] += np.abs(gs[1][relabel]) ** 2 * ReadBit(
                state, site + n_sites
            )
    # Print the expectation values of the number of spin-up and spin-down electrons for the ground state at each site
    # This step is to convert the numpy float type to python float type
    Density_spin_down = [x.item() for x in Density_spin_down]
    Density_spin_up = [x.item() for x in Density_spin_up]
    print(
        f"For ground state {gs_tag}: Density of spin-up electrons: {Density_spin_up}\n"
        f"Density of spin-down electrons: {Density_spin_down}"
    )

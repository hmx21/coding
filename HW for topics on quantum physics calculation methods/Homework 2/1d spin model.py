import numpy as np
from scipy import sparse
from BitOperation import *

# Parameters
N = 8
g = 1
J = 1
h = 1


# Generate the check_list
def Get_Check_List(N) -> tuple[list, dict]:
    """
    Generate a list of all possible states of N spins

    :param N: Number of spins

    :return check_list: A list.\n
        Each element is a list containing 3 numbers:\n
        1. whether the state is a representative configuration, 1 for yes, 0 for no\n
        2. the corresponding representative configuration\n
        3. the translation step R\n
    :return representative_data: A dictionary.\n
        The key is the representative configuration.\n
        The value is the number of states which share the same representative configuration
    """

    # initialize the check_list
    check_list = []
    for i in range(2**N):
        check_list.append([1, i, 0])

    # initialize the representative_data
    representative_data = {}

    # generate the check_list
    for i in range(2**N):
        if check_list[i][0] == 1:
            for j in range(1, N + 1):
                new_state = RotLBit(i, N, j)
                if new_state == i:
                    representative_data[i] = j
                    break
                else:
                    check_list[new_state][0] = 0
                    check_list[new_state][1] = i
                    check_list[new_state][2] = j
    return check_list, representative_data


# Separation of the Hilbert space, according to the eigenvalues of the translation operator
def Get_Momentum_Sector(check_list, representative_data) -> dict:
    """
    Apply the P_k operator to the representative configuration of each state
    to get the eigenvector of translation operator with eigenvalue 'exp(-iK)'
    where P_k = 1/N sum_{j=0}^{N-1} exp(i*2*pi*k*j/N) T^j\n

    :param check_list: A list of the data of all possible states of N spins\n
    :param representative_data: A dictionary of the representative configurations and the number of states sharing the same representative configuration\n

    :return momentum_sector: A dictionary of momentum sectors.\n
        The key is the eigenvalue of the translation operator(represented by k)\n
        The value is a list of tuple of eigenstates in the momentum sector\n
        For each tuple:\n
        The first element is the representative configuration\n
        The second element is a list, with each element representing the coefficient of the corresponding basis state\n
        the basis states share the same representative configuration(note that the basis here has the binary form!!!)\n
    """

    # initialize the momentum_sectors
    momentum_sectors = {}

    # generate the momentum_sectors
    for k in range(N):
        # initialize the momentum_sectors[k]
        momentum_sectors[k] = []
        P_k_unit = np.exp(2j * np.pi * k / N)
        for i in check_list:
            # check if the state is the representative configuration
            if i[0] == 1:
                cycle = N / representative_data[i[1]]
                P_k_cycle = np.exp(2j * np.pi * representative_data[i[1]] * k / N)
                # initialize the eigenstate of T
                eigenstate_of_T = [0] * representative_data[i[1]]
                for j in check_list:
                    if j[1] == i[1]:
                        eigenstate_of_T[j[2]] += sum(
                            [
                                P_k_unit ** j[2] * P_k_cycle**n / N
                                for n in range(int(cycle))
                            ]
                        )
                norm = np.linalg.norm(eigenstate_of_T)
                # check the eigenstate is not zero
                if norm > 1e-10:
                    momentum_sectors[k].append((i[1], eigenstate_of_T))
    return momentum_sectors


# Relabel each momentum subspaces
def Get_Relabeled_Subspaces(momentum_sectors) -> dict:
    """
    Relabel each momentum subspaces

    :param momentum_sectors: A dictionary of momentum sectors.\n
        The key is the eigenvalue of the translation operator(represented by k)\n
        The value is a list of tuple of eigenstates in the momentum sector\n
        For each tuple:\n
        The first element is the representative configuration\n
        The second element is a list, with each element representing the coefficient of the corresponding basis state\n
        the basis states share the same representative configuration\n

    :return relabeled_subspaces: A dictionary.\n
        The key is the eigenvalue of the translation operator(represented by k)\n
        The value is a dictionary.\n
            The key is the relabeled index of eigenstate in the subspace\n
            The value is a tuple.\n
                The first element is the representative configuration\n
                The second element is a list, with each element representing the coefficient of the corresponding basis state\n
                the basis states share the same representative configuration.
    """

    # initialize the relabeled_subspaces
    relabeled_subspaces = {}

    # generate the relabeled_subspaces
    for key, states in momentum_sectors.items():
        relabeled_subspaces[key] = {idx: state for idx, state in enumerate(states)}
    return relabeled_subspaces


def Find_Element_Index(check_list, check_info) -> int:
    """
    Find the index of the element in check_list which has the same state as the input information.\n

    :param check_list: A list of the data of all possible states of N spins\n
    :param check_info: A tuple of the state to be found\n
        The first element is the representative configuration\n
        The second element is the translation step R\n

    :return: The index of the element in check_list which has the same state as the input state
    """
    for idx, i in enumerate(check_list):
        if i[1] == check_info[0] and i[2] == check_info[1]:
            return idx


# Generate the Hamiltonian matrix for each momentum sector
def Get_Subspaces_Hamiltonian(check_list, relabeled_subspaces) -> dict:
    """
    Generate the Hamiltonian matrix for each momentum sector

    :param check_list: A list of the data of all possible states of N spins\n
    :param relabeled_subspaces: A dictionary of momentum sectors.\n
        The key is the eigenvalue of the translation operator(represented by k)\n
        The value is a dictionary.\n
            The key is the relabeled index of eigenstate in the subspace\n
            The value is a tuple.\n
                The first element is the representative configuration\n
                The second element is a list, with each element representing the coefficient of the corresponding basis state\n
                the basis states share the same representative configuration.

    :return hamiltonian_matrices: A dictionary.\n
        The key is the eigenvalue of the translation operator(represented by k)\n
        The value is a sparse matrix representing the Hamiltonian matrix in the basis of the states in the relabeled_subspace
    """

    # initialize the hamiltonian_matrices
    hamiltonian_matrices = {}

    # generate the hamiltonian_matrices
    for momentum in relabeled_subspaces:
        HI = []
        HJ = []
        HV = []
        for col_index, col_state_info in relabeled_subspaces[momentum].items():
            representative_config_col = col_state_info[0]
            col_state = col_state_info[1]
            # Diagonal elements
            value = Get_Diagonal_Element(
                check_list, representative_config_col, col_state
            )
            HI.append(col_index)
            HJ.append(col_index)
            HV.append(value)
            for row_index, row_state_info in relabeled_subspaces[momentum].items():
                representative_config_row = row_state_info[0]
                row_state = row_state_info[1]
                # Off-diagonal elements
                if row_index != col_index:
                    value = Get_Off_Diagonal_Element(
                        check_list,
                        representative_config_col,
                        representative_config_row,
                        col_state,
                        row_state,
                    )
                    HI.append(row_index)
                    HJ.append(col_index)
                    HV.append(value)
        # Store the Hamiltonian matrix for each momentum sector
        hamiltonian_matrices[momentum] = sparse.coo_matrix(
            (HV, (HI, HJ)),
            shape=(
                len(relabeled_subspaces[momentum]),
                len(relabeled_subspaces[momentum]),
            ),
        ).tocsc()

    return hamiltonian_matrices


def Get_Diagonal_Element(check_list, representative_config_col, col_state) -> complex:
    """
    Get the diagonal element of the Hamiltonian matrix\n

    :param check_list: A list of the data of all possible states of N spins\n
    :param representative_config_col: The representative configuration of the state\n
    :param col_state: The coefficient of the basis state\n

    :return: The diagonal element of the Hamiltonian matrix
    """
    # initialize the value
    value = 0

    for translation_step, coefficient in enumerate(col_state):
        # find the index of the basis state
        base_state = Find_Element_Index(
            check_list, (representative_config_col, translation_step)
        )
        value += (
            -J
            * sum(
                [
                    2
                    * (ReadBit(base_state, i) - 1 / 2)
                    * 2
                    * (ReadBit(base_state, (i + 1) % N) - 1 / 2)
                    for i in range(N)
                ]
            )
            * np.linalg.norm(coefficient) ** 2
        )

    # normalization
    value /= np.linalg.norm(col_state) ** 2
    return value


def Get_Off_Diagonal_Element(
    check_list,
    representative_config_col,
    representative_config_row,
    col_state,
    row_state,
) -> complex:
    """
    Get the off-diagonal element of the Hamiltonian matrix\n

    :param check_list: A list of the data of all possible states of N spins\n
    :param representative_config_col: The representative configuration of the column state\n
    :param representative_config_row: The representative configuration of the row state\n
    :param col_state: The coefficient of the column basis state\n
    :param row_state: The coefficient of the row basis state\n

    :return: The off-diagonal element of the Hamiltonian matrix
    """
    # initialize the value
    value = 0

    for translation_step, coefficient in enumerate(row_state):
        # find the index of the basis state
        base_state = Find_Element_Index(
            check_list, (representative_config_row, translation_step)
        )
        coe_conjugate = np.conjugate(coefficient)
        # apply off-diagonal Hamiltonian to the column state
        for i in range(N):
            new_state = FlipBit(representative_config_col, i)
            # compute the contribution of each basis state of the row state
            if new_state == base_state:
                # h term with pauli_x
                value += -h * coe_conjugate
                # g term with pauli_y
                if ReadBit(representative_config_col, i) == 1:
                    value += (
                        g
                        * 1.0j
                        * 2
                        * (ReadBit(representative_config_col, (i + 1) % N) - 1 / 2)
                        * 2
                        * (ReadBit(new_state, (i - 1) % N) - 1 / 2)
                        * coe_conjugate
                    )
                else:
                    value += (
                        -g
                        * 1.0j
                        * 2
                        * (ReadBit(representative_config_col, (i + 1) % N) - 1 / 2)
                        * 2
                        * (ReadBit(new_state, (i - 1) % N) - 1 / 2)
                        * coe_conjugate
                    )
    # normalization
    value /= np.linalg.norm(col_state) * np.linalg.norm(row_state)
    return value


def Find_All_Eigenvalues(hamiltonian_matrices, momentum_sectors) -> tuple[dict, list]:
    """
    Find all eigenvalues of the Hamiltonian matrix in each momentum sector

    :param hamiltonian_matrices: A dictionary of momentum sectors.\n
        The key is the eigenvalue of the translation operator(represented by k)\n
        The value is a sparse matrix representing the Hamiltonian matrix in the basis of the states in the relabeled_subspace\n
    :param momentum_sectors: A dictionary of momentum sectors.\n
        The key is the eigenvalue of the translation operator(represented by k)\n
        The value is a list of tuple of eigenstates in the momentum sector\n
        For each tuple:\n
        The first element is the representative configuration\n
        The second element is a list, with each element representing the coefficient of the corresponding basis state\n
        the basis states share the same representative configuration\n

    :return all_eigenvalues: A dictionary.\n
        The key is the eigenvalue of the translation operator(represented by k)\n
        The value is a list of all eigenvalues in the momentum sector\n
    :return ground_states: A list of the ground states\n
        Each element is a tuple, where the first element is the energy, the second element is the eigenvector, the third element is the tag of the subspace(momentum)
    """
    # initialize the ground_states and all_eigenvalues
    ground_states = []  # Store in a list in case of degenerate ground states
    all_eigenvalues = {}
    # store all the eigenvalues in order to find the ground state
    eigenvalues = []
    # find eigenvalues of the Hamiltonian matrix in each momentum sector
    for momentum in momentum_sectors:
        eigvals, eigstats = sparse.linalg.eigsh(
            hamiltonian_matrices[momentum].toarray(),
            k=len(momentum_sectors[momentum]),
            which="SM",
        )
        eigenvalues.extend(eigvals.tolist())
        all_eigenvalues[momentum] = sorted(eigvals.tolist())
        # Store the ground state energy and corresponding eigenstate, as well as the tag of the subspace(momentum)
        if len(ground_states) == 0:
            ground_states.append((eigvals[0], eigstats[:, 0], momentum))
        elif eigvals[0] < min([v[0] for v in ground_states]):
            ground_states.clear()
            ground_states.append((eigvals[0], eigstats[:, 0], momentum))
        elif eigvals[0] == min([v[0] for v in ground_states]):
            ground_states.append((eigvals[0], eigstats[:, 0], momentum))
    return all_eigenvalues, ground_states


def Find_Ground_State(ground_states, relabeled_subspaces, check_list) -> list:
    """
    Find the ground state of the whole system, translate the eigenvector of subspace into the whole system

    :param ground_states: A list of the ground states\n
        Each element is a tuple, where the first element is the energy, the second element is the eigenvector, the third element is the tag of the subspace(momentum)\n
    :param relabeled_subspaces: A dictionary of momentum sectors.\n
        The key is the eigenvalue of the translation operator(represented by k)\n
        The value is a dictionary.\n
            The key is the relabeled index of eigenstate in the subspace\n
            The value is a tuple.\n
                The first element is the representative configuration\n
                The second element is a list, with each element representing the coefficient of the corresponding basis state\n
                the basis states share the same representative configuration.\n
    :param check_list: A list of the data of all possible states of N spins\n

    :return real_ground_states: A list of the ground states\n
        Each element is a tuple, where the first element is the energy, the second element is the eigenvector\n
    """
    # store the ground state of the whole system
    real_ground_states = []
    for gs in ground_states:
        # initialize one of the ground state of the whole system
        real_gs = np.zeros(2**N, dtype=complex)
        ground_energy = gs[0]
        # store the coefficients of each basis state of the subspace of the ground state
        coefficients = gs[1].tolist()
        momentum = gs[2]
        for i in range(len(coefficients)):
            representative_config = relabeled_subspaces[momentum][i][0]
            eigenstates_in_momentum_sector = relabeled_subspaces[momentum][i][1]
            # translate the eigenvector of subspace into the whole system
            for translation_step, coefficient in enumerate(
                eigenstates_in_momentum_sector
            ):
                base_state = Find_Element_Index(
                    check_list, (representative_config, translation_step)
                )
                real_gs[base_state] += (
                    coefficients[i]
                    * coefficient
                    / np.linalg.norm(eigenstates_in_momentum_sector)
                )
        real_ground_states.append((ground_energy, real_gs))

    return real_ground_states


def Get_State_Sigma(real_gs) -> tuple[complex, complex]:
    """
    Get the expectation value of sigma_x and sigma_z for the ground state

    :param real_gs: The ground state of the whole system

    :return sigma_x: The expectation value of sigma_x
    :return sigma_z: The expectation value of sigma_z
    """

    # store the expectation values of sigma_x and sigma_z for each site
    sigma_xs = []
    sigma_zs = []
    # compute the expectation values of sigma_x and sigma_z for each site
    for site in range(N):
        sigma_x_of_each_site = 0 + 0j
        sigma_z_of_each_site = 0 + 0j
        for state in range(2**N):
            new_state = FlipBit(state, site)
            sigma_x_of_each_site += np.conjugate(real_gs[state]) * real_gs[new_state]
            sigma_z_of_each_site += (
                (np.linalg.norm(real_gs[state]) ** 2)
                * 2
                * (ReadBit(state, site) - 1 / 2)
            )

        sigma_xs.append(sigma_x_of_each_site)
        sigma_zs.append(sigma_z_of_each_site)
    # average the expectation values of sigma_x and sigma_z for all sites
    sigma_x = sum(sigma_xs) / N
    sigma_z = sum(sigma_zs) / N
    return sigma_x, sigma_z


# Get the check_list
check_list, representative_data = Get_Check_List(N)

# Get the momentum sectors
momentum_sectors = Get_Momentum_Sector(check_list, representative_data)

# Get the relabeled subspaces
relabeled_subspaces = Get_Relabeled_Subspaces(momentum_sectors)

# Get the Hamiltonian matrix for each momentum sector
hamiltonian_matrices = Get_Subspaces_Hamiltonian(check_list, relabeled_subspaces)

# Find all the eigenvalues of the Hamiltonian, as well as find the ground state
all_eigenvalues, ground_states = Find_All_Eigenvalues(
    hamiltonian_matrices, momentum_sectors
)

# Find the ground state of the whole system(for given case, the ground state is unique)
real_gs = Find_Ground_State(ground_states, relabeled_subspaces, check_list)[0][1]

# Q1: calculate the dispersion relation of this model
for momentum in all_eigenvalues:
    print(f"for momentum sector {momentum}: {all_eigenvalues[momentum]}")

# Q2: calculate the ground state energy of each site, and the expectation value of sigma_x and sigma_z for the ground state
sigma_x, sigma_z = Get_State_Sigma(real_gs)
print(f"ground state energy of each site: {ground_states[0][0]/N}")
print(f"expectation value of sigma_x: {sigma_x}")
print(f"expectation value of sigma_z: {sigma_z}")

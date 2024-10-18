#!/usr/bin/env python
#
# Simple DMRG tutorial.  This code integrates the following concepts:
#  - Infinite system algorithm
#  - Finite system algorithm
#
# Copyright 2013 James R. Garrison and Ryan V. Mishmash.
# Open source under the MIT license.  Source code at
# <https://github.com/simple-dmrg/simple-dmrg/>

# This code will run under any version of Python >= 2.6.  The following line
# provides consistency between python2 and python3.
from __future__ import print_function, division  # requires Python >= 2.6

# numpy and scipy imports
import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import kron, identity
from scipy.sparse.linalg import eigsh  # Lanczos routine from ARPACK
from scipy.optimize import curve_fit

# We will use python's "namedtuple" to represent the Block and EnlargedBlock
# objects
from collections import namedtuple

Block = namedtuple("Block", ["length", "basis_size", "operator_dict"])
EnlargedBlock = namedtuple("EnlargedBlock", ["length", "basis_size", "operator_dict"])


def is_valid_block(block):
    for op in block.operator_dict.values():
        if op.shape[0] != block.basis_size or op.shape[1] != block.basis_size:
            return False
    return True


# This function should test the same exact things, so there is no need to
# repeat its definition.
is_valid_enlarged_block = is_valid_block

# Model-Syecific code for the Homework xx and yy couple chain
model_d = 2  # single-site basis size
N = 40  # length of chain

Sx1 = np.array([[0, 1], [1, 0]], dtype="d")  # single-site Pauli X
Sy1 = np.array([[0, -1j], [1j, 0]], dtype="complex")  # single-site Pauli Y

H1 = np.array([[0, 0], [0, 0]], dtype="d")  # single-site portion of H is zero


def H2(Sx1, Sy1, Sx2, Sy2, g):  # two-site part of H
    """Given the operators S^x and S^y on two sites in different Hilbert Syaces
    (e.g. two blocks), returns a Kronecker product representing the
    correSyonding two-site term in the Hamiltonian that joins the two sites.
    """
    return -g * kron(Sy1, Sy2) + -1 * kron(Sx1, Sx2)


# conn refers to the connection operator, that is, the operator on the edge of
# the block, on the interior of the chain.  We need to be able to represent S^x
# and S^y on that site in the current basis in order to grow the chain.
initial_block = Block(
    length=1,
    basis_size=model_d,
    operator_dict={
        "H": H1,
        "conn_Sx": Sx1,
        "conn_Sy": Sy1,
    },
)


def enlarge_block(block, g):
    """This function enlarges the provided Block by a single site, returning an
    EnlargedBlock.
    """
    mblock = block.basis_size
    o = block.operator_dict

    # Create the new operators for the enlarged block.  Our basis becomes a
    # Kronecker product of the Block basis and the single-site basis.  NOTE:
    # `kron` uses the tensor product convention making blocks of the second
    # array scaled by the first.  As such, we adopt this convention for
    # Kronecker products throughout the code.
    enlarged_operator_dict = {
        "H": kron(o["H"], identity(model_d))
        + kron(identity(mblock), H1)
        + H2(o["conn_Sx"], o["conn_Sy"], Sx1, Sy1, g),
        "conn_Sx": kron(identity(mblock), Sx1),
        "conn_Sy": kron(identity(mblock), Sy1),
    }

    return EnlargedBlock(
        length=(block.length + 1),
        basis_size=(block.basis_size * model_d),
        operator_dict=enlarged_operator_dict,
    )


def rotate_and_truncate(operator, transformation_matrix):
    """Transforms the operator to the new (possibly truncated) basis given by
    `transformation_matrix`.
    """
    return (
        transformation_matrix.conjugate()
        .transpose()
        .dot(operator.dot(transformation_matrix))
    )


def single_dmrg_step(sys, env, m, g):
    """Performs a single DMRG step using `sys` as the system and `env` as the
    environment, keeping a maximum of `m` states in the new basis.
    """
    assert is_valid_block(sys)
    assert is_valid_block(env)

    # Enlarge each block by a single site.
    sys_enl = enlarge_block(sys, g)
    if sys is env:  # no need to recalculate a second time
        env_enl = sys_enl
    else:
        env_enl = enlarge_block(env, g)

    assert is_valid_enlarged_block(sys_enl)
    assert is_valid_enlarged_block(env_enl)

    # Construct the full superblock Hamiltonian.
    m_sys_enl = sys_enl.basis_size
    m_env_enl = env_enl.basis_size
    sys_enl_op = sys_enl.operator_dict
    env_enl_op = env_enl.operator_dict
    superblock_hamiltonian = (
        kron(sys_enl_op["H"], identity(m_env_enl))
        + kron(identity(m_sys_enl), env_enl_op["H"])
        + H2(
            sys_enl_op["conn_Sx"],
            sys_enl_op["conn_Sy"],
            env_enl_op["conn_Sx"],
            env_enl_op["conn_Sy"],
            g,
        )
    )

    # Call ARPACK to find the superblock ground state.  ("SA" means find the
    # "smallest in amplitude" eigenvalue.)
    (energy,), psi0 = eigsh(superblock_hamiltonian, k=1, which="SA")

    # Construct the reduced density matrix of the system by tracing out the
    # environment
    #
    # We want to make the (sys, env) indices correSyond to (row, column) of a
    # matrix, reSyectively.  Since the environment (column) index updates most
    # quickly in our Kronecker product structure, psi0 is thus row-major ("C
    # style").
    psi0 = psi0.reshape([sys_enl.basis_size, -1], order="C")
    rho = np.dot(psi0, psi0.conjugate().transpose())

    # Diagonalize the reduced density matrix and sort the eigenvectors by
    # eigenvalue.
    evals, evecs = np.linalg.eigh(rho)

    # Calculate the raletive entropy of the reduced density matrix
    entanglement_entropy = -sum([x * np.log(x) for x in evals if x > 1e-15])

    possible_eigenstates = []
    for eval, evec in zip(evals, evecs.transpose()):
        possible_eigenstates.append((eval, evec))
    possible_eigenstates.sort(
        reverse=True, key=lambda x: x[0]
    )  # largest eigenvalue first

    # Build the transformation matrix from the `m` overall most significant
    # eigenvectors.
    my_m = min(len(possible_eigenstates), m)
    transformation_matrix = np.zeros(
        (sys_enl.basis_size, my_m), dtype="complex", order="F"
    )
    for i, (eval, evec) in enumerate(possible_eigenstates[:my_m]):
        transformation_matrix[:, i] = evec

    # truncation_error = 1 - sum([x[0] for x in possible_eigenstates[:my_m]])
    # print("truncation error:", truncation_error)

    # Rotate and truncate each operator.
    new_operator_dict = {}
    for name, op in sys_enl.operator_dict.items():
        new_operator_dict[name] = rotate_and_truncate(op, transformation_matrix)

    newblock = Block(
        length=sys_enl.length, basis_size=my_m, operator_dict=new_operator_dict
    )

    return newblock, energy, entanglement_entropy


# Input additional coupling parameter g
def finite_system_algorithm(L, m_warmup, m_sweep, g):
    assert L % 2 == 0  # require that L is an even number

    # To keep things simple, this dictionary is not actually saved to disk, but
    # we use it to represent persistent storage.
    block_disk = {}  # "disk" storage for Block objects

    # Use the infinite system algorithm to build up to desired size.  Each time
    # we construct a block, we save it for future reference as both a left
    # ("l") and right ("r") block, as the infinite system algorithm assumes the
    # environment is a mirror image of the system.
    block = initial_block
    block_disk["l", block.length] = block
    block_disk["r", block.length] = block
    while 2 * block.length < L:
        # Perform a single DMRG step and save the new Block to "disk"
        (block, energy, entanglement_entropy) = single_dmrg_step(
            block, block, m=m_warmup, g=g
        )
        # print("E/L =", energy / (block.length * 2))
        block_disk["l", block.length] = block
        block_disk["r", block.length] = block

    # Now that the system is built up to its full size, we perform sweeps using
    # the finite system algorithm.  At first the left block will act as the
    # system, growing at the expense of the right block (the environment), but
    # once we come to the end of the chain these roles will be reversed.
    sys_label, env_label = "l", "r"
    sys_block = block
    del block  # rename the variable

    # Use the sweep_num to store the number of sweeps, use the energy to justify whether the energy converges
    total_energy = 0
    delta_energy_rate = None
    sweep_num = -1

    # To store the entanglement entropy
    EE_list = [0 for _ in range(L - 3)]
    EE_dict = {}

    while True:
        # Load the appropriate environment block from "disk"
        env_block = block_disk[env_label, L - sys_block.length - 2]
        if env_block.length == 1:
            # We've come to the end of the chain, so we reverse course.
            sys_block, env_block = env_block, sys_block
            sys_label, env_label = env_label, sys_label

        # Use the delta_energy_rate to check whether the energy converges
        if sys_label == "l" and 2 * sys_block.length == L:
            sweep_num += 1
            if sweep_num > 0:
                total_energy += energy
                if sweep_num > 1:
                    average_energy = total_energy / sweep_num
                    delta_energy_rate = np.abs(energy - average_energy) / np.abs(
                        average_energy
                    )

        # Perform a single DMRG step.
        sys_block, energy, entanglement_entropy = single_dmrg_step(
            sys_block, env_block, m=m_sweep, g=g
        )
        EE_list[sys_block.length - 2] = float(entanglement_entropy)

        # Save the block from this step to disk.
        block_disk[sys_label, sys_block.length] = sys_block

        # If the energy is convert (delta_energy_rate < 1e-6) or we have run the full sweep.
        if (
            sys_label == "l"
            and 2 * sys_block.length == L
            and delta_energy_rate is not None
            and delta_energy_rate < 1e-6
        ) or sweep_num > 100:
            EE_list[-1] = EE_list[0]
            EE_dict[(L, m_sweep, g)] = EE_list
            print(f"(N, m, g)=({L}, {m_sweep}, {g})")
            print(f"Energy: {energy}")
            print(f"L: {[i+2 for i in range(L-3)]}")
            print(f"EE: {EE_list}")
            break  # escape from the "while True" loop

    return EE_dict


# Define the function to fit the entanglement entropy
def S(L, c, c_prime):
    return c / 6 * np.log(N / np.pi * np.sin(np.pi * L / N)) + c_prime


if __name__ == "__main__":
    # Q1: ground state energy and entanglement entropy
    results = {}
    # m = 10, g \in [0.5, 1, 1.5]
    for g in [0.5, 1, 1.5]:
        results.update(finite_system_algorithm(L=N, m_warmup=N, m_sweep=10, g=g))
    # g = 1, m \in [10, 20, 30]
    for m_sweep in [10, 20, 30]:
        if (N, m_sweep, 1) not in results:
            results.update(
                finite_system_algorithm(L=N, m_warmup=N, m_sweep=m_sweep, g=1)
            )
    # Q2: plot the entanglement entropy
    x_values = [i + 2 for i in range(N - 3)]
    plt.figure(figsize=(10, 6))
    for key, value in results.items():
        N, m, g = key
        plt.plot(x_values, value, label=f"(m, g)=({m}, {g})")

    plt.xlabel("L")
    plt.ylabel("Entanglement Entropy")
    plt.legend()
    plt.savefig("entanglement_entropy.png")

    # Q3: for m = 20, g = 1, fit S(L) using S(L) = c/6*log(N/\pi*\sin(\pi*L/N)) + c^{\prime}
    required_m = 20
    required_g = 1
    y_values = results[(N, required_m, required_g)]

    # Fit the entanglement entropy using curve S(L) = c/6*log(N/\pi*\sin(\pi*L/N)) + c^{\prime}
    popt, pcov = curve_fit(S, x_values, y_values)

    # Draw the original curve and the fitted curve
    plt.figure(figsize=(10, 6))
    plt.plot(
        x_values,
        y_values,
        label=f"(m, g)=({required_m}, {required_g})",
    )
    plt.plot(x_values, S(np.array(x_values), *popt), "-", label="Fitted Curve")

    plt.xlabel("L")
    plt.ylabel("Entanglement Entropy")
    plt.legend()
    plt.savefig("entanglement_entropy_fitted.png")

    print(f"Fitted parameters: c = {popt[0]}, c_prime = {popt[1]}")

import matplotlib.pyplot as plt
import numpy as np
import itertools


def Odot(x, y, n):
    # check the anti/commutation relation of two strings
    x1, x2 = x[:n], x[n:]
    y1, y2 = y[:n], y[n:]
    inner_product1 = sum(int(a) * int(b) for a, b in zip(x1, y2)) % 2
    inner_product2 = sum(int(a) * int(b) for a, b in zip(x2, y1)) % 2
    return (inner_product1 + inner_product2) % 2


def Count_Length(s, n):
    # count the length of a Pauli string
    # the length is defined as the number of non-identity Pauli operators
    x1, x2 = s[:n], s[n:]
    count1 = sum(int(bit) for bit in x1)
    count2 = sum(int(bit) for bit in x2)
    repeats = sum(int(a) * int(b) for a, b in zip(x1, x2))
    return count1 + count2 - repeats


def Generate_Pauli_Group(n):
    # generate all binary strings of length 2n(the whole Pauli group on n qubits)
    Pauli_group = ["".join(p)
                   for p in itertools.product("01", repeat=2 * n)]
    return Pauli_group


def Get_Generators(Z_1, Pauli_group):
    # get the generators of the Pauli group on n qubits
    # the process is as follows:
    # 1. find the first Z operator Z_1
    # 2. find the first X operator X_1, the first string obeying odot(X_1, Z_1) = 1
    # 3. find the second Z operator Z_2, the first string obeying odot(Z_2, Z_1) = 0 and odot(Z_2, X_1) = 0
    # ...
    generators = {}
    n = len(Z_1) / 2
    n = int(len(Z_1) / 2)
    generators["Z"] = [None] * n
    generators["X"] = [None] * n
    generators["Z"][0] = Z_1

    # process to get the generators, but we only have the "first" list of generators, depending on the order of the binary strings
    for candidate in Pauli_group:
        if Odot(Z_1, candidate, n) == 1:
            generators["X"][0] = candidate
            break
    for i in range(1, n):
        if generators["Z"][i] == None:
            for candidate in Pauli_group:
                Find_Z = True
                for j in range(i):
                    if (
                        Odot(candidate, generators["Z"][j], n) == 1
                        or Odot(candidate, generators["X"][j], n) == 1
                        or candidate == generators["Z"][j]
                        or candidate == generators["X"][j]
                        or candidate == "0" * 2 * n
                    ):
                        Find_Z = False
                        break
                if Find_Z:
                    generators["Z"][i] = candidate
                    break
        if generators["X"][i] == None:
            for candidate in Pauli_group:
                Find_X = True
                for j in range(i):
                    if (
                        Odot(candidate, generators["Z"][j], n) == 1
                        or Odot(candidate, generators["X"][j], n) == 1
                        or Odot(candidate, generators["Z"][i], n) == 0
                        or candidate == generators["Z"][j]
                        or candidate == generators["X"][j]
                        or candidate == "0" * 2 * n
                    ):
                        Find_X = False
                        break
                if Find_X:
                    generators["X"][i] = candidate
                    break
    return generators


def Generate_Pauli_Strings(generators, n):
    # generate all Pauli strings on n qubits, which has length of n initially
    blocks = {}
    blocks["Z"] = generators["Z"]
    blocks["X"] = generators["X"]
    blocks["Y"] = [
        bin(int(generators["Z"][i], 2) ^ int(
            generators["X"][i], 2))[2:].zfill(2 * n)
        for i in range(n)
    ]
    pauli_strings_repre = []
    pauli_strings = ["".join(p) for p in itertools.product("ZXY", repeat=n)]
    for pauli_string in pauli_strings:
        pauli_repre = "0" * 2 * n
        for i in range(n):
            pauli_repre = bin(int(pauli_repre, 2) ^ int(blocks[pauli_string[i]][i], 2))[
                2:
            ].zfill(2 * n)
        pauli_strings_repre.append(pauli_repre)
    return pauli_strings_repre


def Get_Pauli(P):
    # get the Pauli operators from the Pauli string
    n = int(len(P) / 2)
    X_part = P[:n]
    Z_part = P[n:]
    pauli = ""
    for i in range(n):
        if X_part[i] == "1" and Z_part[i] == "1":
            pauli += "Y"
        elif X_part[i] == "1":
            pauli += "X"
        elif Z_part[i] == "1":
            pauli += "Z"
        else:
            pauli += "I"
    return pauli


def Get_Minimun_Cases(n, generators=None, random_pick=True):
    # get the minimal average length for n qubits, store the results in a list
    Pauli_group = Generate_Pauli_Group(n)
    # store the minimal average length and the corresponding Pauli group
    min_length = n
    minimun_cases = []
    # if we don't have the generators, we need to randomly pick one
    if random_pick:
        for Z_1 in Pauli_group:
            if Z_1 == "0" * 2 * n:
                continue
            else:
                generators = Get_Generators(Z_1, Pauli_group)
                initial_n_length_Pauli = Generate_Pauli_Strings(generators, n)
                average_length = sum(Count_Length(pauli, n) for pauli in initial_n_length_Pauli) / len(
                    initial_n_length_Pauli
                )
                if average_length == min_length and average_length != n:
                    minimun_cases.append(
                        (average_length, initial_n_length_Pauli, generators))
                if average_length < min_length:
                    min_length = average_length
                    minimun_cases.clear()
                    minimun_cases.append(
                        (min_length, initial_n_length_Pauli, generators))
    # if we have the generators, we can directly calculate the average length(but may not the minimal one)
    elif random_pick == False and generators != None:
        initial_n_length_Pauli = Generate_Pauli_Strings(generators, n)
        average_length = sum(Count_Length(pauli, n) for pauli in initial_n_length_Pauli) / len(
            initial_n_length_Pauli
        )
        minimun_cases.append(
            (average_length, initial_n_length_Pauli, generators))
    return minimun_cases


def Print_Minimun_Cases(minimun_cases):
    # print the minimal average length and the corresponding Pauli group
    for i, case in enumerate(minimun_cases):
        # initialization
        length = case[0]
        Pauli_Strings = case[1]
        generator = case[2]
        num_of_qubit = len(generator["X"])
        print(f"=============The #{i+1} minimun case=============")
        print(f"The minimal average length is {length}")
        print("The generators transform as follows:")
        for i in range(num_of_qubit):
            initial_pauli = Get_Pauli(bin(1 << i)[2:].zfill(2 * num_of_qubit))
            final_pauli = Get_Pauli(generator["Z"][i])
            print(f"Z_{i+1}: {initial_pauli}------>{final_pauli}")
            initial_pauli = Get_Pauli(
                bin(1 << (i+num_of_qubit))[2:].zfill(2 * num_of_qubit))
            final_pauli = Get_Pauli(generator["X"][i])
            print(f"X_{i+1}: {initial_pauli}------>{final_pauli}")
        print("The Pauli strings are as follows:")
        for i, pauli in enumerate(Pauli_Strings):
            if i == len(Pauli_Strings) - 1:
                print(Get_Pauli(pauli))
            else:
                print(Get_Pauli(pauli), end="、")


def CNOT_tensor(n):
    # generate the tensor product of CNOT gates
    CNOT = [[1, 0, 0, 0], [0, 1, 0, 0],
            [0, 0, 0, 1], [0, 0, 1, 0]]
    CNOT_tensor = CNOT
    for i in range(n - 2):
        CNOT_tensor = np.kron(CNOT_tensor, CNOT)
    return CNOT_tensor


minimun_cases = Get_Minimun_Cases(3)
Print_Minimun_Cases(minimun_cases)


# Plotting the results
"""
plt.figure(figsize=(12, 6))

plt.subplot(1, 2, 1)
plt.plot(qubit_range, min_lengths, marker="o")
plt.title("Minimal Average Length vs. Number of Qubits")
plt.xlabel("Number of Qubits (n)")
plt.ylabel("Minimal Average Length")

plt.subplot(1, 2, 2)
plt.plot(qubit_range, min_ratios, marker="o")
plt.title("Minimal Ratio vs. Number of Qubits")
plt.xlabel("Number of Qubits (n)")
plt.ylabel("Minimal Ratio")

plt.tight_layout()
plt.savefig("minimal_average_length_and_ratio.png")
plt.show()
"""

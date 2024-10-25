import itertools


def odot(x, y, n):
    x1, x2 = x[:n], x[n:]
    y1, y2 = y[:n], y[n:]
    inner_product1 = sum(int(a) * int(b) for a, b in zip(x1, y2)) % 2
    inner_product2 = sum(int(a) * int(b) for a, b in zip(x2, y1)) % 2
    return (inner_product1 + inner_product2) % 2


def count_repeats(s, n):
    x1, x2 = s[:n], s[n:]
    count1 = sum(int(bit) for bit in x1)
    count2 = sum(int(bit) for bit in x2)
    repeats = sum(int(a) * int(b) for a, b in zip(x1, x2))
    return count1 + count2 - repeats


def generate_binary_strings(n):
    binary_strings = ["".join(p) for p in itertools.product("01", repeat=2 * n)]
    return sorted(binary_strings, key=lambda s: count_repeats(s, n), reverse=True)


def get_generators(Z_1, binary_strings):
    generators = {}
    n = len(Z_1) / 2
    n = int(len(Z_1) / 2)
    generators["Z"] = [None] * n
    generators["X"] = [None] * n
    generators["Z"][0] = Z_1
    for candidate in binary_strings:
        if odot(Z_1, candidate, n) == 1:
            generators["X"][0] = candidate
            break
    for i in range(1, n):
        if generators["Z"][i] == None:
            for candidate in binary_strings:
                Find_Z = True
                for j in range(i):
                    if (
                        odot(candidate, generators["Z"][j], n) == 1
                        or odot(candidate, generators["X"][j], n) == 1
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
            for candidate in binary_strings:
                Find_X = True
                for j in range(i):
                    if (
                        odot(candidate, generators["Z"][j], n) == 1
                        or odot(candidate, generators["X"][j], n) == 1
                        or odot(candidate, generators["Z"][i], n) == 0
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


def generate_pauli_strings(generators, n):
    blocks = {}
    blocks["Z"] = generators["Z"]
    blocks["X"] = generators["X"]
    blocks["Y"] = [generators["Z"][i] + generators["X"][i] for i in range(n)]
    blocks["Y"] = [
        bin(int(generators["Z"][i], 2) ^ int(generators["X"][i], 2))[2:].zfill(2 * n)
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


import matplotlib.pyplot as plt

min_Z_1_s = []
min_lengths = []
min_ratios = []
qubit_range = [1, 2, 3, 4, 5, 6]
for n in qubit_range:
    binary_strings = generate_binary_strings(n)
    print(f"for n = {n} qubits, we obtain the follows:")
    min_length = n
    min_ratio = 1
    min_Z_1 = None
    count = 0
    for Z_1 in binary_strings:
        if Z_1 == "0" * 2 * n:
            continue
        else:
            generators = get_generators(Z_1, binary_strings)
            paulis = generate_pauli_strings(generators, n)
            average_length = sum(count_repeats(pauli, n) for pauli in paulis) / len(
                paulis
            )
            ratio = average_length / n
            if average_length < min_length:
                min_length = average_length
                min_ratio = ratio
                min_Z_1 = Z_1
            print(
                f"for Z1 = {Z_1}, the average length is {average_length}, the ratio is {ratio}"
            )
    min_lengths.append(min_length)
    min_ratios.append(min_ratio)
    min_Z_1_s.append(min_Z_1)
    print(
        f"For n = {n} qubits, the minimal average length is {min_length}, the ratio is {min_ratio}, for Z1 = {min_Z_1}"
    )

print(f"The minimal Z1s are {min_Z_1_s} for qubits {qubit_range}")

# Plotting the results
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

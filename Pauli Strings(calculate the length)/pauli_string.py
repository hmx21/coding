import numpy as np
import itertools


def Binary_Addition_Mod2(x, y):
    # add two binary strings and take the result mod 2
    return "".join(str((int(a) + int(b)) % 2) for a, b in zip(x, y))


def Odot(x, y, n):
    # check the anti/commutation relation of two strings
    x1, x2 = x[:n], x[n:]
    y1, y2 = y[:n], y[n:]
    inner_product1 = sum(int(a) * int(b) for a, b in zip(x1, y2)) % 2
    inner_product2 = sum(int(a) * int(b) for a, b in zip(x2, y1)) % 2
    return (inner_product1 + inner_product2) % 2


def Generate_Pauli_Group(n):
    # generate all binary strings of length 2n(the whole Pauli group on n qubits)
    pauli_group = ["".join(p) for p in itertools.product("01", repeat=2 * n)]
    return pauli_group


def Get_Minimun_Cases(n):
    # get the minimal average length for n qubits(sweeping the whole Pauli group to get the generators), store the results in a list
    Pauli_group = Generate_Pauli_Group(n)
    # store the minimal average length and the corresponding Pauli group
    min_length = n
    minimun_cases = []
    # if we don't have the generators, we need to randomly pick one
    for Z_1 in Pauli_group:
        if Z_1 == "0" * 2 * n:
            continue
        else:
            generators = Get_Random_Generators(Z_1, Pauli_group)
            clifford = Clifford_Operator(generators)
            average_length = clifford.Get_Average_Length(n)
            if average_length == min_length and average_length != n:
                minimun_cases.append((average_length, clifford))
            if average_length < min_length:
                min_length = average_length
                minimun_cases.clear()
                minimun_cases.append((min_length, clifford))
    return minimun_cases


def Get_Random_Generators(Z_1, pauli_group):
    # get the generators of the Pauli group on n qubits
    # the process is as follows:
    # 1. find the first Z operator Z_1
    # 2. find the first X operator X_1, the first string obeying odot(X_1, Z_1) = 1
    # 3. find the second Z operator Z_2, the first string obeying odot(Z_2, Z_1) = 0 and odot(Z_2, X_1) = 0
    # ...
    n = int(len(Z_1) / 2)
    generators = {}
    generators["Z"] = [None] * n
    generators["X"] = [None] * n
    generators["Z"][0] = Z_1

    # process to get the generators, but we only have the "first" list of generators, depending on the order of the binary strings
    for candidate in pauli_group:
        if Odot(Z_1, candidate, n) == 1:
            generators["X"][0] = candidate
            break

    for i in range(1, n):
        if generators["Z"][i] == None:
            for candidate in pauli_group:
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
            for candidate in pauli_group:
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


def H_Effect(sympletic_form, qubit):
    # apply a Hadamard gate to a Pauli operator
    n = int(len(sympletic_form) / 2)
    X_part = list(sympletic_form[:n])
    Z_part = list(sympletic_form[n:])
    if int(X_part[qubit]) + int(Z_part[qubit]) == 1:
        X_part[qubit], Z_part[qubit] = Z_part[qubit], X_part[qubit]
    return "".join(X_part) + "".join(Z_part)


def CNOT_Effect(sympletic_form, control, target):
    # apply a CNOT gate to a Pauli operator
    n = int(len(sympletic_form) / 2)
    X_part = list(sympletic_form[:n])
    Z_part = list(sympletic_form[n:])
    if Z_part[target] == "1":
        Z_part[control] = str((int(Z_part[control]) + 1) % 2)
    if X_part[control] == "1":
        X_part[target] = str((int(X_part[target]) + 1) % 2)
    return "".join(X_part) + "".join(Z_part)


def S_Effect(sympletic_form, qubit):
    # apply a S gate to a Pauli operator
    n = int(len(sympletic_form) / 2)
    X_part = list(sympletic_form[:n])
    Z_part = list(sympletic_form[n:])
    if X_part[qubit] == "1":
        Z_part[qubit] = str((int(Z_part[qubit]) + 1) % 2)
    return "".join(X_part) + "".join(Z_part)


def Gate_Effect(gate, generators):
    # apply a gate to the generators of the Clifford group
    if gate[0] == "H":
        qubit = gate[1]
        for i, x_gate in enumerate(generators["X"]):
            generators["X"][i] = H_Effect(x_gate, qubit)
        for i, z_gate in enumerate(generators["Z"]):
            generators["Z"][i] = H_Effect(z_gate, qubit)
    elif gate[0] == "CNOT":
        control, target = gate[1]
        for i, x_gate in enumerate(generators["X"]):
            generators["X"][i] = CNOT_Effect(x_gate, control, target)
        for i, z_gate in enumerate(generators["Z"]):
            generators["Z"][i] = CNOT_Effect(z_gate, control, target)
    elif gate[0] == "S":
        qubit = gate[1]
        for i, x_gate in enumerate(generators["X"]):
            generators["X"][i] = S_Effect(x_gate, qubit)
        for i, z_gate in enumerate(generators["Z"]):
            generators["Z"][i] = S_Effect(z_gate, qubit)
    return generators


class Pauli_Operator(object):
    def __init__(self, sympletic_form):
        # sympletic_form is a binary string of length 2n
        self.sympletic_form = sympletic_form
        # get the Pauli string from the sympletic form, example: "1101" -> "XY"
        self.pauli_string = self.Get_Pauli()
        # count the length of a Pauli string
        self.length = self.Count_Length()

    def Get_Pauli(self):
        # get the Pauli operators from the sympletic form
        n = int(len(self.sympletic_form) / 2)
        X_part = self.sympletic_form[:n]
        Z_part = self.sympletic_form[n:]
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

    def Count_Length(self):
        # count the length of a Pauli string
        # the length is defined as the number of non-identity Pauli operators
        n = int(len(self.sympletic_form) / 2)
        x1, x2 = self.sympletic_form[:n], self.sympletic_form[n:]
        count1 = sum(int(bit) for bit in x1)
        count2 = sum(int(bit) for bit in x2)
        repeats = sum(int(a) * int(b) for a, b in zip(x1, x2))
        return count1 + count2 - repeats


class Clifford_Operator(object):
    def __init__(self, generators, gate_list=None, n=None):
        # "generators" is a dictionary of the generators of the Clifford group
        # the effect of a Clifford operator is to permute the Pauli strings
        # once we know the transformation of the generators, we can get the transformation of all the Pauli strings
        # thus we have all known the Clifford operator
        self.generators = generators
        self.number_of_qubits = n if n != None else int(
            len(self.generators["Z"]))
        self.initial_generators = self.Get_Initial_Generators()
        # store the effect of the Clifford operator on all the Pauli strings
        self.all_pauli_strings = (
            self.Get_All_Pauli_Strings() if self.generators != None else None
        )
        # the gate is another way to represent the Clifford operator
        # we just write the Clifford operator in terms of H, CNOT, S
        # the gate_list is a list of tuples, each tuple represent a gate and the qubits it acts on
        # the first element of the tuple is the gate, the second element is the qubits
        # the position of the tuple in the list is the order of the gate
        # example: [("H", 0), ("CNOT", (0, 1)), ("S", 1)] means S_2 \otimes CNOT_{12} \otimes H_1
        # in this example, the control qubit of CNOT is 1, the target qubit is 2
        # note that the qubit position is labeled from 0 in the code but from 1 in reality
        self.gate_list = gate_list
        if self.gate_list != None:
            self.generators = self.Get_Generators_From_Gate()
            self.all_pauli_strings = self.Get_All_Pauli_Strings()

    def Get_Initial_Generators(self):
        # get the initial generators of the Clifford group
        if self.generators == None:
            n = self.number_of_qubits
        else:
            n = int(len(self.generators["Z"]))
        initial_generators = {"Z": [], "X": []}
        for i in range(n):
            initial_generators["Z"].append(
                bin((1 << (2 * n - 1)) >> (i + n))[2:].zfill(2 * n)
            )

            initial_generators["X"].append(
                bin((1 << (2 * n - 1)) >> i)[2:].zfill(2 * n)
            )

        return initial_generators

    def Get_All_Pauli_Strings(self):
        # get the effect of the Clifford operator on all the Pauli strings
        n = int(len(self.generators["Z"]))
        all_pauli_strings = {}
        pauli_group_sympletic = Generate_Pauli_Group(n)
        for p in pauli_group_sympletic:
            pauli = Pauli_Operator(p)
            temp = "0" * 2 * n
            X_part = p[:n]
            Z_part = p[n:]
            for i in range(n):
                if X_part[i] == "1":
                    temp = Binary_Addition_Mod2(temp, self.generators["X"][i])
                if Z_part[i] == "1":
                    temp = Binary_Addition_Mod2(temp, self.generators["Z"][i])
            new_pauli = Pauli_Operator(temp)
            all_pauli_strings[pauli] = new_pauli
        return all_pauli_strings

    def Print_Clifford_Map(self, initial_length=None):
        # print the effect of the Clifford operator on all the Pauli strings
        if initial_length != None:
            print(
                f"For the Pauli operators which has initial length of {initial_length}, one can derive the following transformation:"
            )
            for key, value in self.all_pauli_strings.items():
                if key.length == initial_length:
                    print(key.pauli_string, "----->", value.pauli_string)
        else:
            print("The transformation of all the Pauli strings are as follows:")
            for key, value in self.all_pauli_strings.items():
                print(key.pauli_string, "----->", value.pauli_string)

    def Get_Average_Length(self, initial_length):
        # get the average length of the Clifford operator
        total_length = 0
        for key, value in self.all_pauli_strings.items():
            if key.length == initial_length:
                total_length += value.length
        return total_length / (3 ** len(self.generators["Z"]))

    def Get_Generators_From_Gate(self):
        # get the generators of the Clifford group from the gate
        generators = self.initial_generators
        for gate in self.gate_list:
            generators = Gate_Effect(gate, generators)
        return generators

    def Plot_Distribution(self, initial_length):
        # plot the distribution of the length of the Pauli strings
        import matplotlib.pyplot as plt

        length_list = []
        for key, value in self.all_pauli_strings.items():
            if key.length == initial_length:
                length_list.append(value.length)
        plt.hist(length_list, bins=range(2 * self.number_of_qubits + 1))
        plt.xlabel("Length of the Pauli strings")
        plt.ylabel("Frequency")
        plt.title("Distribution of the length of the Pauli strings")
        plt.show()


def Get_CNOT_Tensor_Product(n):
    # get the tensor product of CNOT gate on n qubits
    CNOT_tensor = []
    if n % 2 == 0:
        for i in range(int(n / 2)):
            CNOT_tensor.append(("CNOT", (2 * i, 2 * i + 1)))
    else:
        for i in range(int(n // 2)):
            CNOT_tensor.append(("CNOT", (2 * i, 2 * i + 1)))
        CNOT_tensor.append(("CNOT", (n - 2, n - 1)))
    return CNOT_tensor


def Get_CNOT_Cyclic(n):
    # get the cyclic CNOT gate on n qubits
    CNOT_tensor = []
    for i in range(n):
        CNOT_tensor.append(("CNOT", (i, (i + 1) % n)))
    return CNOT_tensor


if __name__ == "__main__":

    # Result 0
    Hadamard = Clifford_Operator(None, [("H", 0)], 1)
    print("For Hadamard gate, the corresponding transformation is:")
    Hadamard.Print_Clifford_Map()
    Phase = Clifford_Operator(None, [("S", 0)], 1)
    print("For Phase gate, the corresponding transformation is:")
    Phase.Print_Clifford_Map()
    CNOT = Clifford_Operator(None, [("CNOT", (0, 1))], 2)
    print("For CNOT gate, the corresponding transformation is:")
    CNOT.Print_Clifford_Map()

    # Result 1
    for n in range(2, 8):
        CNOT_tensor = Get_CNOT_Tensor_Product(n)
        CNOT_for_result_1 = Clifford_Operator(None, CNOT_tensor, n)
        print(
            f"For {n} qubits, the average length of the transformation by {CNOT_tensor} is {CNOT_for_result_1.Get_Average_Length(n)}"
        )

    # Result 2
    # CNOT_13\otimes CNOT_23\otimes CNOT_12
    gate_list = [
        ("CNOT", (0, 1)),
        ("CNOT", (1, 2)),
        ("CNOT", (2, 0)),
    ]
    number_of_qubits = 3
    CNOT_13_23_12 = Clifford_Operator(None, gate_list, number_of_qubits)
    CNOT_13_23_12.Print_Clifford_Map(number_of_qubits)
    print(
        f"The average length of the transformation by {gate_list} is {CNOT_13_23_12.Get_Average_Length(number_of_qubits)}"
    )
    CNOT_13_23_12.Plot_Distribution(number_of_qubits)

    # More for Result 2
    for i in range(3, 8):
        CNOT_cyclic = Get_CNOT_Cyclic(i)
        CNOT_cyclic_operator = Clifford_Operator(None, CNOT_cyclic, i)
        print(
            f"For {i} qubits, the average length of the transformation by {CNOT_cyclic} is {CNOT_cyclic_operator.Get_Average_Length(i)}"
        )

    # the fake!!! minimun cases
    # min_cases = Get_Minimun_Cases(number_of_qubits)
    # for case in min_cases:
    # print("The minimal average length is", case[0])

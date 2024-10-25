import numpy as np
import numpy.linalg as LA
import scipy.sparse.linalg as LAs
import Sub180221 as Sub
import math, copy


# for H=-\sum_j(Sz_j*Sy_{j+1}Sz_{j+2}+Sz_j*Sz_{j+1}+Sx_j), generate the corresponding MPO
#    M =
#    s0	-sz	0	-sx
#    0	0	sy  sz
#    0	0	0	sz
#    0	0	0	s0
def GetMpo_Obc(Dp):
    # Spin operators
    S0, Sp, Sm, Sz, Sx, Sy = Sub.SpinOper(Dp)

    # generate the MPO
    Dmpo = 4
    Mpo = np.zeros((Dmpo, Dp, Dmpo, Dp), dtype=complex)

    Mpo[0, :, 0, :] = S0
    Mpo[0, :, 1, :] = -2 * Sz
    Mpo[0, :, 3, :] = -2 * Sx
    Mpo[1, :, 2, :] = 2 * Sy
    Mpo[1, :, 3, :] = 2 * Sz
    Mpo[2, :, 3, :] = 2 * Sz
    Mpo[3, :, 3, :] = S0

    return Mpo


# Tensor initialization
def InitMps(Ns, Dp, Ds):
    T = [None] * Ns
    for i in range(Ns):
        Dl = min(Dp**i, Dp ** (Ns - i), Ds)
        Dr = min(Dp ** (i + 1), Dp ** (Ns - 1 - i), Ds)
        T[i] = np.random.rand(Dl, Dp, Dr)
    U = np.eye(np.shape(T[-1])[-1])
    # LQ decompositions start from the right end
    for i in range(Ns - 1, 0, -1):
        U, T[i] = Sub.Mps_LQP(T[i], U)
    return T


# Environment initialization
def InitH(Mpo, T):
    Ns = len(T)
    Dmpo = np.shape(Mpo)[0]
    HL = [None] * Ns
    HR = [None] * Ns
    HL[0] = np.zeros((1, Dmpo, 1))
    HL[0][0, 0, 0] = 1.0
    HR[-1] = np.zeros((1, Dmpo, 1))
    HR[-1][0, -1, 0] = 1.0
    for i in range(Ns - 1, 0, -1):
        HR[i - 1] = Sub.NCon(
            [HR[i], T[i], Mpo, np.conj(T[i])],
            [[1, 3, 5], [-1, 2, 1], [-2, 2, 3, 4], [-3, 4, 5]],
        )
    return HL, HR


# Optimize each T[i]
def OptTSite(Mpo, HL, HR, T, Method=0):
    DT = np.shape(T)
    Dl = np.prod(DT)
    # direct method
    if Method == 0:
        A = Sub.NCon(
            [HL, Mpo, HR],
            [[-1, 1, -4], [1, -5, 2, -2], [-6, 2, -3]],
        )
        # effective Hamiltonian
        A = Sub.Group(A, [[0, 1, 2], [3, 4, 5]])
        Eig, V = LAs.eigsh(A, k=1, which="SA")
        T = np.reshape(V, DT)
        Eig = np.real(Eig[0])
    # iterative method
    if Method == 1:

        def UpdateV(V):
            V = np.reshape(V, DT)
            V = Sub.NCon(
                [HL, V, Mpo, HR],
                [[-1, 3, 1], [1, 2, 4], [3, 2, 5, -2], [4, 5, -3]],
            )
            V = np.reshape(V, [Dl])
            return V

        V0 = np.reshape(T, [Dl])
        MV = LAs.LinearOperator((Dl, Dl), matvec=UpdateV)
        Eig, V = LAs.eigsh(MV, k=1, which="SA", v0=V0)
        T = np.reshape(V, DT)
        Eig = np.real(Eig[0])
    # optimize two sites for one step
    if Method == 2:
        A = Sub.NCon(
            [HL, Mpo, Mpo, HR],
            [[-1, 1, -5], [1, -6, 2, -2], [2, -7, 3, -3], [-8, 3, -4]],
        )
        A = Sub.Group(A, [[0, 1, 2, 3], [4, 5, 6, 7]])
        Eig, V = LAs.eigsh(A, k=1, which="SA")
        # two site MPS, reshape into a matrix, in order to do SVD
        T = np.reshape(V, [DT[0] * DT[1], DT[2] * DT[3]])
        Eig = np.real(Eig[0])

    return T, Eig


# Sweep back and forth to optimize each tensor
def OptT(Mpo, HL, HR, T):
    Ns = len(T)
    Eng0 = np.zeros(Ns)
    Eng1 = np.zeros(Ns)
    for r in range(100):
        # sweep from left to right
        for i in range(Ns - 1):
            T[i], Eng1[i] = OptTSite(Mpo, HL[i], HR[i], T[i], Method=0)
            T[i], U = Sub.Mps_QR0P(T[i])
            HL[i + 1] = Sub.NCon(
                [HL[i], np.conj(T[i]), Mpo, T[i]],
                [[1, 3, 5], [1, 2, -1], [3, 4, -2, 2], [5, 4, -3]],
            )
            T[i + 1] = np.tensordot(U, T[i + 1], (1, 0))
        # sweep from right to left
        for i in range(Ns - 1, 0, -1):
            T[i], Eng1[i] = OptTSite(Mpo, HL[i], HR[i], T[i], Method=0)
            U, T[i] = Sub.Mps_LQ0P(T[i])
            HR[i - 1] = Sub.NCon(
                [HR[i], T[i], Mpo, np.conj(T[i])],
                [[1, 3, 5], [-1, 2, 1], [-2, 2, 3, 4], [-3, 4, 5]],
            )
            T[i - 1] = np.tensordot(T[i - 1], U, (2, 0))
        # stop if the energy converges
        if abs(Eng1[1] - Eng0[1]) < 1.0e-7:
            break
        Eng0 = copy.copy(Eng1)
    return T, Eng1 / float(Ns)


# Sweep back and forth to optimize each tensor(but two site optimization for one step)
def OptT_SVD(Mpo, HL, HR, T):
    Ns = len(T)
    Eng0 = np.zeros(Ns)
    Eng1 = np.zeros(Ns)
    for r in range(100):
        # sweep from left to right
        for i in range(Ns - 1):
            temp, Eng1[i] = OptTSite(
                Mpo, HL[i], HR[i + 1], np.tensordot(T[i], T[i + 1], (2, 0)), Method=2
            )
            # the bond dimension of mps
            bond_dim = np.shape(T[i])[2]
            # SVD decomposition to get the canonical form
            U, S, V = LA.svd(
                temp,
                full_matrices=False,
            )
            U = U[:, 0:bond_dim]
            S = S[0:bond_dim]
            V = V[0:bond_dim, :]
            # update the tensors via U, S, V(pay attention to the direction)
            T[i] = np.reshape(U, np.shape(T[i]))
            T[i + 1] = np.reshape(np.dot(np.diag(S), V), np.shape(T[i + 1]))
            HL[i + 1] = Sub.NCon(
                [HL[i], np.conj(T[i]), Mpo, T[i]],
                [[1, 3, 5], [1, 2, -1], [3, 4, -2, 2], [5, 4, -3]],
            )
        # sweep from right to left
        for i in range(Ns - 1, 0, -1):
            temp, Eng1[i] = OptTSite(
                Mpo, HL[i - 1], HR[i], np.tensordot(T[i - 1], T[i], (2, 0)), Method=2
            )
            # the bond dimension of mps
            bond_dim = np.shape(T[i])[0]
            # SVD decomposition to get the canonical form
            U, S, V = LA.svd(
                temp,
                full_matrices=False,
            )
            U = U[:, 0:bond_dim]
            S = S[0:bond_dim]
            V = V[0:bond_dim, :]
            T[i] = np.reshape(V, np.shape(T[i]))
            T[i - 1] = np.reshape(np.dot(U, np.diag(S)), np.shape(T[i - 1]))
            HR[i - 1] = Sub.NCon(
                [HR[i], T[i], Mpo, np.conj(T[i])],
                [[1, 3, 5], [-1, 2, 1], [-2, 2, 3, 4], [-3, 4, 5]],
            )
        # stop if the energy converges
        if abs(Eng1[1] - Eng0[1]) < 1.0e-7:
            break
        Eng0 = copy.copy(Eng1)
    return T, Eng1 / float(Ns)


# To calculate the magnetization per site
# auxiliary functions
def FlipBit(i, n):
    return i ^ (1 << n)


def ReadBit(i, n):
    return (i & (1 << n)) >> n


# calculate the magnetization in z direction
def MagnetZ(T):
    Ns = len(T)
    Sz = 2 * Sub.SpinOper(np.shape(T[0])[1])[3]
    Mz = np.zeros(Ns)
    # because the T[0] is not a canonical form, we need to calculate it separately.
    T_mixterm = Sub.NCon([T[0], np.conj(T[0])], [[1, 2, -1], [1, 2, -2]])
    for i in range(Ns):
        if i == 0:
            Mz[i] += np.real(
                Sub.NCon([T[i], Sz, np.conj(T[i])], [[1, 2, 4], [2, 3], [1, 3, 4]])
            )
        else:
            Mz[i] += np.real(
                Sub.NCon(
                    [T_mixterm, T[i], Sz, np.conj(T[i])],
                    [[1, 2], [1, 3, 5], [3, 4], [2, 4, 5]],
                )
            )
            # update the T_mixterm
            T_mixterm = Sub.NCon(
                [T_mixterm, T[i], np.conj(T[i])],
                [
                    [
                        1,
                        2,
                    ],
                    [1, 3, -1],
                    [2, 3, -2],
                ],
            )
    return Mz


# calculate the magnetization in x direction
def MagnetX(T):
    Ns = len(T)
    Sx = 2 * Sub.SpinOper(np.shape(T[0])[1])[4]
    Mx = np.zeros(Ns)
    # because the T[0] is not a canonical form, we need to calculate it separately.
    T_mixterm = Sub.NCon([T[0], np.conj(T[0])], [[1, 2, -1], [1, 2, -2]])
    for i in range(Ns):
        if i == 0:
            Mx[i] += np.real(
                Sub.NCon([T[i], Sx, np.conj(T[i])], [[1, 2, 4], [2, 3], [1, 3, 4]])
            )
        else:
            Mx[i] += np.real(
                Sub.NCon(
                    [T_mixterm, T[i], Sx, np.conj(T[i])],
                    [[1, 2], [1, 3, 5], [3, 4], [2, 4, 5]],
                )
            )
            # update the T_mixterm
            T_mixterm = Sub.NCon(
                [T_mixterm, T[i], np.conj(T[i])],
                [
                    [
                        1,
                        2,
                    ],
                    [1, 3, -1],
                    [2, 3, -2],
                ],
            )
    return Mx


# Exact diagonalization for this model
def Get_Hamiltonian(Ns):
    H = np.zeros((2**Ns, 2**Ns), dtype=complex)
    for row in range(2**Ns):
        # Sx term
        for site in range(Ns):
            column = FlipBit(row, site)
            H[row][column] += -1
        # Sz_j*Sz_{j+1} term
        for site in range(Ns - 1):
            H[row][row] += (
                -1 * 2 * (ReadBit(row, site) - 0.5) * 2 * (ReadBit(row, site + 1) - 0.5)
            )
        # Sz_j*Sy_{j+1}Sz_{j+2} term
        for site in range(Ns - 2):
            column = FlipBit(row, site + 1)
            if ReadBit(row, site + 1) == 1:
                H[row][column] += (
                    -1.0j
                    * 2
                    * (ReadBit(row, site) - 0.5)
                    * 2
                    * (ReadBit(row, site + 2) - 0.5)
                )
            else:
                H[row][column] += (
                    1.0j
                    * 2
                    * (ReadBit(row, site) - 0.5)
                    * 2
                    * (ReadBit(row, site + 2) - 0.5)
                )
    return H


# calculate the magnetization in z direction via exact diagonalization
def MagnetZ_ED(V, Ns):
    Mz = np.zeros(Ns, dtype=complex)
    for site in range(Ns):
        for state in range(2**Ns):
            Mz[site] += (np.abs(V[state]) ** 2 * 2 * (ReadBit(state, site) - 0.5))[0]
    Mz = np.real(Mz)
    return Mz


# calculate the magnetization in x direction via exact diagonalization
def MagnetX_ED(V, Ns):
    Mx = np.zeros(Ns, dtype=complex)
    for site in range(Ns):
        for state in range(2**Ns):
            Mx[site] += (np.conj(V[state]) * V[FlipBit(state, site)])[0]
    Mx = np.real(Mx)
    return Mx


# Question: calculate the energy per site and  magnetization per site
if __name__ == "__main__":
    Ns = 10
    Dp = 2
    Ds = [4, 6]

    # MPS variational optimization
    for Ds in Ds:
        print(f"For the case of Ns = {Ns}, Ds = {Ds}, we can get:")
        # 2 sites per step
        Mpo = GetMpo_Obc(Dp)
        T = InitMps(Ns, Dp, Ds)
        HL, HR = InitH(Mpo, T)
        T, Energy = OptT_SVD(Mpo, HL, HR, T)
        Mx = MagnetX(T)
        Mz = MagnetZ(T)
        print(
            "===========================2-site optimization==========================="
        )
        print(f"Energy per site: {Energy}")
        print(f"Magnetization in x direction per site: {Mx}")
        print(f"Magnetization in z direction per site: {Mz}")
        # 1 site per step
        Mpo = GetMpo_Obc(Dp)
        T = InitMps(Ns, Dp, Ds)
        HL, HR = InitH(Mpo, T)
        T, Energy = OptT(Mpo, HL, HR, T)
        Mx = MagnetX(T)
        Mz = MagnetZ(T)
        print(
            "===========================1-site optimization==========================="
        )
        print(f"Energy per site: {Energy}")
        print(f"Magnetization in x direction per site: {Mx}")
        print(f"Magnetization in z direction per site: {Mz}")

    # exact diagonalization
    H = Get_Hamiltonian(Ns)
    E, GS = LAs.eigsh(H, k=1, which="SA")
    Mx = MagnetX_ED(GS, Ns)
    Mz = MagnetZ_ED(GS, Ns)
    print("===========================Exact diagonalization===========================")
    print(f"Energy per site: {E / Ns}")
    print(f"Magnetization in x direction per site: {Mx}")
    print(f"Magnetization in z direction per site: {Mz}")

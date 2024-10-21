import numpy as np
import scipy.sparse.linalg as LA

# Exact diagonalization of a 3x3 lattice model


# Bit manipulation functions
def FlipBit(i, n):
    return i ^ (1 << n)


def ReadBit(i, n):
    return (i & (1 << n)) >> n


# Model parameters
Jx = 0.3
Jy = 1.4
hz = 0.5

N = 9
n = 2**N

# Construct the Hamiltonian
H_Neighbor = [(0, 3), (3, 6), (1, 4), (4, 7), (2, 5), (5, 8)]
V_Neighbor = [(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3), (6, 7), (7, 8), (8, 6)]

H = np.zeros([n, n])
for row in range(n):
    for i, j in H_Neighbor:
        column = FlipBit(FlipBit(row, j), i)
        H[row][column] += -Jx
    for i, j in V_Neighbor:
        column = FlipBit(FlipBit(row, j), i)
        symbol = ReadBit(row, i) + ReadBit(row, j)
        H[row][column] += -Jy if symbol == 1 else Jy
    for i in range(N):
        symbol = ReadBit(row, i)
        H[row][row] += hz if symbol == 1 else -hz

E_exact, W = LA.eigsh(H, 20, which="SA")
print(E_exact)

# Matrix product operator (MPO) representation of the Hamiltonian

# Define the matrix blocks
sigma_x = [[0, 1], [1, 0]]
sigma_y = [[0, -1.0j], [1.0j, 0]]
O = [[0, 0], [0, 0]]
I = [[1, 0], [0, 1]]
_Jx_sigma_x = [[0, -Jx], [-Jx, 0]]
_Jy_sigma_y = [[0, -1.0j * -Jy], [1.0j * -Jy, 0]]
hz_sigma_z = [[hz, 0], [0, -hz]]


def OP_prod(A, B):
    return np.kron(A, B)


def MPO_prod(A, B, D):
    C = [[None for _ in range(D)] for _ in range(D)]
    for i in range(D):
        for j in range(D):
            for k in range(D):
                if C[i][j] is None:
                    C[i][j] = OP_prod(A[i][k], B[k][j])
                else:
                    C[i][j] = C[i][j] + OP_prod(A[i][k], B[k][j])
    return C


# Construct the MPO

D = 7

MPO = [None for _ in range(N)]
for i in range(N):
    MPO[i] = [[None for _ in range(D)] for _ in range(D)]
    for j in range(D):
        for k in range(D):
            MPO[i][j][k] = O

for i in range(N):
    MPO[i][0][0] = I
    MPO[i][0][6] = hz_sigma_z
    MPO[i][6][6] = I

for i in [0, 1, 3, 4, 6, 7]:
    MPO[i][0][1] = _Jx_sigma_x

for i in range(6):
    MPO[i][0][2] = _Jy_sigma_y

for i in [1, 4, 7]:
    MPO[i][1][3] = I

for i in [1, 2, 4, 5, 7, 8]:
    MPO[i][1][6] = sigma_x

for i in range(1, 7):
    MPO[i][2][4] = I

for i in [2, 5, 8]:
    MPO[i][3][6] = sigma_x

for i in range(2, 8):
    MPO[i][4][5] = I

for i in range(3, N):
    MPO[i][5][6] = sigma_y


# Construct the Hamiltonian using the MPO

H_MPO = MPO[0]

for i in range(1, N):
    H_MPO = MPO_prod(H_MPO, MPO[i], D)

H_MPO = np.real(H_MPO[0][D - 1])

E_MPO, W_MPO = LA.eigsh(H_MPO, 20, which="SA")
print(E_MPO)


# Compare the exact diagonalization result with the MPO result
print(np.allclose(E_exact, E_MPO))

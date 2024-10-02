# from scipy import sparse


# set nth bit to 1
def SetBit(i, n):
    return i | (1 << n)


# clear nth bit to 0
def ClearBit(i, n):
    return i & ~(1 << n)


# flip nth bit
def FlipBit(i, n):
    return i ^ (1 << n)


# read nth bit
def ReadBit(i, n):
    return (i & (1 << n)) >> n


# count how many 1 bits in i
def PopCntBit(i):
    return bin(i).count("1")


# pick up n bits from kth bit
def PickBits(i, k, n):
    return (i & ((2**n - 1) << k)) >> k


# circular bit shift left
def RotLBit(i, L, n):
    return (PickBits(i, 0, L - n) << n) + (i >> (L - n))


# circular bit shift right
def RotRBit(i, L, n):
    return (i >> n) + ((PickBits(i, 0, n) << L - n))


# Sample code of AF Heisenberg model on a 1D chain with open boundary condition
# Ns=4

# dimension of the Hamiltonian matrix
# Nl=2**Ns

# get the list of hopping terms(only nearest neighbor hopping)
# def GetHopList(Ns):
#     HopList = []
#     for i in range(Ns-1):
#         HopList.append([i,i+1])
#     return HopList

# HopList = GetHopList(Ns)

# get the index of the matrix element, HI and HJ are the row and column index of the matrix element, HV is the value of the matrix element
# HI=[]
# HJ=[]
# HV=[]

# get the off-diagonal matrix element of the Heisenberg model
# for i0 in range(Nl):

#     # ih is the index of the hopping term, something like (1,2)
#     for ih in range(len(HopList)):
#         Pos0=HopList[ih][0]
#         Pos1=HopList[ih][1]

#         # the hopping term exists only when the two sites have different spins
#         if ReadBit(i0,Pos0)!=ReadBit(i0,Pos1):

#             # i1 is the row index of the matrix element, i0 is the column index of the matrix element
#             i1=FlipBit(i0,Pos0)
#             i1=FlipBit(i1,Pos1)
#             HI.append(i0)
#             HJ.append(i1)
#             HV.append(0.5)

#         # get the diagonal matrix element of the Heisenberg model
#     HI.append(i0)
#     HJ.append(i0)
#     HV.append((ReadBit(i0,Pos0)-0.5)*(ReadBit(i0,Pos1)-0.5))

# print the Hamiltonian matrix
# Hamr = sparse.coo_matrix((HV, (HI, HJ)), shape=(2**Ns, 2**Ns)).tocsc()
# eigval, eigvec = sparse.linalg.eigsh(Hamr, k=1, which='SA', return_eigenvectors=True)

# with the periodic boundary condition
# def GetHopListPBC(Ns):
#     HopList = []
#     for i in range(Ns):
#         HopList.append([i, (i+1) % Ns])
#     return HopList

# HopListPBC = GetHopListPBC(Ns)

# HIPBC = []
# HJPBC = []
# HVPBC = []

# for i0 in range(Nl):
#     for ih in range(len(HopListPBC)):
#         Pos0 = HopListPBC[ih][0]
#         Pos1 = HopListPBC[ih][1]
#         if ReadBit(i0, Pos0) != ReadBit(i0, Pos1):
#             i1 = FlipBit(i0, Pos0)
#             i1 = FlipBit(i1, Pos1)
#             HIPBC.append(i0)
#             HJPBC.append(i1)
#             HVPBC.append(0.5)
#     HIPBC.append(i0)
#     HJPBC.append(i0)
#     HVPBC.append((ReadBit(i0, Pos0) - 0.5) * (ReadBit(i0, Pos1) - 0.5))

# HamrPBC = sparse.coo_matrix((HVPBC, (HIPBC, HJPBC)), shape=(2**Ns, 2**Ns)).tocsc()
# eigval2, eigvec2 = sparse.linalg.eigsh(HamrPBC, k=1, which='SA', return_eigenvectors=True)

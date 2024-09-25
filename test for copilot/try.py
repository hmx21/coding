# Define Pauli matrices
import numpy as np

Sx = np.array([[0, 1], [1, 0]], dtype=complex)
Sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
Sz = np.array([[1, 0], [0, -1]], dtype=complex)

# Define spin-1/2 operators for three spins
# Define identity matrix for single spin
I = np.eye(2, dtype=complex)

# Define spin-1/2 operators for three spins in an 8-dimensional Hilbert space
S1 = [np.kron(np.kron(S, I), I) for S in [Sx, Sy, Sz]]
S2 = [np.kron(np.kron(I, S), I) for S in [Sx, Sy, Sz]]
S3 = [np.kron(np.kron(I, I), S) for S in [Sx, Sy, Sz]]

# Hamiltonian H = S1·S2 + S2·S3
H = sum(np.dot(S1[i], S2[i]) for i in range(3)) + sum(np.dot(S2[i], S3[i]) for i in range(3))

# Check the shape of the Hamiltonian
print("Shape of the Hamiltonian:", H.shape)

# Calculate eigenvalues
eigenvalues, _ = np.linalg.eigh(H)
print("Eigenvalues of the Hamiltonian:", eigenvalues)


# Check if a number is prime
def is_prime(num):
    if num < 2:
        return False
    for i in range(2, int(np.sqrt(num)) + 1):
        if num % i == 0:
            return False
    return True

# give me an array of prime numbers
def prime_array(n):
    primes = []
    i = 2
    while len(primes) < n:
        if is_prime(i):
            primes.append(i)
        i += 1
    return primes

# print the first 10 prime numbers
print(prime_array(10))
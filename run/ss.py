import numpy as np

def arnoldi_iteration(A, b, k):
    """
    Performs k steps of the Arnoldi iteration to compute an orthonormal basis
    for the k-dimensional Krylov subspace and an upper Hessenberg matrix H.

    Parameters:
    A (np.array): The square matrix for which to compute the Krylov subspace.
    b (np.array): The initial vector (must be of the same dimension as A's rows).
    k (int): The number of Arnoldi iterations to perform.

    Returns:
    Q (np.array): An m x (k+1) array where columns are the orthonormal basis vectors.
    H (np.array): A (k+1) x k upper Hessenberg matrix.
    """
    m = A.shape[0]
    if b.shape[0] != m:
        raise ValueError("Dimension of b must match rows of A.")
    if k <= 0 or k > m:
        raise ValueError("k must be a positive integer less than or equal to m.")

    Q = np.zeros((m, k + 1), dtype=A.dtype)
    H = np.zeros((k + 1, k), dtype=A.dtype)

    # Normalize the initial vector b and set it as the first column of Q
    Q[:, 0] = b / np.linalg.norm(b)

    for i in range(k):
        # Compute the new candidate vector v = A * q_i
        v = A @ Q[:, i]

        # Orthogonalize v against previous basis vectors in Q
        for j in range(i + 1):
            H[j, i] = np.vdot(Q[:, j], v)  # Inner product
            v -= H[j, i] * Q[:, j]

        # Compute the norm of the remaining vector v
        H[i + 1, i] = np.linalg.norm(v)

        # Check for breakdown (if norm is close to zero)
        if H[i + 1, i] < 1e-12:  # Using a small tolerance for numerical stability
            # Algorithm has converged or encountered a breakdown, return truncated matrices
            return Q[:, :i+1], H[:i+1, :i]

        # Normalize v and add it to Q as the next basis vector
        Q[:, i + 1] = v / H[i + 1, i]

    return Q[:, :k+1], H[:k+1, :k]

if __name__ == "__main__":
    # Example Usage:
    A = np.array([[1, 2, 0],
                  [2, 1, 2],
                  [0, 2, 1]], dtype=float)
    b = np.array([1, 0, 0], dtype=float)
    k_iterations = 3

    Q, H = arnoldi_iteration(A, b, k_iterations)

    print("Orthonormal Basis Q:")
    print(Q)
    print("\nUpper Hessenberg Matrix H:")
    print(H)

    # Verify Q is orthonormal (Q^T * Q should be identity)
    print("\nQ^T @ Q:")
    print(Q.T @ Q)

    # Verify A * Q_k = Q_{k+1} * H_k (approximately)
    # A @ Q[:, :k_iterations] should be approximately equal to Q[:, :k_iterations+1] @ H
    print("\nA @ Q[:, :k_iterations]:")
    print(A @ Q[:, :k_iterations])
    print("\nQ[:, :k_iterations+1] @ H:")
    print(Q[:, :k_iterations+1] @ H)

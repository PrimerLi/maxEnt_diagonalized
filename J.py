import default

def J(alpha, K_used, A, omega, sigma):
    import numpy as np
    import numpy.linalg

    domega = omega[1] - omega[0]

    matrix = K_used
    nuse = len(sigma)
    Nomega = len(omega)

    for i in range(nuse):
        matrix[i] = matrix[i]/sigma[i]**2
        
    result = -np.transpose(K_used).dot(matrix)

    for i in range(Nomega):
        result[i, i] = result[i, i] - alpha*domega/A[i]
    return result

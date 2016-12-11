import kernel

def chi(G_used, K_used, A, omega, sigma):
    import numpy as np
    import numpy.linalg

    nuse = len(G_used)
    Nomega = len(omega)

    vector_orginal = G_used - K_used.dot(A)
    vector = np.zeros(2*nuse)
    for i in range(nuse):
        vector[i] = vector_original[i]/sigma[i]**2
        vector[i+nuse] = vector_original[i+nuse]/sigma[i]**2

    result = np.transpose(vector_original).dot(vector)

    return result

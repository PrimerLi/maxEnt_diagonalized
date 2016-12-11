import default
import kernel

def f(alpha, G_used, K_used, A, omega, sigma):
    import numpy as np
    import numpy.linalg 

    domega = omega[1] - omega[0]
    Nomega = len(omega)
    result = np.zeros(Nomega)
    for nw in range(Nomega):
        result[nw] = -alpha*domega*(1 + np.log(A[nw]/default.D(omega[nw])))

    vector_right = G_used - K_used.dot(A)
    nuse = len(sigma)
    for i in range(nuse):
        vector_right[i] = vector_right[i]/sigma[i]**2
        vector_right[i+nuse] = vector_right[i+nuse]/sigma[i]**2
    result = result + np.transpose(K_used).dot(vector_right)

    return result

#!/usr/bin/env python

import default
import kernel
import f
import J
import chi
import entropy
import newton
import printFile
import nan

def read_spectral(fileName):
    import numpy as np
    omega = []
    A = []
    ifile = open(fileName, "r")
    for index, string in enumerate(ifile):
        a = string.split()
        omega.append(float(a[0]))
        A.append(float(a[1]))
    ifile.close()
    return omega, np.asarray(A)

def readFiles(Greal, Gimag):
    import os
    import sys
    import numpy as np
    import sys
    
    omega_n = []
    G_real = []
    G_imag = []
    
    try:
        ifile = open(Greal, "r")
    except:
        sys.exit(Greal + " does not exist. ")
    for index, string in enumerate(ifile):
        a = string.split()
        omega_n.append(float(a[0]))
        G_real.append(float(a[1]))
    ifile.close()

    try:
        ifile = open(Gimag, "r")
    except:
        sys.exit(Gimag + " does not exist. ")
    for index, string in enumerate(ifile):
        a = string.split()
        G_imag.append(float(a[1]))
    ifile.close()

    G_real = np.asarray(G_real)
    G_imag = np.asarray(G_imag)
    
    return omega_n, G_real, G_imag


def main():
    import os
    import sys
    import numpy as np
    import numpy.linalg

    if (len(sys.argv) == 1):
        print "a0 = sys.argv[1], b0 = sys.argv[2]. alpha = a0*exp(-i*b0). "
        return -1
    if (len(sys.argv) == 2):
        a0 = float(sys.argv[1])
        b0 = 0.05
    if (len(sys.argv) == 3):
        a0 = float(sys.argv[1])
        b0 = float(sys.argv[2])
    
    Greal = "G_cc_real.txt"
    Gimag = "G_cc_imag.txt"
    omega_n, G_real, G_imag = readFiles(Greal, Gimag)
    Niom = len(omega_n)
    
    Nomega = 101
    omega_lower = -5
    omega_upper = -omega_lower
    domega = (omega_upper - omega_lower)/float(Nomega-1)
    omega = np.zeros(Nomega)
    for i in range(Nomega):
        omega[i] = omega_lower + i*domega
    A_initial = np.zeros(Nomega)
    model = np.zeros(Nomega)
    for i in range(Nomega):
        model[i] = default.D(omega[i])
    printFile.printFile(omega, model, "model.txt")
    if (not os.path.exists("A_initial.txt")):
        for i in range(len(A_initial)):
            A_initial[i] = default.D(omega[i])
        printFile.printFile(omega, A_initial, "A_initial.txt")
    else:
        omega, A_initial = read_spectral("A_initial.txt")
        Nomega = len(omega)
    
    C_real = np.zeros((Niom, Niom))

    ifile = open("CM_cc_real.txt", "r")
    for (index, string) in enumerate(ifile):
        a = string.split()
        rowIndex = int(a[0])-1
        colIndex = int(a[1])-1
        C_real[rowIndex, colIndex] = float(a[2])
    ifile.close()

    eigenvalues, eigenvectors = numpy.linalg.inv(C_real)
    Omatrix = np.zeros(C_real.shape)
    for i in range(Niom):
        for j in range(Niom):
            Omatrix[i, j] = eigenvectors[i][j]

    K_matrix_real = kernel.K_matrix_real(omega_n, omega)
    K_matrix_imag = kernel.K_matrix_imag(omgea_n, omega)
    
    G_real_rotated = np.transpose(Omatrix).dot(G_real)
    G_imag_rotated = np.transpose(Omatrix).dot(G_imag)
    K_matrix_real_rotated = np.transpose(Omatrix).dot(K_matrix_real)
    K_matrix_imag_rotated = np.transpose(Omatrxi).dot(K_matrix_imag)
    
    sigma = []
    eps = 1.0e-12
    ofile = open("eps.txt", "w")
    ofile.write(str(eps) + "\n")
    ofile.close()
    for i in range(len(eigenvalues)):
        if (eigenvalues[i] > eps):
            sigma.append(np.sqrt(eigenvalues[i]))

    if (not os.path.exists("sigma.txt")):
        ofile = open("sigma.txt", "w")
        for i in range(len(sigma)):
            ofile.write(str(sigma[i]) + "\n")
        ofile.close()

    nsue = len(sigma)
    G_used = np.zeros(nuse*2)
    K_used = np.zeros((2*nuse, Nomega))
    for i in range(nuse):
        G_used[i] = G_real_rotated[i]
        G_used[i+nuse] = G_imag_rotated[i]
    for i in range(nuse):
        for j in range(Nomega):
            K_used[i, j] = K_matrxi_real_rotated[i,j]
            K_used[i+nuse, j] = K_matrix_imag_rotated[i,j]

    if (not os.path.exists("G_used.txt")):
        ofile = open("G_used.txt", "w")
        for i in range(len(G_used)):
            ofile.write(str(G_used[i]) + "\n")
        ofile.close()
    if (not os.path.exists("K_used.txt", "w")):
        printFile.printMatrix(K_used, "K_used.txt")

    if (True):
        alpha = []
        for i in range(30):
            alpha.append(a0*np.exp(-i*b0))
        
        ofile = open("alpha.txt", "a")
        for i in range(len(alpha)):
            A_updated = newton.newton(alpha[i], G_used, K_used, omega, A_initial, sigma)
            if (nan.array_isnan(A_updated)):
                omega, A_initial = read_spectral("A_initial.txt")
                continue
            output = "A_updated_alpha_" + str(alpha[i]) + ".txt"
            printFile.printFile(omega, A_updated, output)
            os.system("cp " +  output + " A_initial.txt")
            ofile.write(str(alpha[i]) + "\n")
            print "alpha = ", alpha[i]
        ofile.close()
    else:
        alpha = 0.01
        print "alpha = ", alpha
        A_updated = newton(alpha, G_used, K_used, omega, A_initial, sigma)
        if (nan.array_isnan(A_updated)):
            printFile.printFile(omega, A_updated, "A_updated_alpha_" + str(alpha) + ".txt")
            os.system("cp A_updated_alpha_" + str(alpha) + ".txt" + "A_initial.txt")
    return 0

main()

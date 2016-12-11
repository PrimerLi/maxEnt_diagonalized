import chi
import entropy
import numpy as np

def readFile(fileName):
    import os
    import sys
    import numpy as np

    quantity = []
    try:
        ifile = open(fileName, "r")
    except:
        sys.exit(fileName + " does not exist. \n")
    for index, string in enumerate(ifile):
        quantity.append(float(string))
    return np.asarray(quantity)

def readA(fileName):
    omega = []
    A = []
    ifile = open(fileName, "r")
    for i, string in enumerate(ifile):
        a = string.split()
        omega.append(float(a[0]))
        A.append(float(a[1]))
    ifile.close()
    return omega, A

def main():
    import os
    import sys
    import numpy.linalg

    sigma = readFile("sigma.txt")
    alpha = readFile("alpha.txt")
    G_used = readFile("G_used.txt")
    ifile = open("model.txt", "r")
    Nomega = len(enumerate(ifile))
    ifile.close()
    K_used = np.zeros((len(G_used), Nomega))
    ifile = open("K_used.txt", "r")
    for index, string in enumerate(ifile):
        a = string.split()
        row = int(a[0])
        col = int(a[1])
        K_used[row, col] = float(a[2])
    ifile.close()
    
    spectrals = []
    probability = []
    chi_values = []
    entropy_values = []
    for i in range(len(alpha)):
        print alpha[i]
        fileName = "A_updated_alpha_" + str(alpha[i]) + ".txt"
        omega, A = readA(fileName)
        A = np.asarray(A)
        spectrals.append(A)
        chi_values.append(chi.chi(G_used, K_used, A, omega, sigma))
        entropy_values.append(entropy.entropy(omega, A))
        probability.append(np.exp(alpha[i]*entropy_values[i] - 0.5*chi_values[i]))

    spectral_mean = np.zeros(len(spectrals[0]))
    for i in range(1, len(spectrals)):
        spectral_mean[:] = spectral_mean[:] + ((-alpha[i] + alpha[i-1])*probability[i]*spectrals[i])[:]
    s = 0.0
    for i in range(1, len(alpha)):
        s = s + (-alpha[i] + alpha[i-1])*probability[i]
    print "s = ", s
    ofile = open("Chi_log_alpha.txt", "w")
    for i in range(len(alpha)):
        ofile.write(str(np.log(alpha[i])) + "   " + str(chi_values[i]) + "\n")
    ofile.close()
    ofile = open("entropy_log_alpha.txt", "w")
    for i in range(len(alpha)):
        ofile.write(str(np.log(alpha[i])) + "    " + str(entropy_values[i]) + "\n")
    ofile.close()
    ofile = open("P_log_alpha.txt", "w")
    for i in range(len(alpha)):
        ofile.write(str(np.log(alpha[i])) + "   " + str(probability[i]) + "\n")
    ofile.close()
    ofile = open("bryan.txt", "w")
    for i in range(len(omega)):
        ofile.write(str(omega[i]) + "    " + str(spectral_mean[i]/s) + "\n")
    ofile.close()
    ofile = open("bryan_2pi.txt", "w")
    for i in range(len(omega)):
        ofile.write(str(omega[i]) + "    " + str(spectral_mean[i]/(2*np.pi*s)) + "\n")
    ofile.close()
    return 0

main()

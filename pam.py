from functions import *
import numpy as np
import matplotlib.pyplot as plt

# Ερώτημα Α
def calculate_SNRbdB(M, Pb):
    f = (Pb * M * np.log2(M)) / (2 * (M - 1)) # συνάρτηση f(x)
    SNRb = ((qfuncinv(f))**2 * (M**2 - 1)) / (6 * np.log2(M))
    return convert_to_dB(SNRb) # επέστρεψε το αποτέλεσμα σε dB

# Προσδιορίζει αν η συνάρτηση calculate_SNRbdB δουλεύει σωστά
def validate(M, Pb):
    # υπολογίζει το SNR με βάση τα M και Pb που περνάμε ως ορίσματα
    SNRbdB = calculate_SNRbdB(M, Pb)

    # χρησιμοποιεί το SNR για να βρει το νέο Pb
    y = np.sqrt((6 * convert_from_dB(SNRbdB) * np.log2(M)) / (M**2 - 1))
    pb =  2 * ((M - 1) / (M * np.log2(M))) * qfunc(y)

    print(f"SNR = {SNRbdB}", end = " ")

    # αν τα δύο Pb είναι ίσα τότε το SNR είναι έγκυρο
    if np.isclose(Pb, pb):
        print("is valid.")
    else:
        print("is not valid.")

# Επαλήθευση για Μ = 16 και Pb = 10^(-10)
validate(M = 16, Pb = 1e-10)

# Ερώτημα Β
def plot_SNRbdB(x):
    Pb = 10**(-x-2)
    M_values = [2**m for m in range(1, 11)]
    SNRbdB_values = [calculate_SNRbdB(M, Pb) for M in M_values]

    plt.plot(M_values, SNRbdB_values, 'o')
    plt.title(f'Pb = {Pb}')
    plt.xlabel('M')
    plt.ylabel('f(M)')
    
    # τα M εμφανίζονται σε δυνάμεις του 2
    plt.xscale('log', base = 2)
    plt.xticks(M_values, [str(m) for m in M_values])

    # κατακόρυφες γραμμές για κάθε M
    plt.grid(axis='x', linestyle='--', linewidth=0.8)

    plt.show()

# Ο αριθμός μητρώου μου είναι 2021088
plot_SNRbdB(8)

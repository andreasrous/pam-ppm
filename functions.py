from scipy import special as sp
import numpy as np
from scipy.special import erfc
import matplotlib.pyplot as plt

# η συνάρτηση Q, υλοποιημένη απ' την commlib
def qfunc(x):
    return 0.5 * erfc(x / np.sqrt(2))

# η αντίστροφη της Q
def qfuncinv(y):
    return np.sqrt(2) * sp.erfinv(1 - 2 * y)
    
# μετατρέπει το SNRb σε dB
def convert_to_dB(SNRb):
    return 10 * np.log10(SNRb)

# μετατρέπει το SNRb από dB στην αρχική του τιμή
def convert_from_dB(SNRbdB):
     return 10 ** (SNRbdB /10)

# μετατρέπει οποιοδήποτε όνομα (ή string) σε χαρακτήρες ASCII των 8-bit
def name_to_binary(name):
    # πίνακας με τα bits που αντιπροσωπεύουν το όνομα
    binary_representation = []

    for char in name:
        # μετέτρεψε κάθε χαρακτήρα στον 8-bit ASCII κώδικά του
        binary_str = format(ord(char), '08b')
        # μετέτρεψε το string που προκύπτει σε πίνακα ακεραίων από 0 και 1
        binary_array = [int(bit) for bit in binary_str]
        # προσάρτησε τον πίνακα στο αποτέλεσμα
        binary_representation.extend(binary_array)

    return np.array(binary_representation, dtype=int)

# PAM symbol constellation - από την commlib
def pam_constellation(M, beta = 1):
    m = np.arange(1, M + 1).astype(int)
    return (2 * m - M - 1) * beta

# από την commlib
def gray_code(m):
    if m == 1:
        g = ['0', '1']
    elif m > 1:
        gs = gray_code(m-1)
        gsr = gs[::-1]
        gs0 = ['0' + x for x in gs]
        gs1 = ['1' + x for x in gsr]
        g = gs0 + gs1
    return g

# από την commlib
def pam_gray_map(M, beta = 1):
    gc = gray_code(np.log2(M))
    Am = pam_constellation(M, beta = beta)
    pam_map = []

    for i, cw in enumerate(gc):
        pam_map.append([i, gc[i], Am[i]])
    
    return pam_map

# build a dictionary like map for faster encoding - από την commlib
def pam_gray_forward_map(M, beta = 1):
    pam_map = pam_gray_map(M, beta = beta)
    symbols = [x[2] for x in pam_map]
    bits = [x[1] for x in pam_map]
    forward_map = {}
    
    for i, symbol in enumerate(symbols):
        key = bits[i]
        forward_map[key] = symbol
        
    return forward_map

# από την commlib - τροποποιημένη
def array_to_str(a, m):
    astr = [str(x) for x in a]
    while len(astr) < m: # συμπληρώνει (αν χρειαστεί) με μηδενικά την τελευταία ομάδα bits 
        astr.insert(0, '0') # τα έξτρα μηδενικά μπαίνουν στην αρχή, άρα το π.χ. 11 γίνεται 011 για m = 3
    return ''.join(astr)

# bits to symbols - από την commlib
def bits_to_symbols(bits, fmap, return_bits = False):

    M = len(fmap)
    m = np.log2(M).astype(int)
    Nbits = bits.size
    
    symbols = []
    bitgroups = []
    i = 0
    j = 0
    
    while i < Nbits:
        key = array_to_str(bits[i : i + m], m)
        bitgroups.append(key)
        symbols.append(fmap[key])
        i += m
        j += 1
    if not return_bits:    
        return np.array(symbols)    
    else:
        return np.array(symbols), bitgroups

# από την commlib - τροποποιημένη για ppm
def plot_signal(t, x, TS, plot_type = None , close_all = False,
                xlabel = 't', ylabel = 'x(t)', figure_no = None,
                xlim = None, ylim = None, show_grid = False, title = None):
    
    if close_all:
        plt.close('all')
     
    if figure_no is None:
        plt.figure()
    else:
        plt.figure(figure_no)
    
    if plot_type == None:
        plt.plot(t, x, color = 'black', linewidth=0.8)
    else:
        plt.plot(t, x, plot_type, color = 'black', linewidth=0.8)

    plt.fill_between(t, 0, x, color='#3399ff')

    for i in range(int(np.ceil(t[-1] / TS))):
        plt.axvline(x=i * TS, color='gray', linestyle='--', linewidth=0.5)

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    
    if xlim is not None:
        plt.xlim(xlim)
    
    if ylim is not None:
        plt.ylim(ylim)
    
    if show_grid:
        plt.grid()

    if title is not None:
        plt.title(title)

# plot ppm waveform and show associated bits - από την commlib
# τροποποιημένη για ppm
def plot_ppm(t, x, M, bits, TS, tstart = 0, dy = 0, plot_type = None, title = None):
    
    plot_signal(t, x, TS, xlabel = 't', ylabel = 'p(t)', ylim = (-1, 3),
                plot_type = plot_type, title = title)
    
    Nbits = 0
    i = 0
    m = np.log2(M).astype(int)
    tcurrent = tstart
    
    while (i < Nbits) or (tcurrent <= np.max(t)):
       key = array_to_str(bits[i : i + m], m)
       plt.text(tcurrent + 0.5 * TS, 0.25 + dy, key, ha='center', va='center')

       i += m
       tcurrent += TS

# pulse position modulation : fast version - από την commlib
# τροποποιημένη για ppm
def ppm_waveform(ak, TS, M, fmap, tinitial = 0, tguard = 0.0):
    samples = 10 * M // 2
    Dt = TS / samples                                # sampling period    
    Nguard = np.round(tguard / Dt)                   # guard points                         
    Ntot = 2 * Nguard + samples * ak.size            # total number of points   
    
    x = np.zeros(Ntot.astype(int))
    t = np.arange(tinitial, tinitial + Ntot * Dt, Dt)

    i = np.floor((t - tinitial) / TS).astype(int)
    majority_value_replace(i, samples)
    j = np.where(np.logical_and(i >= 0, i < ak.size))

    x[j] = ak[i[j]]
    
    # αντιστοιχεί τα σύμβολα gray σε τιμές k όπου: 0 ≤ k ≤ M-1
    symbols_to_values = map_symbols_to_values(fmap)

    # μετατρέπει κάθε σύμβολο gray στον πίνακα x στην αντίστοιχη τιμή του k
    for idx in range(0, len(x)):
        x[idx] = symbols_to_values[x[idx]]

    # επεξεργάζεται τον πίνακα x και t, samples στοιχεία τη φόρα
    for idx in range(0, len(x), samples):
        chunk = x[idx:idx + samples]
        t_chunk = t[idx:idx + samples]
        k = chunk[0]
        v = t_chunk[0]
        pulse_start = round_number((k * TS / M) + v) # υπολογίζει τον χρόνο όπου ο παλμός πρέπει να γίνει 1
        pulse_end = round_number(((k + 1) * TS / M) + v) # υπολογίζει τον χρόνο όπου ο παλμός τελειώνει

        # βρίσκει τους δείκτες του πίνακα t που αντιστοιχούν στα pulse_start και pulse_end
        start_index = np.argmax(t >= pulse_start)
        end_index = np.argmax(t >= pulse_end)

        # διορθώνει το end_index στην περίπτωση των τελευταίων samples στοιχείων του πίνακα t
        if (end_index == 0):
            end_index = len(x)

        # θέτει τις τιμές μεταξύ start_index και end_index σε 1, και τις υπόλοιπες 0
        x[start_index:end_index] = 1
        x[idx:start_index] = 0
        x[end_index:idx + samples] = 0

    return t, x

# αντιστοιχεί τα σύμβολα gray σε τιμές k όπου: 0 ≤ k ≤ M-1
# η αντιστοίχιση γίνεται με βάση τον κώδικα Gray: διαδοχικοί παλμοί απέχουν κατά ένα bit
def map_symbols_to_values(fmap):
    unique_symbols = list(set(fmap.values())) # παίρνει τα σύμβολα gray
    unique_symbols.sort() # τα ταξινομεί σε αύξουσα σειρά π.χ. (-3, -1, 1, 3) για M = 4

    # αντιστοιχεί τα σύμβολα σε τιμές k ώστε διαδοχικοί παλμοί να απέχουν κατά ένα bit
    # π.χ. {-3: 0, -1: 1, 1: 2, 3: 3} ή {00: 0, 01: 1, 11: 2, 10: 3}
    # ΠΡΟΣΟΧΗ ότι το σύμβολα 01 και 11 απέχουν κατά ένα bit,
    # επομένως αν το 01 έχει k = 1, το 11 πρέπει να έχει k = 2
    symbol_to_value_map = {symbol: i for i, symbol in enumerate(unique_symbols)}
    return symbol_to_value_map

# διορθώνει την τιμή i της ppm_waveform, όπου μερικές τιμές εμφανίζονταν πάνω απο samples φορές
# π.χ. αν samples = 5 μπορεί να είχαμε [8. 8. 8. 8. 8. 8. 9. 9. 9. 9.]
# όπου το 8 εμφανίζεται έξι φορές αλλά το 9 τέσσερις (πρόβλημα της np.floor)
def majority_value_replace(array, chunk_size):
    # επεξεργάζεται τον πίνακα array, chunk_size στοιχεία την φορά
    for chunk in np.array_split(array, len(array) // chunk_size):
        # βρίσκει ποιο στοιχείο επαναλαμβάνεται τις περισσότερες φορές στον υποπίνακα
        # και αντικαθιστά όλα τα στοιχεία του υποπίνακα με αυτό το στοιχείο
        majority_value = np.bincount(chunk).argmax()
        chunk[:] = majority_value

# στρογγυλοποιεί τους αριθμούς για πιο ακριβείς υπολογισμούς
def round_number(num):
    return float(f'{num: .3e}') # π.χ. το 1.50000002e-10 θα γίνει 1.500e-10
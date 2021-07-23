import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import h5py as h5

# TeX labels
plt.rc('text', usetex=True)
plt.rc('font', family='TeX Gyre Adventor', size=14)
plt.rc('pdf', fonttype=42)

def plot_excitations(filename, title, saveimage, cmax=0.005):
    plt.figure()
    with h5.File(filename, 'r') as f:
        number_spins = 6 * f.attrs['L'] ** 2
        kpath = f["kpath"][:]
        Sqomega = f["Sqω"][:]
        
    plt.imshow(np.abs(Sqomega)[1:100, :] / number_spins ** 2,
               origin='lower',
               aspect='auto')
        
    plt.clim(0, cmax)
    plt.colorbar()
    
    plt.xticks([0, 100, 200, 300, 400], ['$\\Gamma$', '$(2\\pi, 0)$', '$(4\\pi, 0)$', '$(2\\pi, 2\\pi)$', '$\\Gamma$'])
    plt.yticks([0, 50, 100], ['$0$', '$\\pi$', '$2\\pi$'])
    
    plt.ylabel('$\\omega / J$')
    plt.title(title)
    
    plt.tight_layout()
    plt.savefig(saveimage)

# nice UUD magnon branches
plot_excitations(filename='../data/results/uud_20x20.h5',
                 title='$|S(\\vec q, \\omega)|$ for the UUD phase',
                 saveimage='../plots/sqomega_uud.png')
plt.show()

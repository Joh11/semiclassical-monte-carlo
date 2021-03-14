import numpy as np
import matplotlib.pyplot as plt
import h5py as h5

def load_and_mean():
    with h5.File("kagome-fig4-t40.h5", "r") as f:
        # find the number of samples
        nsamples = 0
        while True:
            try:
                f[f"{nsamples + 1}"]
                nsamples += 1
            except KeyError:
                break

        print(f'Found {nsamples} samples ...')
        Sqt = np.abs(np.transpose(f[f'1/Sqt'], (2, 1, 0)))
        Sqomega = np.abs(np.transpose(f[f'1/Sqω'], (2, 1, 0)))
        vs = np.transpose(f[f'1/vs'], (4, 3, 2, 1, 0))

        for i in range(2, nsamples + 1):
            print(f'Doing {i} ...')
            Sqt += np.abs(np.transpose(f[f'{i}/Sqt'], (2, 1, 0)))
            Sqomega += np.abs(np.transpose(f[f'{i}/Sqω'], (2, 1, 0)))
            vs += np.transpose(f[f'{i}/vs'], (4, 3, 2, 1, 0))

        Sqt /= nsamples
        Sqomega /= nsamples
        vs /= nsamples
    return vs, Sqt, Sqomega

vs, Sqt, Sqomega = load_and_mean()

def single_spin(vs, i, j, s):
    "Returns the time evolution of a single spin"
    return vs[:, s, i, j, :]

# plot Sqt
# ========
plt.figure()
normSqt = Sqt / Sqt[:, :, 0].reshape(Sqt.shape[0], -1, 1)
for h in range(0, 10, 3):
    plt.scatter(0.1 * np.arange(Sqt.shape[-1]), normSqt[h, h],
                label=f'|q| = {np.sqrt(2) * np.fft.fftfreq(Sqt.shape[0])[h]:.2f}')

plt.xlabel('t []')
plt.ylabel(r'$S(\vec q, t) / S(\vec q, 0)$')
plt.legend()
plt.show()
# plt.savefig("Sqt_L30.png")

# plot spin
# =========
# i, j, s = 1, 0, 0
# S = single_spin(vs, i, j, s)

# plt.plot(np.arange(vs.shape[-1]), S[0], label='x')
# plt.plot(np.arange(vs.shape[-1]), S[1], label='y')
# plt.plot(np.arange(vs.shape[-1]), S[2], label='z')

# plt.legend()

# plt.show()


# plot S(Q, 0)
# ========
# plt.figure()

# t = 40
# plt.imshow(np.log(Sqt[:, :, t]),
#            origin='lower',
#            extent=(-0.5, 0.5, -0.5, 0.5),
#            cmap='inferno')
# plt.colorbar()

# plt.show()

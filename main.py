import sys
import numpy as np
import matplotlib.pyplot as plt

from modlamp.sequences import Random
from modlamp.descriptors import PeptideDescriptor

AM = 'ACDEFGHIKLMNPQRSTVWY'


class EstimateDoS:
    num_exp = 0

    def __init__(self, numiter, plen, E_max=0, E_min=-3, num_discretize=50):
        self.numiter = numiter
        self.plen = plen
        self.lib = Random(1, lenmin=plen, lenmax=plen)
        self.lib.generate_sequences(proba='rand')
        self.E_max = E_max
        self.E_min = E_min
        self.num_discretize = num_discretize
        self.S = np.array([1.0] * num_discretize)

    def incr_exp(self):
        self.num_exp += 1

    def cal_prop(self, seqs):
        self.incr_exp()
        desc = PeptideDescriptor(str(seqs), 'TM_tend')
        desc.calculate_moment()
        mom = desc.descriptor[0][0]
        return -mom

    def histogram_idx(self, E):
        if E < self.E_min:
            return 0
        elif E >= self.E_max:
            return self.num_discretize - 1
        else:
            return int((E - self.E_min) / (self.E_max - self.E_min) * self.num_discretize)

    def plot(self):
        plt.clf()
        plt.close()
        E = []
        g = []
        for i in range(self.num_discretize):
            E.append((self.E_max - self.E_min) * (i + 0.5) / self.num_discretize + self.E_min)
            g.append(np.exp(self.S[i] - np.max(self.S)))
        plt.plot(E, g)
        plt.savefig('gE.png')

    def wang_landau(self):
        f = 1.0
        x_current = self.lib.sequences[0]
        E_current = self.cal_prop(x_current)
        index_current = self.histogram_idx(E_current)

        for f_itt in range(self.numiter):
            print("f_itt = ", f_itt)
            H = np.zeros(self.num_discretize)

            while True:
                x_new = x_current
                select_dim = np.random.randint(0, self.plen)
                am_candidate = AM.replace(x_current[select_dim], '')
                x_new = x_current[0:select_dim] + np.random.choice(list(am_candidate)) + x_current[select_dim + 1:]
                E_new = self.cal_prop(x_new)

                index_new = self.histogram_idx(E_current)
                if self.S[index_current] > self.S[index_new] or np.exp(self.S[index_current] - self.S[index_new]) > np.random.rand():
                    x_current = x_new
                    E_current = E_new
                    index_current = index_new

                self.S[index_current] += f
                H[index_current] += 1

                if np.min(H[np.nonzero(H)]) / np.max(H[np.nonzero(H)]) > 0.95 and np.sum(H) > 10000:
                    break

            print("f =", f)
            print("total sampling number =", int(np.sum(H)))

            f = f / 2

if __name__ == '__main__':
    numiter = int(sys.argv[1])
    plen = 5
    estimate_dos = EstimateDoS(numiter, plen)
    estimate_dos.wang_landau()
    estimate_dos.plot()

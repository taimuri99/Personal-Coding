import numpy as np
import math


class Spin(object):
    def __init__(self, pos, M):
        self.pos = pos
        self.M = M


class SpinGlass(object):
    """
    Defines a spin glass over a 2D lattice with periodic boundary conditions.
    Initial states is paramagnetic.
    self.T = "Thermal energy of the system", equal to Temperature*kB
    self.side = "Side length"
    self.config = "Spin configuration"
    self.M = "Total magnetization"
    self.E = "Total energy"
    """

    def __init__(self, T, side):
        self._T = T
        self._side = side
        self._config = {}
        # paramagnetic init
        np.random.seed()
        for x in range(side):
            for y in range(side):
                M = int(np.random.randint(0, 2) * 2 - 1)
                self._config[(x, y)] = Spin((x, y), M)
        self._M = sum([s.M for s in self._config.values()])
        self._E = self.hamiltonian(self._config.values())

    @property
    def T(self):
        """Temperature"""
        return self._T

    @T.setter
    def T(self, Temp):
        self._T = Temp

    @property
    def side(self):
        """Side of the lattice in units of length"""
        return self._side

    @property
    def config(self):
        """The spin configuration in a side*side matrix"""
        config_mat = np.zeros(self._side * self._side).reshape((self._side, self._side))
        for (x, y), s in self._config.items():
            config_mat[x, y] = s.M
        return config_mat

    @property
    def M(self):
        return self._M

    @property
    def E(self):
        return self._E

    def hamiltonian(self, spins):
        energy = 0.
        for s in spins:
            # Ising model has no kinetic energy term. There is no coupling with external fields.
            for n in self.first_neigh(s):
                # Two body Interaction with first neighbors
                energy += - s.M * n.M / 2.
        return energy

    def first_neigh(self, s):
        neighbors = []
        for d in [0, 1]:
            for n in [-1, +1]:
                n_pos = list(s.pos)
                n_pos[d] = (n_pos[d] + n + self._side) % self._side
                neighbors.append(self._config[tuple(n_pos)])
        return neighbors

    def metropolis_flipper(self):
        for s in self._config.values():
            delta_E = -2.*2. * self.hamiltonian([s])  # the "minus" signs act as if the spin has been flipped
            acceptance = math.exp(-delta_E / self._T)
            if acceptance >= 1:
                s.M *= -1
                self._M += 2 * s.M
                self._E += delta_E
            elif acceptance > np.random.random():
                s.M *= -1
                self._M += 2 * s.M
                self._E += delta_E


    def __iter__(self):
        return iter(self._config.values())

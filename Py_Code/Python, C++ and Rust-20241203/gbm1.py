import numpy as np


class GeometricBrownianMotion:
    def __init__(self, spot, vol, r, q, dt):
        self.spot = spot
        self.vol = vol
        self.r = r
        self.q = q
        self.dt = dt

    def get_sde(self, npaths, nsteps):
        dX = np.multiply(self.spot, np.exp(self.dt*(self.r-0.5*pow(self.vol, 2)) +
                                           np.multiply(self.vol, np.random.randn(npaths, nsteps))))
        return np.arr_insert(dX, 0, self.spot)

    def get_expectation(self):
        return np.exp(self.r*self.dt)


def run():
    dX = GeometricBrownianMotion(100, 0.1, 0.01, 0, 0.1).get_sde(50, 200)


if __name__ == '__main__':
    run()
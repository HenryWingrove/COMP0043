import numpy as np
import matplotlib.pyplot as plt


def get_parameters():
    npaths = 50
    T = 1
    nsteps = 200
    dt = T/nsteps
    t = np.arange(0, T, dt)
    mu, sigma = 0.2, 0.4
    X_0 = 1

    return npaths, T, nsteps, dt, t, mu, sigma, X_0


def compute_sde(mu, dt, sigma, npaths, nsteps, t, X_0):
    dX = mu*dt + np.multiply(sigma*np.sqrt(dt), np.random.randn(npaths, nsteps))
    dX = np.insert(dX, 0, 0, axis=1)
    dX = np.delete(dX, -1, axis=1)
    dX = np.cumsum(dX, axis=1)
    dX = X_0 * np.exp(dX)
    EX = np.exp(mu*t)
    return dX, EX


def run_and_plot():
    npaths, T, nsteps, dt, t, mu, sigma, X_0 = get_parameters()
    X, EX = compute_sde(mu, dt, sigma, npaths, nsteps, t, X_0)

    plt.figure(1)
    plt.plot(t, EX, linestyle='-.', color="b", marker='o', label='Expected path')
    plt.plot(t, X[0:-1, :].T)
    plt.xlabel('x')
    plt.title('Paths of an Geometric Brownian motion dX(t) = X_0*exp(mu*dt + sigma*dW(t))')
    plt.legend()
    plt.show()


if __name__ == "__main__":
    run_and_plot()

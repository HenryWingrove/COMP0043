from math import (
    pi
)
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt


def get_parameters():
    mu = 0.2
    sigma = 0.1
    a = -0.4
    b = 0.8
    ngrid = 120
    nsample = np.power(10, 6)

    return mu, sigma, a, b, ngrid, nsample


def get_grid(a, b, ngrid):
    x = np.linspace(a, b, ngrid+1)
    deltax = x[1]-x[0]

    return x, deltax


def compute_pdf_and_cdf(x, sigma, mu):
    f1 = np.divide(1, np.sqrt(2*pi)*sigma) * np.exp(-np.multiply((x-mu)/sigma, 2)/2)
    f = stats.norm.pdf(x, mu, sigma)
    F = stats.norm.cdf(x, mu, sigma)

    return f1, f, F


def sample_norm_dist(mu, sigma, nsample):
    return mu+sigma*np.random.randn(nsample)


def run_and_plot():
    mu, sigma, a, b, ngrid, nsample = get_parameters()
    x, deltax = get_grid(a, b, ngrid)
    f1, f, F = compute_pdf_and_cdf(x, sigma, mu)

    plt.figure(1)
    plt.plot(x, f1, linestyle='-.', color="b", marker='o', label='pdf')
    plt.xlabel('x')
    plt.title('Normal distribution with mu = 0.2 and sigma = 0.1')
    plt.legend()
    plt.show()

    plt.figure(2)
    plt.plot(x, f, linestyle='-.', color="b", marker='o', label='PDF')
    plt.plot(x, F, linestyle=':', color="r", marker='D', label='CDF')
    plt.xlabel('x')
    plt.title('Normal distribution with mu = 0.2 and sigma = 0.1')
    plt.legend()
    plt.show()


if __name__ == "__main__":
    run_and_plot()

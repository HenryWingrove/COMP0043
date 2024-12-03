import numpy as np
from scipy import stats
import matplotlib.pyplot as plt


def get_parameters():
    mu = 0.2
    sigma = 0.1
    a = 0.0
    b = 2
    ngrid = 120
    nsample = np.power(10, 6)

    return mu, sigma, a, b, ngrid, nsample


def get_grid(a, b, ngrid):
    x = np.linspace(a, b, ngrid+1)

    return x


def compute_pdf_and_cdf(x, sigma, mu):
    f = stats.lognorm.pdf(x, mu, sigma)
    F = stats.lognorm.cdf(x, mu, sigma)

    return f, F


def sample_lognorm_dist(mu, sigma, nsample):
    return np.exp(mu+sigma*np.random.randn(nsample))


def run_and_plot():
    mu, sigma, a, b, ngrid, nsample = get_parameters()
    x = get_grid(a, b, ngrid)
    f, F = compute_pdf_and_cdf(x, sigma, mu)

    plt.figure(1)
    plt.plot(x, f, linestyle='-.', color="b", marker='o', label='PDF')
    plt.plot(x, F, linestyle=':', color="r", marker='D', label='CDF')
    plt.xlabel('x')
    plt.title('Lognormal distribution with mu = 0.2 and sigma = 0.1')
    plt.legend()
    plt.show()


if __name__ == "__main__":
    run_and_plot()

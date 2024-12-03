import numpy as np
from scipy import stats
import matplotlib.pyplot as plt


def get_parameters():
    n = 5
    non_central_par = 2
    a = 0.0
    b = 20.0
    ngrid = 100
    nsample = np.power(10, 6)

    return n, non_central_par, a, b, ngrid, nsample


def get_grid(a, b, ngrid):
    x = np.linspace(a, b, ngrid+1)

    return x


def compute_pdf_and_cdf(x, n, non_central_par):
    f = stats.chi2.pdf(x, n, non_central_par)
    F = stats.chi2.cdf(x, n, non_central_par)

    return f, F


def run_and_plot():
    n, non_central_par, a, b, ngrid, nsample = get_parameters()
    x = get_grid(a, b, ngrid)
    f, F = compute_pdf_and_cdf(x, n, non_central_par)

    plt.figure(1)
    plt.plot(x, f, linestyle='-.', color="b", marker='', label='PDF')
    plt.plot(x, F, linestyle=':', color="r", marker='', label='CDF')
    plt.xlabel('x')
    plt.title('Chi-squared distribution with mu = 0.2 and sigma = 0.1')
    plt.legend()
    plt.show()


if __name__ == "__main__":
    run_and_plot()

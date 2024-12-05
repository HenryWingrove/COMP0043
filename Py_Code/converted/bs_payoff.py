import numpy as np
from scipy.stats import norm
from scipy.fft import fft, ifft, fftshift, ifftshift
import matplotlib.pyplot as plt
import time

# Black-Scholes Formula
def black_scholes(S0, K, r, T, sigma, q):
    mu_abm = r - q - 0.5 * sigma**2  # Drift coefficient
    d2 = (np.log(S0 / K) + mu_abm * T) / (sigma * np.sqrt(T))
    d1 = d2 + sigma * np.sqrt(T)
    call = S0 * np.exp(-q * T) * norm.cdf(d1) - K * np.exp(-r * T) * norm.cdf(d2)
    put = K * np.exp(-r * T) * norm.cdf(-d2) - S0 * np.exp(-q * T) * norm.cdf(-d1)
    return call, put

# Direct Integration in Log-Price Space
def direct_integration(S0, K, r, T, sigma, q, xwidth, ngridx):
    mu_abm = r - q - 0.5 * sigma**2
    Nx = ngridx // 2
    dx = xwidth / ngridx
    x = dx * np.arange(-Nx, Nx)
    pdf = 1 / (np.sqrt(2 * np.pi * T) * sigma) * np.exp(-(x - mu_abm * T)**2 / (2 * sigma**2 * T))
    call_payoff = np.maximum(S0 * np.exp(x) - K, 0) * pdf
    put_payoff = np.maximum(K - S0 * np.exp(x), 0) * pdf
    call_price = np.exp(-r * T) * np.sum(call_payoff) * dx
    put_price = np.exp(-r * T) * np.sum(put_payoff) * dx
    return call_price, put_price

# Fourier Transform-Based Pricing
def fourier_pricing(S0, K, r, T, sigma, q, xwidth, ngrid, alpha):
    N = ngrid // 2
    dx = xwidth / ngrid
    x = dx * np.arange(-N, N)
    b = xwidth / 2
    dxi = np.pi / b
    xi = dxi * np.arange(-N, N)

    mu_abm = r - q - 0.5 * sigma**2
    psi = lambda xi: 1j * mu_abm * xi - 0.5 * (sigma * xi)**2

    G = np.exp(alpha * x) * np.maximum(S0 * np.exp(x) - K, 0)
    G_fft = fftshift(fft(ifftshift(G)))

    Psi = np.exp(psi(xi + 1j * alpha) * T)
    result = np.real(G_fft * np.conj(Psi))
    call_price = np.exp(-r * T) / np.pi * np.trapz(result[N:], dx=dxi)
    return call_price

# Monte Carlo Simulation
def monte_carlo(S0, K, r, T, sigma, q, nblocks, npaths):
    mu_abm = r - q - 0.5 * sigma**2
    call_prices = []
    put_prices = []

    for _ in range(nblocks):
        X = mu_abm * T + sigma * np.random.randn(npaths) * np.sqrt(T)
        S = S0 * np.exp(X)
        call_payoff = np.maximum(S - K, 0)
        put_payoff = np.maximum(K - S, 0)
        call_prices.append(np.mean(call_payoff) * np.exp(-r * T))
        put_prices.append(np.mean(put_payoff) * np.exp(-r * T))

    call_price = np.mean(call_prices)
    put_price = np.mean(put_prices)
    call_stdev = np.std(call_prices) / np.sqrt(nblocks)
    put_stdev = np.std(put_prices) / np.sqrt(nblocks)
    return call_price, put_price, call_stdev, put_stdev

# Main Function
def main():
    # Contract parameters
    T = 1.0
    K = 1.1
    S0 = 1.0
    r = 0.05
    q = 0.02
    sigma = 0.4

    # Direct integration parameters
    xwidth = 8
    ngridx = 2**16

    # Fourier parameters
    ngrid = 2**6
    alphac = -6
    alphap = 6

    # Monte Carlo parameters
    nblocks = 10000
    npaths = 2000

    # Analytical Solution
    start = time.time()
    bs_call, bs_put = black_scholes(S0, K, r, T, sigma, q)
    print(f"Black-Scholes Call: {bs_call:.10f}, Put: {bs_put:.10f}, Time: {time.time() - start:.5f}s")

    # Direct Integration
    start = time.time()
    di_call, di_put = direct_integration(S0, K, r, T, sigma, q, xwidth, ngridx)
    print(f"Direct Integration Call: {di_call:.10f}, Put: {di_put:.10f}, Time: {time.time() - start:.5f}s")

    # Fourier Transform Method
    start = time.time()
    ft_call = fourier_pricing(S0, K, r, T, sigma, q, xwidth, ngrid, alphac)
    print(f"Fourier Transform Call: {ft_call:.10f}, Time: {time.time() - start:.5f}s")

    # Monte Carlo Simulation
    start = time.time()
    mc_call, mc_put, mc_call_stdev, mc_put_stdev = monte_carlo(S0, K, r, T, sigma, q, nblocks, npaths)
    print(f"Monte Carlo Call: {mc_call:.10f}, Put: {mc_put:.10f}, Call Stdev: {mc_call_stdev:.10f}, Put Stdev: {mc_put_stdev:.10f}, Time: {time.time() - start:.5f}s")

if __name__ == "__main__":
    main()
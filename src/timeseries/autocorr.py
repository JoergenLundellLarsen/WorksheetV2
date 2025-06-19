import numpy as np
from statsmodels.tsa.stattools import acf
import matplotlib.pyplot as plt

def plot_autocorrelation(data, lags, title):
    acorr = acf(data, nlags=lags, fft=True)
    plt.figure(figsize=(10, 5))
    plt.stem(range(len(acorr)), acorr)
    plt.title(title)
    plt.xlabel("Lag")
    plt.ylabel("Autocorrelation")
    plt.grid(True)
    plt.show()

def compute_autocorrelation(data, lags=50, plot=True, title="Autocorrelation"):
    if plot:
        plot_autocorrelation(data, lags, title)
    return acf(data, nlags=lags, fft=True)

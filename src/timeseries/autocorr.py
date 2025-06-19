import numpy as np
from statsmodels.tsa.stattools import acf
import matplotlib.pyplot as plt

def plot_autocorrelation(data, lags, title, filename=None):
    """
    Plot the autocorrelation function of a time series.
    
    Parameters:
        data (array-like): Time series values.
        lags (int): Number of lags.
        title (str): Plot title.
        filename (str): If set, save the plot to file.
    """
    acorr = acf(data, nlags=lags, fft=True)
    plt.figure(figsize=(10, 5))
    plt.stem(range(len(acorr)), acorr)
    plt.title(title)
    plt.xlabel("Lag")
    plt.ylabel("Autocorrelation")
    plt.grid(True)
    if filename:
        plt.savefig(filename)
    plt.close()

def compute_autocorrelation(data, lags=50, plot=True, title="Autocorrelation", filename=None):
    """
    Compute autocorrelation, optionally plot it.

    Returns:
        np.ndarray: Autocorrelation values.
    """
    if plot:
        plot_autocorrelation(data, lags, title, filename)
    return acf(data, nlags=lags, fft=True)

def estimate_autocorr_time(data, cutoff=0.05):
    """
    Estimate the autocorrelation time.

    Parameters:
        data (array-like): Input time series.
        cutoff (float): Threshold for truncating the sum.

    Returns:
        int: Estimated autocorrelation time.
    """
    acorr = acf(data, fft=True, nlags=len(data)//2)
    act = 1
    for i in range(1, len(acorr)):
        if acorr[i] < cutoff:
            break
        act += 2 * acorr[i]
    return int(np.ceil(act))

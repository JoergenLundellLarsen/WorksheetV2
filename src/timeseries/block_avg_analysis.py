import numpy as np
import matplotlib.pyplot as plt

def block_average(data, num_blocks):
    """
    Perform block averaging on a time series.

    Parameters:
        data (array-like): Input data.
        num_blocks (int): Number of blocks to divide the data into.

    Returns:
        tuple: (mean, standard error, array of block means)
    """
    data = np.array(data)
    block_size = len(data) // num_blocks
    block_means = np.zeros(num_blocks)

    for i in range(num_blocks):
        start = i * block_size
        end = start + block_size
        block_means[i] = np.mean(data[start:end])

    mean = np.mean(block_means)
    error = np.std(block_means, ddof=1) / np.sqrt(num_blocks)
    return mean, error, block_means

def plot_block_averaging(data, max_blocks, title, filename=None):
    """
    Plot block averaged mean and error vs number of blocks.

    Parameters:
        data (array-like): Input data.
        max_blocks (int): Max number of blocks to test.
        title (str): Plot title.
        filename (str): Save path for the plot.
    """
    block_numbers = np.arange(1, max_blocks + 1)
    means = []
    errors = []

    for n in block_numbers:
        mean, error, _ = block_average(data, n)
        means.append(mean)
        errors.append(error)

    plt.figure(figsize=(10, 5))
    plt.errorbar(block_numbers, means, yerr=errors, fmt='o-', capsize=5)
    plt.title(title)
    plt.xlabel('Number of Blocks')
    plt.ylabel('Mean ± Std. Error')
    plt.grid(True)

    if filename:
        plt.savefig(filename)
    plt.close()

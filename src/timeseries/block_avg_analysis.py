import numpy as np
import matplotlib.pyplot as plt

def block_average(data, num_blocks):
    data = np.array(data)
    block_size = len(data) // num_blocks
    block_means = np.zeros(num_blocks)

    for i in range(num_blocks):
        start = i * block_size
        end = start + block_size
        block_means[i] = np.mean(data[start:end])

    overall_mean = np.mean(block_means)
    standard_error = np.std(block_means, ddof=1) / np.sqrt(num_blocks)

    return overall_mean, standard_error, block_means

def plot_block_averaging(data, max_blocks, title, filename=None):
    block_numbers = np.arange(1, max_blocks + 1)
    means = []
    errors = []

    for num_blocks in block_numbers:
        mean, error, _ = block_average(data, num_blocks)
        means.append(mean)
        errors.append(error)

    plt.figure(figsize=(10, 5))
    plt.errorbar(block_numbers, means, yerr=errors, fmt='o-', capsize=5)
    plt.title(title)
    plt.xlabel('Number of Blocks')
    plt.ylabel('Mean ± Standard Error')
    plt.grid(True)

    if filename:
        plt.savefig(filename)
    plt.close()

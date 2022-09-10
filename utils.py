# %%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import os

import config

"""
    I don't know what a utils folder does, but I've taken it to mean random
    but useful functions that reappear loads of places but break the flow
    of understanding.
"""


def plot_log_axes(x, y, filename=None, N=None, xlabel='Energy [eV]',
                  ylabel=r'Flux [cm$^{-2}$particle$^{-1}$]', other=False,
                  legend=None, title=''):
    fig, ax = plt.subplots()
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if type(y) is list:
        if type(x) is list:
            for i in range(len(y)):
                ax.plot(x[i], y[i])
        else:
            for i in range(len(y)):
                ax.plot(x, y[i])
    else:
        ax.plot(x, y)
    if legend:
        ax.legend(legend)
    if len(title) > 0:
        ax.set_title(title)

    # File saving
    if filename:
        if not N:
            raise Exception(
                "You must specify how many particles you have simulated")

        run_env = config.RUN_ENV
        if other:
            filepath = os.path.join(
                run_env, "graphs", f"other/{filename}.png")
            if not os.path.exists(os.path.dirname(filepath)):
                os.makedirs(os.path.dirname(filepath))
            plt.savefig(filepath, dpi=500)
            return
        if type(x) is list:
            subdir = 'comparison/'
        else:
            subdir = 'standard/'
        # Saving file to output
        filepath = os.path.join(run_env, 'graphs',
                                f'e{N}/{subdir}/{filename}.png')
        if not os.path.exists(os.path.dirname(filepath)):
            os.makedirs(os.path.dirname(filepath))
        plt.savefig(filepath, dpi=500)


# Ancient history, can ignore this function until you need it, but you'll need to edit it
def plot_histogram(x, y, filename, N, xlabel='Energy [eV]', ylabel=r'Flux \[cm$^{-2}$s$^2$\]', norm=False):
    fig, ax = plt.subplots()
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    left_edges = x[:-1]
    widths = (x[1:] - x[:-1]) * 0.85
    ax.bar(left_edges, y, width=widths, align='edge')

    # Saving file to output
    norm_path = 'normalised/' if norm else 'standard'
    filepath = f'graphs/e{N}/{norm_path}/{filename}_hist.png'
    if not os.path.exists(os.path.dirname(filepath)):
        os.makedirs(os.path.dirname(filepath))
    plt.savefig(filepath, dpi=1000)


def load_tally(filepath):
    df = pd.read_csv(filepath)
    df = df[["energy low [eV]", "energy high [eV]", "mean"]]
    df["F/dE"] = df['mean'] / (df["energy high [eV]"] - df["energy low [eV]"])
    df["F/dU"] = df['mean'] / \
        (np.log(df["energy high [eV]"] / df["energy low [eV]"]))
    df["mid_bins"] = (df['energy low [eV]'] + df['energy high [eV]']) / 2
    df["dU"] = np.log(df["energy high [eV]"] / df["energy low [eV]"])
    df["dE"] = df["energy high [eV]"] - df["energy low [eV]"]
    return df


if __name__ == '__main__':
    ...


# %%

# %%

"""
Measure the time complexity.
This program will measure the time complexity of running random sequences
of lengths {5, 10, 25, 50, 100, 200, 500} for n steps to get an average.

author: Reis Gadsden
version: 5/4/2021
git: https://github.com/reismgadsden/predict_sec_struct_with_cfg

class: CS-5531-101 @ Appalachian State University
instructor: Dr. Mohammad Mohebbi
"""

import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime # used to generate dates for file saving
import timeit # used to time the methods (time did not work on my machine)
import predict_with_cfg as pwc
import random


def main():
    sizes = [5, 10, 25, 50, 100, 250, 500]
    timing = []
    n = 10
    bases = "ACGU"
    for size in sizes:
        run_time = 0
        for i in range(0, n):
            sequence = ""
            for j in range(0, size):
                sequence += bases[random.randint(0, 3)]

            start = timeit.default_timer()
            output = pwc.predict_stochastic_secondary_structure(sequence)
            end = timeit.default_timer()
            run_time += end - start
            print(i)

        timing.append(run_time / n)

    plot_bars(timing, sizes, n)
    plt.clf()
    plt.cla()

    plot_bars_ylog(timing, sizes, n)
    plt.clf()
    plt.cla()

    plot_line(timing, sizes, n)
    plt.clf()

    plot_line_ylog(timing, sizes, n)
    plt.clf()


def plot_bars(timing, sizes, iterations):
    x = np.arange(len(sizes))

    fig, ax = plt.subplots()

    bar_graph = ax.bar(x, timing)
    ax.set_ylabel("Average run time (in seconds) over " + str(iterations) + " total iterations")
    ax.set_xticks(x, sizes)
    ax.set_title("Stochastic Context Free Grammar Base Pairing Run-Time")
    plt.xlabel("Length of Sequence")

    # apply tight layout just to make sure everything fits nicely
    fig.tight_layout()

    # save an image of the graph to the output folder
    # using the date to make sure every graph has unique name, name goes down to the second
    current_date = str(datetime.today()).strip().replace(" ", "_").replace("-", "_").replace(":", "_")
    plt.savefig("output/bar_" + current_date[0: current_date.index(".")] + ".png")


def plot_bars_ylog(timing, sizes, iterations):
    x = np.arange(len(sizes))

    fig, ax = plt.subplots()

    bar_graph = ax.bar(x, timing)
    ax.set_ylabel("Average run time (in seconds, log10) over " + str(iterations) + " total iterations")
    ax.set_xticks(x, sizes)
    ax.set_title("Stochastic Context Free Grammar Base Pairing Run-Time")
    plt.xlabel("Length of Sequence")
    plt.yscale("log")

    # apply tight layout just to make sure everything fits nicely
    fig.tight_layout()

    # save an image of the graph to the output folder
    # using the date to make sure every graph has unique name, name goes down to the second
    current_date = str(datetime.today()).strip().replace(" ", "_").replace("-", "_").replace(":", "_")
    plt.savefig("output/bar_ylog_" + current_date[0: current_date.index(".")] + ".png")


def plot_line(timing, sizes, iterations):
    plt.plot(sizes, timing)
    plt.xlabel("Length of Sequence")
    plt.ylabel("Average run time (in seconds) over " + str(iterations) + " total iterations")
    plt.title("Stochastic Context Free Grammar Base Pairing Run-Time")
    plt.xticks(sizes)

    # save an image of the graph to the output folder
    # using the date to make sure every graph has unique name, name goes down to the second
    current_date = str(datetime.today()).strip().replace(" ", "_").replace("-", "_").replace(":", "_")
    plt.savefig("output/line_" + current_date[0: current_date.index(".")] + ".png")


def plot_line_ylog(timing, sizes, iterations):
    plt.plot(sizes, timing)
    plt.xlabel("Length of Sequence")
    plt.ylabel("Average run time (in seconds, log10) over " + str(iterations) + " total iterations")
    plt.title("Stochastic Context Free Grammar Base Pairing Run-Time")
    plt.xticks(sizes)
    plt.yscale("log")

    # save an image of the graph to the output folder
    # using the date to make sure every graph has unique name, name goes down to the second
    current_date = str(datetime.today()).strip().replace(" ", "_").replace("-", "_").replace(":", "_")
    plt.savefig("output/line_ylog_" + current_date[0: current_date.index(".")] + ".png")


if __name__ == "__main__":
    start = timeit.default_timer()
    main()
    end = timeit.default_timer()
    print(str(end-start))

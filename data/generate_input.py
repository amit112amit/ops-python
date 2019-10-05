"""
Create input files for melting simulations.
"""

import pandas as pd
import numpy as np
import itertools

ALPHA = 250000
DELTA = 15
VITERMAX = 2000000
PRINTSTEP = 2000

# Set up area interpolation
DATA = pd.read_csv('T7_OPS_Asphericity.dat', sep='\t').sort_values('Gamma')
GAMMAVALS = DATA['Gamma'].to_numpy()
AREAVALS = DATA['Area'].to_numpy()

# The input file HEADER line
HEADER = 'Alpha\tBeta\tGamma\tPercentStrain\tAreaConstr\tViterMax\tPrintStep\n'


def write_input_file(temp, gamma, index):
    """
    temp: temperature
    gamma: Gamma or the FvK number
    index: Index number to be suffixed to the file name i.e. schedule-index.dat
    """
    area = np.interp(gamma, GAMMAVALS, AREAVALS)
    full = '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:d}'.format(
        ALPHA, 1.0 / temp, gamma, DELTA, area, VITERMAX, PRINTSTEP)
    file_name = 'schedule-{0}.dat'.format(index)
    with open(file_name, 'w') as fobj:
        fobj.write(HEADER)
        fobj.write(full)


if __name__ == '__main__':

    """
    The temperature scheme is as follows:
        i. If γ < 10^2, 1/β ∈ (0.0, 2.0) --> linearly spaced
        ii. If 10^2 < γ < 10^3, 1/β ∈ (0.0, 1.0) --> linearly spaced
        iii. If γ > 10^3, 1/β ∈ (0.001, 0.25) --> log spaced
    """
    numgamma = 30
    numtemp = 20

    gamma = 2.0 * np.logspace(-1, 4, numgamma)
    temperature = np.empty((numgamma, numtemp))
    
    for g, t in zip(gamma, temperature):
        if g < 1e2:
            t[:] = np.linspace(0.1, 2.0, 20)
        elif g < 1e3:
            t[:] = np.linspace(0.05, 1.0, 20)
        else:
            lowtemp = 0.001
            hightemp = 0.25
            t[:] = np.exp(np.linspace(np.log(lowtemp), np.log(hightemp), 20))

    index = itertools.count(1)

    for i, g in enumerate(gamma):
        for t in temperature[i]:
            write_input_file(t, g, next(index))

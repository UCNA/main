# Author: Simon Slutsky

# Quickie script to make correlation plots between variables of linearity fits
import numpy as np
import matplotlib.pyplot as plt
#from ROOT import TCanvas

if __name__ == "__main__":

    plt.ion() #turn on interactive mode

    # import data. 
    data = np.genfromtxt("imagesNoInitialize/FitResults.txt", skip_header=1, delimiter = "\t", names = ['Run','Channel', 'Wavelength','p0','p0Err','p1','p1Err','p2','p2Err', 'Chi2'])


    cutChannel = data['Channel'] == 3
    cutWave = data['Wavelength'] == 405
    cutCond = cutChannel & cutWave
    data_cut = data[cutCond]

    x = data_cut['p0']
    y = data_cut['p1']
    z = data_cut['p2']

    figxy, axxy = plt.subplots()
    figxz, axxz = plt.subplots()
    figyz, axyz = plt.subplots()
    axxy.scatter(x,y)
    axxz.scatter(x,z)
    axyz.scatter(y,z)

    plt.show(block=True)

# pd_led_linearity.py
# Author: Simon Slutsky
# Created: 10/01/2013
# 
# Plot fitted linearity constants output from pd_led_pmt_batch file 
# Useful syntax: carp = data[data['Channel'] == 4]


import numpy as np
import matplotlib.pyplot as plt
#from ROOT import TCanvas

if __name__ == "__main__":

    plt.ion() #turn on interactive mode

    # import data. 
    # TODO: improve this to read the header string from the file
    data = np.genfromtxt("images2ndorder/FitResults.txt", skip_header=1, delimiter = "\t", names = ['Run','Channel', 'Wavelength','p0','p0Err','p1','p1Err', 'Chi2'])

    fitpars = ['p0', 'p1']
    
    run = data['Run']
    data_cut_405 = list()
    data_cut_465 = list()
    data_cut_err_405 = list()
    data_cut_err_465 = list()
    for wave in (405, 465):
        for chan in range(0,8): # 8 channels
            cutChannel = data['Channel'] == chan
            cutWave = data['Wavelength'] == wave
            cutCond = cutChannel & cutWave
            _data_cut = data[cutCond]
            if wave == 405:
                data_cut_405.append(_data_cut)
            if wave == 465:
                data_cut_465.append(_data_cut)
    

    #fig = plt.figure()
    fig0, (ax0, ax1) = plt.subplots(nrows=2)
    fig1, (ax2, ax3) = plt.subplots(nrows=2)
    fig2, (ax4, ax5) = plt.subplots(nrows=2)
    figChi, (ax6, ax7) = plt.subplots(nrows=2)

    plt.rc('axes', color_cycle=['r', 'g', 'b', 'y'])
    #plt.rc('axes', color_cycle=['c', 'm', 'y', 'k']) # looks terrible
    #ax0.set_color_cycle(['c', 'm', 'y', 'k'])
    #ax1.set_color_cycle(['r', 'g', 'b', 'y'])
    
    ax0.set_title("p0   405nm")
    ax1.set_title("p0   465nm")    
    ax2.set_title("p1   405nm")
    ax3.set_title("p1   465nm")    
    ax4.set_title("p2   405nm")
    ax5.set_title("p2   465nm")    

    ax1.set_xlabel('Run Number')
    ax3.set_xlabel('Run Number')
    ax5.set_xlabel('Run Number')
    #ax1.set_ylabel('your y label...')

    chan = "Ch. "
    _chan = ""
    for i in range(0,8):
        _chan = chan + str(i)
        # 405nm
        ax0.errorbar(data_cut_405[i]['Run'], data_cut_405[i]['p0'], yerr=data_cut_405[i]['p0Err'],
                 linestyle='None', marker='o', markersize=4, label=_chan)
        ax2.errorbar(data_cut_405[i]['Run'], data_cut_405[i]['p1'], yerr=data_cut_405[i]['p1Err'],
                 linestyle='None', marker='o', markersize=4, label=_chan)
     #   ax4.errorbar(data_cut_405[i]['Run'], data_cut_405[i]['p2'], yerr=data_cut_405[i]['p2Err'],
#                 linestyle='None', marker='o', markersize=4, label=_chan) 
        ax6.plot(data_cut_405[i]['Run'], data_cut_405[i]['Chi2'], 
                 linestyle='None', marker='o', markersize=4, label=_chan)
        # 465nm
        ax1.errorbar(data_cut_465[i]['Run'], data_cut_465[i]['p0'], yerr=data_cut_465[i]['p0Err'],
                     linestyle='None', marker='o', markersize=4, label=_chan)
        ax3.errorbar(data_cut_465[i]['Run'], data_cut_465[i]['p1'], yerr=data_cut_465[i]['p1Err'], 
                 linestyle='None', marker='o', markersize=4, label=_chan)
    #    ax5.errorbar(data_cut_465[i]['Run'], data_cut_465[i]['p2'], yerr=data_cut_465[i]['p2Err'],
 #                linestyle='None', marker='o', markersize=4, label=_chan)
        ax7.plot(data_cut_465[i]['Run'], data_cut_465[i]['Chi2'], 
                 linestyle='None', marker='o', markersize=4, label=_chan)


    
    leg0 = ax0.legend()
    leg1 = ax1.legend()
    leg2 = ax2.legend()
    leg3 = ax3.legend()
    leg4 = ax4.legend()
    leg5 = ax5.legend()
    leg6 = ax6.legend()
    leg7 = ax7.legend()
    
    ax0.set_xlim([20800, 24000])
    ax1.set_xlim([20800, 24000])
    ax2.set_xlim([20800, 24000])
    ax3.set_xlim([20800, 24000])
    ax4.set_xlim([20800, 24000])
    ax5.set_xlim([20800, 24000])

    ax0.set_ylim([-10, 50])
    ax1.set_ylim([-0, 10])
    ax2.set_ylim([-0.15, 0.075])
    ax3.set_ylim([-0.015, 0.005])
    ax4.set_ylim([-1e-4, 5e-4])
    ax5.set_ylim([-1e-6, 1e-5])
    plt.show(block=True)     #block=True keeps the plot window open when in interactive mode 



# too hard to do it manually

#run = data['Run']
#_run0 = np.where(data['Channel'] == 0, run, 0]
#run0 = _run0[_run0.nonzero() #strip zeroes from the array

#p0 = data['p0']
#_p00 = np.where(data['Channel'] == 0, p0, 0]
#p00 = _p00[_p00.nonzero()

#for i in range (0, len(run)):
#    if y[i] > 1e8 or y[i] < -1e2:
#        x[i] = 0
#        y[i] = 0
#        print i

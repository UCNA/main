# dirty script to count number of repeated runs in fit file
import os


if __name__ == "__main__":
    beginchar = "grep -cP '"
#    endchar = "\t' /data1/saslutsky/LEDPulser/images_06_23_2015_16way_separate_wavelength_coeff_16306_19316/FitResults_Combined.txt"
    endchar = "\t' /data1/saslutsky/LEDPulser/images_06_10_2015_16way_separate_wavelength_coeff_20254_23173/FitResults_Combined.txt"

#    for run in range (16306, 19317):
    for run in range (20254, 23173):
        os.system(beginchar + str(run) + endchar)
        
    
        

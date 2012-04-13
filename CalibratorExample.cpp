#include "CalDBSQL.hh"
#include "EnergyCalibrator.hh"
#include "PMTGenerator.hh"
#include <stdio.h>

/// an example of how to use the energy calibration classes
int main(int argc, char *argv[]) {
	
	// get pointer to default calibration DB
	CalDBSQL* CDB = CalDBSQL::getCDB();
	// calibrator for run number 15931
	RunNum rn = 15931;
	PMTCalibrator PCal(rn,CDB);
	
	// we can get lots of information about the run calibrations from the PMTCalibrator (see EnergyCalibrator.hh header):
	printf("Run %i: PMT W2 sees %.1f PE for 500keV quenched energy at position (5,8)\n",rn,PCal.nPE(WEST, 1, 500, 5., 8.));
	printf("Run %i: Total of %.1f PE/MeV at East position (7,-4.3)\n",rn,PCal.nPE(EAST, nBetaTubes, 1000, 7., -4.3));
	printf("The light transport factor for PMT E3 at (-12,33) is eta=%.2f\n",PCal.eta(EAST, 2, -12., 33.));
	printf("Run %i: an ADC reading of 145 for PMT W1 at (-6,-7) corresponds to Equenched=%.1f\n\n",rn,PCal.calibratedEnergy(WEST, 0, -6., -7., 145.).x);
	
	// object for simulating PMT response
	PMTGenerator PGen;
	// use previously defined PMTCalibrator for response
	PGen.setCalibrator(&PCal);
	// set side to generate response at
	PGen.setSide(EAST);
	// set position to generate response at (decay tube coordinates x,y in mm)
	PGen.setPosition(15., -7.3);
	// advanced usage: set actual position (15,-7.3) in scintillator and offset to (17.1,-6.2) in wirechamber
	PGen.setPosition(15., -7.3, 2.1, 1.1);
	
	// generate events for a range of energies
	for(float eq = 10; eq < 500; eq += 10) {
		ScintEvent evt = PGen.generate(eq);
		printf("Eq=%gkeV:\tADCs = (%.1f,\t%.1f,\t%.1f,\t%.1f),\tEobs = %.1f",
			   eq,evt.adc[0],evt.adc[1],evt.adc[2],evt.adc[3],evt.energy.x);
		if(PGen.triggered())
			printf(",\ttriggered.\n");
		else
			printf(",\tNO 2-of-4 trigger.\n");
	}
	
	
	return 0;
}
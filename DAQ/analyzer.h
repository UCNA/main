/********************************************************************\

  Name:         analyzer.h
  Created by:   Stefan Ritt

  Contents:     Analyzer global include file

  $Log: analyzer.h,v $
  Revision 1.2  1998/10/12 12:18:58  midas
  Added Log tag in header


\********************************************************************/
                                                        
/*-- Parameters ----------------------------------------------------*/

/* number of channels */
#define N_ADC              32 
#define N_TDC              32 
#define N_SCLR              32

/*-- Histo ID bases ------------------------------------------------*/

#define ADCCALIB_ID_BASE 1000
#define ADCSUM_ID_BASE   2000
#define PDCCALIB_ID_BASE 3000
#define S830_ID_BASE    4000
#define TDCCALIB_ID_BASE 5000
#define PDC2CALIB_ID_BASE 6000
#define PDC3CALIB_ID_BASE 8000
//#define ACCUM_ID_BASE   7000

double blinding_factor1;
double blinding_factor2;



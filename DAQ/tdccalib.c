/********************************************************************\

  Name:         tdccalib.c
  Created by:   Stefan Ritt

  Contents:     Example analyzer module for ADC calibration. Looks
                for ADC0 bank, subtracts pedestals and applies gain
                calibration. The resulting values are appended to 
                the event as an CADC bank ("calibrated ADC"). The
                pedestal values and software gains are stored in
                adccalib_param structure which was defined in the ODB
                and transferred to experim.h.

  $Log: adccalib.c,v $
  Revision 1.3  2002/05/10 05:22:47  pierre
  add MANA_LITE #ifdef

  Revision 1.2  1998/10/12 12:18:58  midas
  Added Log tag in header


\********************************************************************/
                                                    
/*-- Include files -------------------------------------------------*/

/* standard includes */
#include <stdio.h>
#include <time.h>

/* midas includes */
#include "midas.h"
#include "experim.h"
#include "analyzer.h"

/* cernlib includes */
#ifdef OS_WINNT
#define VISUAL_CPLUSPLUS
#endif

#ifdef __linux__
#define f2cFortran
#endif

#ifndef MANA_LITE
#include <cfortran.h>
#include <hbook.h>
#endif

/*-- Parameters ----------------------------------------------------*/

TDC_CALIBRATION_PARAM tdccalib_param;
extern EXP_PARAM      exp_param;
extern RUNINFO        runinfo;

/*-- Module declaration --------------------------------------------*/

INT tdc_calib(EVENT_HEADER*,void*);
INT tdc_calib_init(void);
INT tdc_calib_bor(INT run_number);
INT tdc_calib_eor(INT run_number);

TDC_CALIBRATION_PARAM_STR(tdc_calibration_param_str);

ANA_MODULE tdc_calib_module = {
  "TDC calibration",             /* module name           */  
  "Stefan Ritt",                 /* author                */
  tdc_calib,                     /* event routine         */
  tdc_calib_bor,                 /* BOR routine           */
  tdc_calib_eor,                 /* EOR routine           */
  tdc_calib_init,                /* init routine          */
  NULL,                          /* exit routine          */
  &tdccalib_param,               /* parameter structure   */
  sizeof(tdccalib_param),        /* structure size        */
  tdc_calibration_param_str,     /* initial parameters    */
};

/*-- init routine --------------------------------------------------*/

INT tdc_calib_init(void)
{
  /* book histos */
  tdc_calib_bor(0);
  return SUCCESS;
}

/*-- BOR routine ---------------------------------------------------*/

#define TDC_N_BINS  1028
#define TDC_X_LOW      -0.5
#define TDC_X_HIGH  4095.5

INT tdc_calib_bor(INT run_number)
{
int    i;
char   str[80];

#ifdef MANA_LITE
 printf("manalite: tdc_calib_bor: HBOOK disable\n");
#else
  /* book TDC histos */
/*
  for (i=0; i<5; i++)
    {
    if (HEXIST(TDCCALIB_ID_BASE+i))
      HDELET(TDCCALIB_ID_BASE+i);
    sprintf(str, "TDC%02d", i);
    HBOOK1(TDCCALIB_ID_BASE+i, str, TDC_N_BINS, 
           (float)TDC_X_LOW, (float)TDC_X_HIGH, 0.f); 
    }
*/
#endif

  return SUCCESS;

}

/*-- eor routine ---------------------------------------------------*/

INT tdc_calib_eor(INT run_number)
{
  return SUCCESS;
}

/*-- event routine -------------------------------------------------*/

INT tdc_calib(EVENT_HEADER *pheader, void *pevent)
{
  WORD    i, n_tdc,ch;
  DWORD   *blockdata,*pdata,nch,tdc_en;
  float  *cadc;

  /* look for TDC bank, return if not present */
 /* 
  n_tdc = bk_locate(pevent, "R_TD", &blockdata);
  if (n_tdc == 0 || n_tdc > N_TDC+2)
    return 1;

  bk_create(pevent, "TDC0", TID_DWORD, &pdata);

  for(i=0; i<N_TDC; i++) 
    pdata[i]=0;

  nch=(blockdata[0] & 0x03F00)>>8;

  for(i=1; i<=nch; i++) {
    ch=(blockdata[i] & 0x3F0000)>>16;
    pdata[ch]=blockdata[i] & 0x0FFF;
  }    

  for(i=0; i<5; i++) {
      HF1(TDCCALIB_ID_BASE+i, (float)pdata[i], 1.f);
  }

  bk_close(pevent,pdata+N_TDC);
  
*/
  /* create calibrated ADC bank */
  //bk_create(pevent, "CADC", TID_FLOAT, &cadc);

  /* zero cadc bank */
  //for (i=0 ; i<N_ADC ; i++)
  //cadc[i] = 0.f;

  /* subtract pedestal */
  //for (i=0 ; i<n_adc ; i++)
  //cadc[i] = (float) ((double)pdata[i] - adccalib_param.pedestal[i] + 0.5);

  /* apply software gain calibration */
  //for (i=0 ; i<n_adc ; i++)
  //cadc[i] *= adccalib_param.software_gain[i];

#ifdef MANA_LITE
  //printf("manalite: adc_calib: HBOOK disable\n");
#else
  /* fill ADC histos if above threshold */
  //for (i=0 ; i<n_adc ; i++)
    //if (cadc[i] > (float) adccalib_param.histo_threshold)
  // HF1(ADCCALIB_ID_BASE+i, cadc[i], 1.f);
#endif

  /* close calculated bank */
  //bk_close(pevent, cadc+n_adc);

  return SUCCESS;
}

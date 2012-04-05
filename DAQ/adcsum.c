/********************************************************************\

  Name:         adcsum.c
  Created by:   Stefan Ritt

  Contents:     Example analyzer module for ADC summing. This module
                looks for a bank named CADC and produces an ASUM
                bank which contains the sum of all ADC values. The
                ASUM bank is a "structured" bank. It has been defined
                in the ODB and transferred to experim.h.

  $Log: adcsum.c,v $
  Revision 1.3  2002/05/10 05:22:59  pierre
  add MANA_LITE #ifdef

  Revision 1.2  1998/10/12 12:18:58  midas
  Added Log tag in header


\********************************************************************/
                                                        
/*-- Include files -------------------------------------------------*/

/* standard includes */
#include <stdio.h>
#include <math.h>

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

#ifndef PI
#define PI 3.14159265359
#endif

/*-- Parameters ----------------------------------------------------*/

ADC_SUMMING_PARAM adc_summing_param;

/*-- Module declaration --------------------------------------------*/

INT adc_summing(EVENT_HEADER*,void*);
INT adc_summing_init(void);
INT adc_summing_bor(INT run_number);

ADC_SUMMING_PARAM_STR(adc_summing_param_str);

ANA_MODULE adc_summing_module = {
  "ADC summing",                 /* module name           */  
  "Stefan Ritt",                 /* author                */
  adc_summing,                   /* event routine         */
  NULL,                          /* BOR routine           */
  NULL,                          /* EOR routine           */
  adc_summing_init,              /* init routine          */
  NULL,                          /* exit routine          */
  &adc_summing_param,            /* parameter structure   */
  sizeof(adc_summing_param),     /* structure size        */
  adc_summing_param_str,         /* initial parameters    */
};

/*-- init routine --------------------------------------------------*/

INT adc_summing_init(void)
{
#ifdef MANA_LITE
  printf("manalite: adc_summing_init: HBOOK disable\n");
#else
  /* book sum histo */
  HBOOK1(ADCSUM_ID_BASE, "ADC sum", 500, 0.f, 10000.f, 0.f); 
#endif
  return SUCCESS;
}

/*-- event routine -------------------------------------------------*/

INT adc_summing(EVENT_HEADER *pheader, void *pevent)
{
  INT          i, n_adc;
  float        *cadc;
  ASUM_BANK    *asum;

  /* look for CADC bank, return if not present */
  n_adc = bk_locate(pevent, "CADC", &cadc);
  if (n_adc == 0)
    return 1;

  /* create ADC sum bank */
  bk_create(pevent, "ASUM", TID_STRUCT, &asum);

  /* sum all channels above threshold */
  asum->sum = 0.f;
  for (i=0 ; i<n_adc ; i++)
    asum->sum += cadc[i];

#ifdef MANA_LITE
 printf("manalite: adc_summing: HBOOK disable\n");
#else
  /* fill sum histo */
  HF1(ADCSUM_ID_BASE, asum->sum, 1.f);
#endif

  /* close calculated bank */
  bk_close(pevent, asum+1);
  
  return SUCCESS;
}

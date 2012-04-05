/*   Name:         adccalib.c */
/*   Created by:   Stefan Ritt */

/*   Contents:     Example analyzer module for ADC calibration. Looks */
/*                 for ADC0 bank, subtracts pedestals and applies gain */
/*                 calibration. The resulting values are appended to  */
/*                 the event as an CADC bank ("calibrated ADC"). The */
/*                 pedestal values and software gains are stored in */
/*                 adccalib_param structure which was defined in the ODB */
/*                 and transferred to experim.h. */

/*   $Log: adccalib.c,v $
/*   Revision 1.1.1.1  2009/01/20 18:43:29  jliu
/*   Initial import of the online module.
/* */
/*   Revision 1.3  2002/05/10 05:22:47  pierre */
/*   add MANA_LITE #ifdef */

/*   Revision 1.2  1998/10/12 12:18:58  midas */
/*   Added Log tag in header */



                                                    
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

ADC_CALIBRATION_PARAM adccalib_param;
extern EXP_PARAM      exp_param;
extern RUNINFO        runinfo;
extern GLOBAL_PARAM   global_param;

FILE   *file_rate;

/*-- Module declaration --------------------------------------------*/

enum Side {
	EAST = 0,
	WEST = 1
};


INT adc_calib(EVENT_HEADER*,void*);
INT adc_calib_init(void);
INT adc_calib_bor(INT run_number);
INT adc_calib_eor(INT run_number);

static long clk_start = 0;
long clk_of_this_trig = 0;
long clk_of_previous_trig = 0;
static long scaler_saved[N_SCLR];
static long counts_since_last_calc[N_SCLR];
static long trigger_counter_last[2]; //east, west
static long trigger_counter[2]; //east, west

float rates_since_last_calc[N_SCLR];
time_t now;
char   str[256];

const double anode_cut_e = 250;
const double anode_cut_w = 600;

void remove_isolated_wires(float* wireADC)
{
  if(wireADC[0]>0&&wireADC[1]==0)
    wireADC[0]=0;
  int ii;
  for(ii=1;ii<=14;ii++) {
    if(wireADC[ii]>0&&wireADC[ii-1]==0
       &&wireADC[ii+1]==0)
      wireADC[ii]=0;
  }
  if(wireADC[15]>0&&wireADC[14]==0)
    wireADC[15]=0;
  return;
}



ADC_CALIBRATION_PARAM_STR(adc_calibration_param_str);

ANA_MODULE adc_calib_module = {
  "ADC calibration",             /* module name           */  
  "Stefan Ritt",                 /* author                */
  adc_calib,                     /* event routine         */
  adc_calib_bor,                 /* BOR routine           */
  adc_calib_eor,                 /* EOR routine           */
  adc_calib_init,                /* init routine          */
  NULL,                          /* exit routine          */
  &adccalib_param,               /* parameter structure   */
  sizeof(adccalib_param),        /* structure size        */
  adc_calibration_param_str,     /* initial parameters    */
};

/*-- init routine --------------------------------------------------*/

INT adc_calib_init(void)
{
  /* book histos */
  adc_calib_bor(0);
  return SUCCESS;
}

/*-- BOR routine ---------------------------------------------------*/

#define ADC_N_BINS  4096
#define ADC_X_LOW      -0.5
#define ADC_X_HIGH  4095.5
#define TDC_N_BINS  1028
#define TDC_X_LOW      -0.5
#define TDC_X_HIGH  4095.5


//My stupid pedestal finding routine based on my limited hbook knowledge.
//There might be a better way to do this. Jianglai 07-30-2008
float get_ped(const int hid)
{
  float content[ADC_N_BINS];
  HUNPAK(hid, content, "HIST", 0);
  
  int ii;
  int iimax = 1;
  float maxcontent = 0;
  //don't include 1st bin and 2nd half and the histogram in pedestal finding
  for(ii=2;ii<=ADC_N_BINS/2;ii++){
    //printf("%d\t%f\n",ii,content[ii]);
    if(content[ii]>maxcontent) {
      iimax = ii;
      maxcontent = content[ii];
    }
  }
  float x1, x2, y1, y2, ave;
  HIJXY(hid,iimax,1,x1,y1);
  HIJXY(hid,iimax+1,1,x2,y2);
  ave = 0.5*(x1+x2);
  //printf("%f %f %f\n", x1, x2, ave);
  return ave;
}

//float pmt1,pmt2,pmt3,pmt4;

INT adc_calib_bor(INT run_number)
{
int    i;
char   str[80];

#ifdef MANA_LITE
 printf("manalite: adc_calib_bor: HBOOK disable\n");
#else

 //HBOOK1(9, "He Det #1", 10000,0,100.,0.f);
 //HBOOK1(99, "UCN detector", 10000,0,100.,0.f);
 //HBOOK1(96, "Monitor 1", 10000,0,100.,0.f);
 //HBOOK1(97, "Monitor 2", 10000,0,100.,0.f);
 //HBOOK1(98, "UW detector", 10000,0,100.,0.f);
 //HBOOK1(999, "NaI GMS detector", 1000,0,2.,0.f);
 //HBOOK1(9900, "Flapper moving", 1000,0,21.,0.f);
  //HBOOK1(60, "Cathode Multiplicity (Horizontal wires)", 15,-.5,14.5,0.f);
  //HBOOK1(61, "Cathode Multiplicity (Vertical wires)", 15,-.5,14.5,0.f);

  /* book UCN monitor timing histograms */

 HBOOK1(4001, "UCN Monitor 1 Timing", 3000, 0., 300., 0.f);
 HBOOK1(4002, "UCN Monitor 2 Timing", 3000, 0., 300., 0.f);
 HBOOK1(4003, "UCN Monitor 3 Timing", 3000, 0., 300., 0.f);
 HBOOK1(4004, "UCN Monitor 4 Timing", 3000, 0., 300., 0.f);
 HBOOK1(4005, "UCN Monitor 5 Timing", 3000, 0., 300., 0.f);
 HBOOK1(4006, "UCN Monitor 6 Timing", 3000, 0., 300., 0.f);
 HBOOK1(7000, "UCN Iron Foil Timing", 840, 0., 4200., 0.f);

 HBOOK1(9991, "EAST QADC SUM ANODE", 16000, 0., 16000., 0.f);
 HBOOK1(9992, "WEST QADC SUM ANODE", 16000, 0., 16000., 0.f);
 HBOOK1(9993, "EAST TWO-FOLD TIMING", 180, 0., 180., 0.f);
 HBOOK1(9994, "WEST TWO-FOLD TIMING", 180, 0., 180., 0.f);

 HBOOK1(9995, "low Sn",  16000, 0., 16000., 0.f);
 HBOOK1(9996, "high Sn", 16000, 0., 16000., 0.f);
 HBOOK1(9997, "Sr",      16000, 0., 16000., 0.f);
 HBOOK1(9998, "Bi",      16000, 0., 16000., 0.f);

 HBOOK1(1077, "TOTAL QADC SUM chris pulser ped", 6000, 0., 6000., 0.f);

 HBOOK1(4008, "SWITCHER MONITOR ABSOLUTE TIMING", 250, 0., 250., 0.f);

  /* book UCN monitor PADC histograms */

 HBOOK1(4011, "UCN Monitor 1 PADC", ADC_N_BINS,
        (float)ADC_X_LOW, (float)ADC_X_HIGH, 0.f);
 HBOOK1(4012, "UCN Monitor 2 PADC", ADC_N_BINS,
        (float)ADC_X_LOW, (float)ADC_X_HIGH, 0.f);
 HBOOK1(4013, "UCN Monitor 3 PADC", ADC_N_BINS,
        (float)ADC_X_LOW, (float)ADC_X_HIGH, 0.f);
 HBOOK1(4014, "UCN Monitor 4 PADC", ADC_N_BINS,
        (float)ADC_X_LOW, (float)ADC_X_HIGH, 0.f);
 HBOOK1(4015, "UCN Monitor 5 PADC", ADC_N_BINS,
        (float)ADC_X_LOW, (float)ADC_X_HIGH, 0.f);
 HBOOK1(4016, "TAC channel 5 ", ADC_N_BINS,
        (float)ADC_X_LOW, (float)ADC_X_HIGH, 0.f);
 HBOOK2(4020,"UCN mon 5 ADC vs TAC",ADC_N_BINS/10.,
        (float)ADC_X_LOW, (float)ADC_X_HIGH, ADC_N_BINS/10.,
        (float)ADC_X_LOW, (float)ADC_X_HIGH, 0.f);

 HBOOK1(7777, "East Ref PADC", 100, (float)ADC_X_LOW, (float)ADC_X_HIGH, 0.f);
 
  /* book QADC histos */
  //for (i=0; i<N_ADC; i++)
 for(i=0; i<32; i++)
    {
    if (HEXIST(ADCCALIB_ID_BASE+i))
      HDELET(ADCCALIB_ID_BASE+i);
    sprintf(str, "QADC%02d", i);
    HBOOK1(ADCCALIB_ID_BASE+i, str, ADC_N_BINS, 
           (float)ADC_X_LOW, (float)ADC_X_HIGH, 0.f);
    }
 
/* Beta PMT QADCs with anode cut */
for(i=0; i<8; i++)
    {
    if (HEXIST(ADCCALIB_ID_BASE+900+i))
      HDELET(ADCCALIB_ID_BASE+900+i);
    sprintf(str, "QADC_ANODECUT_%02d", i);
    HBOOK1(ADCCALIB_ID_BASE+900+i, str, ADC_N_BINS, 
           (float)ADC_X_LOW, (float)ADC_X_HIGH, 0.f); 
    }

  /* book QADC histos for time after 1ms*/
  //for (i=0; i<N_ADC; i++)
 for(i=0; i<32; i++) {
    if (HEXIST(ADCCALIB_ID_BASE+100+i))
      HDELET(ADCCALIB_ID_BASE+100+i);
    sprintf(str, "QADC_1_%02d", i);
    HBOOK1(ADCCALIB_ID_BASE+100+i, str, ADC_N_BINS, 
           (float)ADC_X_LOW, (float)ADC_X_HIGH, 0.f); 
    
}

 if(HEXIST(ADCCALIB_ID_BASE+99))
   HDELET(ADCCALIB_ID_BASE+99);
 HBOOK1(ADCCALIB_ID_BASE+99, "QADCSUM_E", 4*ADC_N_BINS,
	(float)ADC_X_LOW, (float)(4*ADC_X_HIGH+1.5), 0.f);

 if(HEXIST(ADCCALIB_ID_BASE+199))
   HDELET(ADCCALIB_ID_BASE+199);
 HBOOK1(ADCCALIB_ID_BASE+199, "QADCSUM_1", 4*ADC_N_BINS,
	(float)ADC_X_LOW, (float)(4*ADC_X_HIGH+1.5), 0.f);

 if(HEXIST(ADCCALIB_ID_BASE+88))
   HDELET(ADCCALIB_ID_BASE+88);
 HBOOK1(ADCCALIB_ID_BASE+88, "QADCSUM_W", 4*ADC_N_BINS,
	(float)ADC_X_LOW, (float)(4*ADC_X_HIGH+1.5), 0.f);


 if (HEXIST(ADCCALIB_ID_BASE+200)) {
   HDELET(ADCCALIB_ID_BASE+200);
 }
 HBOOK1(ADCCALIB_ID_BASE+200, "sis3600", 4000, -0.5, 3999.5, 0.f);
 
  /* book TDC histos */
 for (i=0; i<32; i++)
   {
     if (HEXIST(TDCCALIB_ID_BASE+i))
       HDELET(TDCCALIB_ID_BASE+i);
     sprintf(str, "TDC%02d", i);
     HBOOK1(TDCCALIB_ID_BASE+i, str, TDC_N_BINS,
	    (float)TDC_X_LOW, (float)TDC_X_HIGH, 0.f);
   }
 
 /* book PADC histos */
 for (i=0; i<32; i++)
   { 
     if (HEXIST(PDCCALIB_ID_BASE+i))
       HDELET(PDCCALIB_ID_BASE+i);
     sprintf(str, "PADC%02d", i);
     HBOOK1(PDCCALIB_ID_BASE+i, str, ADC_N_BINS,
	    (float)ADC_X_LOW, (float)ADC_X_HIGH, 0.f);
   }
 
 /*book PADC histos for time after 1ms*/
 for(i=0; i<32; i++) 
   {
     if(HEXIST(PDCCALIB_ID_BASE+100+i))
       HDELET(PDCCALIB_ID_BASE+100+i);
     sprintf(str, "PADC_1_%02d", i);
     HBOOK1(PDCCALIB_ID_BASE+100+i, str, ADC_N_BINS,
	    (float)ADC_X_LOW, (float)ADC_X_HIGH, 0.f);
   }

 /* book PADC2 histos */
 for (i=0; i<32; i++)
   { 
     if (HEXIST(PDC2CALIB_ID_BASE+i))
       HDELET(PDC2CALIB_ID_BASE+i);
     sprintf(str, "PADC2_%02d", i);
     HBOOK1(PDC2CALIB_ID_BASE+i, str, ADC_N_BINS,
	    (float)ADC_X_LOW, (float)ADC_X_HIGH, 0.f);
   }


 /* book PADC3 histos */
 for (i=0; i<32; i++) 
   {
     if (HEXIST(PDC3CALIB_ID_BASE+i))
       HDELET(PDC3CALIB_ID_BASE+i);
     sprintf(str, "PADC3_%02d", i);
     HBOOK1(PDC3CALIB_ID_BASE+i, str, ADC_N_BINS,
	    (float)ADC_X_LOW, (float)ADC_X_HIGH, 0.f);
   }
     
  /* book individual PMT histograms: LED and Bi events */

 HBOOK1(9001, "E1 LED", 4100, 0, 4100, 0.f);
 HBOOK1(9002, "E2 LED", 4100, 0, 4100, 0.f);
 HBOOK1(9003, "E3 LED", 4100, 0, 4100, 0.f);
 HBOOK1(9004, "E4 LED", 4100, 0, 4100, 0.f);
 HBOOK1(9005, "W1 LED", 4100, 0, 4100, 0.f);
 HBOOK1(9006, "W2 LED", 4100, 0, 4100, 0.f);
 HBOOK1(9007, "W3 LED", 4100, 0, 4100, 0.f);
 HBOOK1(9008, "W4 LED", 4100, 0, 4100, 0.f);

 HBOOK1(9101, "E1 Bi", 4100, 0, 4100, 0.f);
 HBOOK1(9102, "E2 Bi", 4100, 0, 4100, 0.f);
 HBOOK1(9103, "E3 Bi", 4100, 0, 4100, 0.f);
 HBOOK1(9104, "E4 Bi", 4100, 0, 4100, 0.f);
 HBOOK1(9105, "W1 Bi", 4100, 0, 4100, 0.f);
 HBOOK1(9106, "W2 Bi", 4100, 0, 4100, 0.f);
 HBOOK1(9107, "W3 Bi", 4100, 0, 4100, 0.f);
 HBOOK1(9108, "W4 Bi", 4100, 0, 4100, 0.f);

  //HBOOK2(3330, "Cath sum Vs. Anode", 100,0,4000,100,0,4000,0.f);
  HBOOK1(3300, "Cathode sum East", 4000, 0, 16000, 0.f);
  HBOOK1(3301, "Cathode Sum West", 4000, 0, 16000, 0.f);
  //HBOOK2(8888, "QDC4 sum Vs. QDC5", 100,0,2000,100,0,2000,0.f);
  //HBOOK2(9999, "QDC6 sum Vs. QDC7", 100,0,2000,100,0,2000,0.f);

  HBOOK1(60,"cathEX",120,-60.0,60.0,0.f);
  HBOOK1(61,"cathEY",120,-60.0,60.0,0.f);
  HBOOK2(66,"cathEY Vs. cathEX",160,-80,80,160,-80,80,0.f);

  HBOOK1(70,"cathWX",120,-60.0,60.0,0.f);
  HBOOK1(71,"cathWY",120,-60.0,60.0,0.f);
  HBOOK2(76,"cathWY Vs. cathWX",160,-80,80,160,-80,80,0.f);

//char nt_row[4][20]={"pmt1","pmt2","pmt3","pmt4"};
 //HBOOKN(ADCCALIB_ID_BASE+99,"QADC ntuple",4,"ONLN",5000,nt_row);
 //HBNT(ADCCALIB_ID_BASE+99,"QADC ntuple"," ");
 //HBNAME(ADCCALIB_ID_BASE+99,"",&pmt1,"pmt1:R");
 //HBNAME(ADCCALIB_ID_BASE+99,"",&pmt2,"pmt2:R");
 //HBNAME(ADCCALIB_ID_BASE+99,"",&pmt3,"pmt3:R");
 //HBNAME(ADCCALIB_ID_BASE+99,"",&pmt4,"pmt4:R");
#endif

  //initialize scaler counters
  clk_start = 0;
  clk_of_this_trig = 0;
  for(i=0;i<N_SCLR;i++){
    scaler_saved[i] = 0;
    counts_since_last_calc[i] = 0;
    rates_since_last_calc[i] = 0;
  }

  for(i=0;i<2;i++){
    trigger_counter_last[i]=trigger_counter[i]=0;
  }
  return SUCCESS;
}

/*-- eor routine ---------------------------------------------------*/

INT adc_calib_eor(INT run_number)
{
  return SUCCESS;
}

/*-- event routine -------------------------------------------------*/

INT adc_calib(EVENT_HEADER *pheader, void *pevent)
{
  WORD    i, n_adc,ch,m,n_tdc,n_pdc,n_s820,multiplicity,n_s830;
  DWORD   *pdata,adcsum,*blockdata,nch,*psis,qdc_en,tdc_en,pdc1_en,
    pdc2_en,pdc3_en,s820_en;
  
  DWORD   *flipper;
  DWORD   *deltaT;
  
  DWORD   *bankHeaderFooter, qdc_hf, tdc_hf, pdc1_hf, pdc2_hf, pdc3_hf;

  DWORD   delta_trigger_time;//difference in timing scaler counts
  //float fdata[4];
  DWORD trigger_time; 

  const float cathX[16]={59.02,51.15,43.28,35.41,27.54,19.67,11.80,3.93,-3.93,
			 -11.80,-19.67,-27.54,-35.41,-43.28,-51.15,-59.02};
  const float cathY[16]={59.02,51.15,43.28,35.41,27.54,19.67,11.80,3.93,-3.93,
			 -11.80,-19.67,-27.54,-35.41,-43.28,-51.15,-59.02};


  //Important: adc_calib is being called for every event!
  //In order to have all pedestals updated correctly, they have to 
  //defined as static variables here, otherwise even time we enter 
  //the subroutine, all peds are redefined to the initial values.
  static float pedE1[16]={700.,685.,785.,720.,725.,760.,690.,665.,880.,
			  565.,595.,760.,745.,750.,655.,610.}; // july 15, 2008
  static float pedE2[16]={780.,765.,760.,760.,735.,830.,735.,785.,975.,
			  765.,765.,775.,780.,830.,810.,730.}; // july 15, 2008
  
  //Jianglai implemented W MWPC realtime display. July 28, 2008
  
  static float pedW1[16] = {370, 843, 629, 727, 620, 757, 715, 667, 579,
			    666, 655, 739, 706, 754, 725, 618};
  /*   float pedW2[] = {753, 653, 776, 764, 595, 758, 831, 750, 716, */
  /* 		   725, 791, 647, 754, 725, 673, 714};  */
  
  static float pedW2[16] = {716, 668, 707, 703, 630, 709, 796, 719, 718, 
			    704, 703, 675, 726, 699, 675, 795};
  
  float cathE1[16],cathE2[16],cathE1max=0,cathE2max=0,posEX,posEY;
  float cathW1[16],cathW2[16],cathW1max=0,cathW2max=0,posWX,posWY;

  WORD i1max,i2max,*mult;
  float tempw=0,tempwx=0,anode,*position;
  
  float *clock;
  float time_since_last_calc_in_sec=0;
 
  float pmtQadc[N_ADC]; 
  float adcsum_east;
  float adcsum_west;
  float anode_east;
  float anode_west;
  float cathode_sum[2];		//< cathode sum on each side
  float mwpc_cut_var[2];	//< variable for doing wirechamber cuts
  float mwpc_cut_level[2];	//< cut level for mwpc cut
  float east_twofold;
  float west_twofold;

  //for electronics error checking
  int nheader = 0;
  int nfooter = 0;
  int nvalid = 0;  

  //for electronics error checking alarm
  static int header_footer_error_sofar=0;
  static int nevent_counter_error_sofar=0;

  //float parX[3]={1,1,1};
  //float parY[3]={1,1,1};
  //float step[3],sigpar[3],chi2,pmin[3],pmax[3];

  qdc_en=0;tdc_en=0;pdc1_en=0;pdc2_en=0;pdc3_en=0;s820_en=0;

/*   /\*look for S820 bank*\/ */
/*   n_s820=bk_locate(pevent,"S820",&pdata); */
/*   trigger_time=pdata[0]; */

  //Do V820, physically a V830. simply skip the header.
  n_s820 = bk_locate(pevent, "R_ST0", &blockdata);
  if (!(n_s820 == 0 || n_s820 > N_SCLR+1)) { //33 channels, including header
    bk_create(pevent, "S820", TID_DWORD, &pdata);

    //number of scaler channels in this event
    nch = (blockdata[0]>>18)&0x3F;
    //printf("%d\n",(blockdata[0]>>18)&0x3F);
  
    //reset all 32 words that will go into the ntuple
    for(i=0; i<nch; i++) 
      pdata[i]=0;
    
    for(i=1; i<=nch; i++) {//from 1 to 32
      ch = i-1; //skip the header word, read the next 32
      //printf("%d\t%d\n",i,blockdata[i]);
      pdata[ch]=blockdata[i]&0xFFFFFFFF;
    }
    //assign time here!!!
    trigger_time = pdata[0];
  } else {
    printf("Can not find R_ST0 bank from the raw data!!!\n");
    //use the time from the previous read then. Do nothing here
  }
  bk_close(pevent,pdata+N_SCLR);


  /* look for QADC bank, return if not present */
  n_adc = bk_locate(pevent, "R_QD", &blockdata);
  nheader = 0;
  nfooter = 0;
  nvalid = 0; 
  if (!(n_adc == 0 || n_adc > N_ADC+2)) {
    //return 1;

    bk_create(pevent, "QADC", TID_DWORD, &pdata);
    
    for(i=0; i<n_adc; i++) {
      DWORD wordID = (blockdata[i] >> 24) & 0x7;
      if(wordID == 0x2) ++nheader;
      if(wordID == 0x4) ++nfooter;
      if(wordID == 0x0) ++nvalid;
    } // for i

    for(i=0; i<N_ADC; i++) 
      pdata[i]=0;
    
    nch=(blockdata[0] & 0x03F00)>>8;
    
    for(i=1; i<=nch; i++) {
      ch=(blockdata[i] & 0x3F0000)>>16;
      if(blockdata[i] & 1<<12) pdata[ch]=4095; //overflow
    //else if(blockdata[i] & 1<<13) pdata[ch]=0; //underflow
      else
	pdata[ch]=blockdata[i] & 0x0FFF;
    }    
    qdc_en=blockdata[i] & 0xFFFFFF;
    

    for(i=0; i<N_ADC; i++) {
      HF1(ADCCALIB_ID_BASE+i, (float)pdata[i], 1.f);
      pmtQadc[i]=(float)pdata[i];
    }
    //east
    adcsum=pdata[0]+pdata[1]+pdata[2]+pdata[3];
    adcsum_east = (float)adcsum;
    HF1(ADCCALIB_ID_BASE+99, (float)adcsum, 1.f);
    //west
    adcsum=pdata[4]+pdata[5]+pdata[6]+pdata[7];
    adcsum_west = (float)adcsum;
    HF1(ADCCALIB_ID_BASE+88, (float)adcsum, 1.f);
    //HF2(8888,pdata[4],pdata[5],1.f); 
    //HF2(9999,pdata[6],pdata[7],1.f); 
   
    bk_close(pevent,pdata+N_ADC);  
  }

  qdc_hf = ((nheader<<4)&0xf0 | (nfooter&0xf));

  m = bk_locate(pevent, "SIS0", &psis);
  //bit 17 for flipper state
  unsigned int flipper_state = (psis[0]>>16)&0x1;
  unsigned int upper16 = (psis[0]>>16)&0xF;
  //printf("flipper state = %x\n",flipper_state);
  psis[0] = psis[0]&0xFFFF; //only store the lower 16 bits into the ntuple
  if(m!=0) //return 1;
    HF1(ADCCALIB_ID_BASE+200, psis[0], 1.f);

  if(psis[0]==260)   HF1(4001, (float)trigger_time/1.e6, 1.f);
  //if(psis[0]==516)  HF1(7000, (float)trigger_time/1.e6, 1.f);
  if(psis[0]==516)   HF1(4002, (float)trigger_time/1.e6, 1.f);
  if(psis[0]==1028)  HF1(4003, (float)trigger_time/1.e6, 1.f);
  if(psis[0]==2052) HF1(4004, (float)trigger_time/1.e6, 1.f); // for afp test 6/10/08
  if(psis[0]==64)  HF1(4005, (float)trigger_time/1.e6, 1.f);
  if(psis[0]==128) HF1(4006, (float)trigger_time/1.e6, 1.f);

  // Filling histograms with GMS information for the History function
  
  if(psis[0]==32) { 
    for(i = 0 ; i<4; i++){
      global_param.gms_value[i] +=(float)pdata[i];
      global_param.nevents[0]++;
    }
  }
  if(psis[0]==64) { 
     for(i = 0 ; i<4; i++){
       global_param.gms_value[i+4] += (float)pdata[i+4];
       global_param.nevents[1]++;
     }
  }
      
  //if(psis[0]==256) HF1(4006, (float)trigger_time/1.e6, 1.f);

  //if(psis[0]==8 || psis[0]==0) HF1(9900, (float)trigger_time/1.e6, 1.f);
  //if(psis[0]==16) {
  // HF1(99, (float)trigger_time/1.e6, 1.f); 
  //}
  //if(psis[0]==128) HF1(98, (float)trigger_time/1.e6, 1.f);
  //if(psis[0]==64) HF1(96, (float)trigger_time/1.e6, 1.f);
  //if(psis[0]==32) HF1(97, (float)trigger_time/1.e6, 1.f);
 
  //if(psis[0]==4) HF1(999, (float)trigger_time/1.208e6, 1.f);

  //Move the reading of S8300 to the last

  //implemented S Clayton's TDC error finder
  n_tdc = bk_locate(pevent, "R_TD", &blockdata);
  nheader = 0;
  nfooter = 0;
  nvalid = 0;  
  if (!(n_tdc == 0 || n_tdc > N_TDC+2)) {
    // return 1;

    bk_create(pevent, "TDC0", TID_DWORD, &pdata);
    
    for(i=0; i<n_tdc; i++) {
      DWORD wordID = (blockdata[i] >> 24) & 0x7;
      if(wordID == 0x2) ++nheader;
      if(wordID == 0x4) ++nfooter;
      if(wordID == 0x0) ++nvalid;
    } // for i

    for(i=0; i<N_TDC; i++)
      pdata[i]=0;
    
    nch=(blockdata[0] & 0x03F00)>>8;

    for(i=1; i<=nch; i++) {
      ch=(blockdata[i] & 0x3F0000)>>16;
      pdata[ch]=blockdata[i] & 0x0FFF;
    }
    tdc_en=blockdata[i] & 0xFFFFFF;
    
    
    for(i=0; i<32; i++) {
      HF1(TDCCALIB_ID_BASE+i, (float)pdata[i], 1.f);
    }
    east_twofold = (float)pdata[16]*180./4096.;
    west_twofold = (float)pdata[17]*180./4096.;
    
    
    bk_close(pevent,pdata+N_TDC);
  }

  tdc_hf = (nheader<<4)&0xf0 | (nfooter&0xf);

  
  n_pdc = bk_locate(pevent, "R_PD", &blockdata);
  nheader = 0;
  nfooter = 0;
  nvalid = 0;  
  if (!(n_pdc == 0 || n_pdc > N_ADC+2)) {
    // return 1;
    
    bk_create(pevent, "PADC", TID_DWORD, &pdata);

    for(i=0; i<n_pdc; i++) {
      DWORD wordID = (blockdata[i] >> 24) & 0x7;
      if(wordID == 0x2) ++nheader;
      if(wordID == 0x4) ++nfooter;
      if(wordID == 0x0) ++nvalid;
    } // for i
    
    for(i=0; i<N_ADC; i++)
      pdata[i]=0;

    nch=(blockdata[0] & 0x03F00)>>8;
    //printf("nch = %d\n",nch);
    for(i=1; i<=nch; i++) {
      ch=(blockdata[i] & 0x3F0000)>>16;
      if(blockdata[i] & 1<<12) pdata[ch]=4095; //overflow
      else if(blockdata[i] & 1<<13) pdata[ch]=0; //underflow
      else
	pdata[ch]=blockdata[i] & 0x0FFF;
      //printf("ch = %d, value = %d\n",ch, pdata[ch]);
    }
    pdc1_en=blockdata[i] & 0xFFFFFF;

    for(i=0; i<32; i++) 
      HF1(PDCCALIB_ID_BASE+i, (float)pdata[i], 1.f);


    //fill cathode w data
    for(i=0;i<16;i++) {
      cathW1[i]=pdata[i]-pedW1[i];
      if(pdata[i]>4050.) cathW1[i]=0;
    }
    remove_isolated_wires(cathW1);
    for(i=0;i<16;i++) {
      cathW2[i]=pdata[i+16]-pedW2[i];
      if(pdata[i+16]>4050.) cathW2[i]=0;
    } 
    remove_isolated_wires(cathW2);

    cathode_sum[WEST]=0;
    for(i=0; i<32; i++)  cathode_sum[WEST]+=pdata[i];
    HF1(3301,cathode_sum[WEST], 1.f);
    
    //HF2(3330,pdata[0],cathSum/28.,1.f);
    //HF2(8888,adcsum,pdata[0],1.f);
   
    bk_close(pevent,pdata+N_ADC);
  }

  pdc1_hf = (nheader<<4)&0xf0 | (nfooter&0xf);

  /*padc2*/
  n_pdc = bk_locate(pevent, "R_P2", &blockdata);
  nheader = 0;
  nfooter = 0;
  nvalid = 0;    
  if (!(n_pdc == 0 || n_pdc > N_ADC+2)) {
    // return 1;
    
    bk_create(pevent, "PDC2", TID_DWORD, &pdata);

    for(i=0; i<n_pdc; i++) {
      DWORD wordID = (blockdata[i] >> 24) & 0x7;
      if(wordID == 0x2) ++nheader;
      if(wordID == 0x4) ++nfooter;
      if(wordID == 0x0) ++nvalid;
    } // for i
    
    for(i=0; i<N_ADC; i++)
      pdata[i]=0;

    nch=(blockdata[0] & 0x03F00)>>8;
    
    for(i=1; i<=nch; i++) {
      ch=(blockdata[i] & 0x3F0000)>>16;
      if(blockdata[i] & 1<<12) pdata[ch]=4095; //overflow
      else if(blockdata[i] & 1<<13) pdata[ch]=0; //underflow
      else
	pdata[ch]=blockdata[i] & 0x0FFF;
    }

    pdc2_en=blockdata[i] & 0xFFFFFF;

    for(i=0; i<32; i++) 
      HF1(PDC2CALIB_ID_BASE+i, (float)pdata[i], 1.f);


    //HRESET(60,' '); HRESET(61,' ');
    for(i=0;i<16;i++) {
      cathE1[i]=pdata[i]-pedE1[i];
      if(pdata[i]>4050.) cathE1[i]=0;
      //HF1E(61,cathY[i],cathE1[i],50.);
    }
    remove_isolated_wires(cathE1);
    for(i=0;i<16;i++) {
      cathE2[i]=pdata[i+16]-pedE2[i];
      if(pdata[i+16]>4050.) cathE2[i]=0;
      //HF1E(60,cathX[i],cathE2[i],50.);
    }
    remove_isolated_wires(cathE2);


    cathode_sum[EAST]=0;
    for(i=0; i<32; i++)  cathode_sum[EAST] += pdata[i];
    HF1(3300,cathode_sum[EAST], 1.f);

    bk_close(pevent,pdata+N_ADC);
  }

  pdc2_hf = (nheader<<4)&0xf0 | (nfooter&0xf);

  /*padc3*/
  n_pdc = bk_locate(pevent, "R_P3", &blockdata);
  nheader = 0;
  nfooter = 0;
  nvalid = 0;    
  if (!(n_pdc == 0 || n_pdc > N_ADC+2)) {
    // return 1;
    
    bk_create(pevent, "PDC3", TID_DWORD, &pdata);

    for(i=0; i<n_pdc; i++) {
      DWORD wordID = (blockdata[i] >> 24) & 0x7;
      if(wordID == 0x2) ++nheader;
      if(wordID == 0x4) ++nfooter;
      if(wordID == 0x0) ++nvalid;
    } // for i
    
    for(i=0; i<N_ADC; i++)
      pdata[i]=0;
    
    nch=(blockdata[0] & 0x03F00)>>8;
    
    for(i=1; i<=nch; i++) {
      ch=(blockdata[i] & 0x3F0000)>>16;
      if(blockdata[i] & 1<<12) pdata[ch]=4095; //overflow
      else if(blockdata[i] & 1<<13) pdata[ch]=0; //underflow
      else
	pdata[ch]=blockdata[i] & 0x0FFF;
    }
    
    pdc3_en=blockdata[i] & 0xFFFFFF;
    
    for(i=0; i<32; i++){ 
      HF1(PDC3CALIB_ID_BASE+i, (float)pdata[i], 1.f);
      
      //fill 7777 histogram if this is an LED trigger
      //if(i==2&&psis[0]>127&&psis[0]<132) HF1(7777, (float)pdata[i], 1.f);
      if( i==2 && (psis[0]==129||psis[0]==130||psis[0]==131||psis[0]==161||psis[0]==162||psis[0]==163) ) HF1(7777, (float)pdata[i], 1.f);

    }

    /* fill UCN monitor PADC histograms for time after 1 ms */

    if (trigger_time>1000) {
      if (psis[0]== 260)   HF1(4011,  (float)pdata[8]  , 1.f);
      if (psis[0]== 516)   HF1(4012,  (float)pdata[9]  , 1.f);
      if (psis[0]== 1028)  HF1(4013,  (float)pdata[10] , 1.f);
      if (psis[0]== 2052) HF1(4014,  (float)pdata[11] , 1.f); // for afp test 6/10/2008
    }
    
    anode      = (float)pdata[0];
    anode_east = (float)pdata[0];
    anode_west = (float)pdata[4];
    mwpc_cut_var[EAST] = anode_east; //cathode_sum[EAST];
    mwpc_cut_var[WEST] = anode_west; //cathode_sum[WEST];
    mwpc_cut_level[EAST] = 200; //10500;
    mwpc_cut_level[WEST] = 150; //10500;
 
    /* fill "online" QADC sum with MWPC anode cut */

    //dynamically defining the cathode sum cut
    if ( (psis[0]==32)  ) {
	HF1(1077,adcsum_east+adcsum_west,1.f);
    }
    if ( (psis[0]==1||psis[0]==33) && mwpc_cut_var[EAST]>mwpc_cut_level[EAST] ) {
	HF1(9991,adcsum_east,1.f);
    	for(i=0; i<4; i++)
		HF1(ADCCALIB_ID_BASE+900+i,pmtQadc[i],1.f);
    }
    if ( (psis[0]==2||psis[0]==34) && mwpc_cut_var[WEST]>mwpc_cut_level[WEST] ) {
	HF1(9992,adcsum_west,1.f);	
    	for(i=4; i<8; i++)
		HF1(ADCCALIB_ID_BASE+900+i,pmtQadc[i],1.f);
    }
    if ( psis[0]<4 ) HF1(9993,east_twofold,1.f);
    if ( psis[0]<4 ) HF1(9994,west_twofold,1.f);

    /* fill "online" individual PMT QADC histograms: LED events */
      //if ( psis[0]==163 ) {
      if (psis[0]==129 || psis[0]==130 || psis[0]==131 || psis[0]==161 || psis[0]==162 || psis[0]==163) {
      for (i=0; i<8; i++) {
        HF1 (9001+i, pmtQadc[i], 1.f);
      }
    }

    /* fill "online" individual PMT QADC histograms: Bi events */
    if ( psis[0]==32 ) {
      for (i=0; i<8; i++) {
        HF1 (9101+i, pmtQadc[i], 1.f);
      }
    }



    /* book GMS PADC histograms for time after 1 ms */

    if (trigger_time>1000) {
      //if (psis[0]==32) HF1(4014, (float)pdata[2], 1.f);
      if (psis[0]==64) {
	//if (psis[0]==0) {
	HF1(4015, (float)pdata[6], 1.f);
	//HF1(4016, (float)pdata[17], 1.f);
	//HF2(4020, (float)pdata[16], (float)pdata[17], 1.f);

      }
      //if (psis[0]==128) HF1(4016, (float)pdata[10], 1.f);
    }
    bk_close(pevent,pdata+N_ADC);
  }

  pdc3_hf = (nheader<<4)&0xf0 | (nfooter&0xf);

  bk_create(pevent, "EVNB", TID_DWORD, &pdata);
  pdata[0]=qdc_en;
  pdata[1]=tdc_en;
  pdata[2]=pdc1_en;
  pdata[3]=pdc2_en;
  pdata[4]=pdc3_en;
  bk_close(pevent,pdata+5);

  
  bk_create(pevent, "BKHF", TID_DWORD, &bankHeaderFooter);
  bankHeaderFooter[0] = qdc_hf;
  bankHeaderFooter[1] = tdc_hf;
  bankHeaderFooter[2] = pdc1_hf;
  bankHeaderFooter[3] = pdc2_hf;
  bankHeaderFooter[4] = pdc3_hf;
  bk_close(pevent,bankHeaderFooter+5);
  
  bk_create(pevent, "MULT", TID_WORD, &mult);
  multiplicity=0;
  for(i=0;i<16;i++) 
    if(cathE1[i]>150) multiplicity++;
  
  //if(anode_east>anode_cut_e)
    //HF1(60,(float)multiplicity,1.f);
  mult[0]=multiplicity;
  multiplicity=0;
  for(i=0;i<16;i++) 
    if(cathE2[i]>150) multiplicity++;
  //if(anode_east>anode_cut_e)
    //HF1(61,(float)multiplicity,1.f);
  mult[1]=multiplicity;
  bk_close(pevent,mult+2);


  //now do error checking on event number
  if(qdc_hf!=17||tdc_hf!=17||pdc1_hf!=17||pdc2_hf!=17||pdc3_hf!=17)
    header_footer_error_sofar++;  
  if((qdc_en-tdc_en)*(pdc1_en-tdc_en)*(pdc2_en-tdc_en)*(pdc3_en-tdc_en)!=0) 
    nevent_counter_error_sofar++;
  
  /*get east position:*/
  cathE1max=0;cathE2max=0;
  for(i=0;i<16;i++) {
    if(cathE1max<cathE1[i]) {
      cathE1max=cathE1[i];
      i1max=i;
    }
    if(cathE2max<cathE2[i]) {
      cathE2max=cathE2[i];
      i2max=i;
    }
  }

  posEX=-1000.;
  posEY=-1000.;
  if((psis[0]==1||psis[0]==33)&& mwpc_cut_var[EAST]>mwpc_cut_level[EAST] && anode_east<4000. && cathE1max>50.&&cathE2max>50.) {
    tempwx=0;tempw=0;
    for(i=0;i<16;i++) {
      if (cathE1[i]>50.) {
          tempwx = tempwx + cathE1[i]*cathY[i];
          tempw  = tempw  + cathE1[i];
      }
    }
    posEY=tempwx/tempw;
    tempwx=0;tempw=0;
    for(i=0;i<16;i++) {
      if (cathE2[i]>50.) {
          tempwx = tempwx + cathE2[i]*cathX[i];
          tempw  = tempw  + cathE2[i];
      }
    }
    posEX=tempwx/tempw;
    if (posEX>-20.&&posEX<0.&&posEY>20.&&posEY<40.) {
      HF1(9995,adcsum_east,1.f);
    }
    if (posEX>-50.&&posEX<-20.&&posEY>-20.&&posEY<20.) {
      HF1(9996,adcsum_east,1.f);
    }
    if (posEX>-15.&&posEX<5.&&posEY>-20.&&posEY<20.) {
      HF1(9997,adcsum_east,1.f);
    }
    if (posEX>10.&&posEX<40.&&posEY>-20.&&posEY<20.) {
      HF1(9998,adcsum_east,1.f);
    }

    HF1(60,posEX,1.f);
    HF1(61,posEY,1.f);
    HF2(66,posEX,posEY,1.f);
  }


  /*get WEST position:*/
  cathW1max=0;cathW2max=0;
  for(i=0;i<16;i++) {
    if(cathW1max<cathW1[i]) {
      cathW1max=cathW1[i];
      i1max=i;
    }
    if(cathW2max<cathW2[i]) {
      cathW2max=cathW2[i];
      i2max=i;
    }
  }


  posWX=-1000.;
  posWY=-1000.;
  if((psis[0]==2||psis[0]==34)&& mwpc_cut_var[WEST] > mwpc_cut_level[WEST] && anode_west<4000.&&cathW1max>50.&&cathW2max>50.
     ) {
    tempwx=0;tempw=0;
    for(i=0;i<16;i++) {
      if (cathW1[i]>50.) {
          tempwx = tempwx + cathW1[i]*cathY[i];
          tempw  = tempw  + cathW1[i];
	  
	  //printf("cath y %.1f pos %.1f \n", cathW1[i], cathY[i]);
      }
    }
    posWY=tempwx/tempw;

    

    tempwx=0;tempw=0;
    for(i=0;i<16;i++) {
      if (cathW2[i]>50.) {
          tempwx = tempwx + cathW2[i]*cathX[i]*(-1.);
          // -1 corrects for west cathX channels negative of east cathX channels
          tempw  = tempw  + cathW2[i];
	  //printf("cath x %.1f pos %.1f \n", cathW2[i], -cathX[i]);
      }
    }
    posWX=tempwx/tempw;
    //printf("wx = %f\n", posWX);


    if (posWX>-20.&&posWX<0.&&posWY>20.&&posWY<40.) {
      HF1(9995,adcsum_east,1.f);
    }
    if (posWX>-50.&&posWX<-20.&&posWY>-20.&&posWY<20.) {
      HF1(9996,adcsum_east,1.f);
    }
    if (posWX>-15.&&posWX<5.&&posWY>-20.&&posWY<20.) {
      HF1(9997,adcsum_east,1.f);
    }
    if (posWX>10.&&posWX<40.&&posWY>-20.&&posWY<20.) {
      HF1(9998,adcsum_east,1.f);
    }

    HF1(70,posWX,1.f);
    HF1(71,posWY,1.f);
    HF2(76,posWX,posWY,1.f);
  }



  bk_create(pevent, "POSI", TID_FLOAT, &position);
  position[0]=posEX;
  position[1]=posEY;
  position[2]=posWX;
  position[3]=posWY;
  bk_close(pevent,position+4);

  //Do V830, simply skip the header
  n_s830 = bk_locate(pevent, "R_S830", &blockdata);
  if (!(n_s830 == 0 || n_s830 > N_SCLR+1)) { //33 channels, including header
    bk_create(pevent, "S830", TID_DWORD, &pdata);

    //number of scaler channels in this event
    nch = (blockdata[0]>>18)&0x3F;
    //printf("%d\n",(blockdata[0]>>18)&0x3F);
  
    //reset all 32 words that will go into the ntuple
    for(i=0; i<nch; i++) 
      pdata[i]=0;
    
    for(i=1; i<=nch; i++) {//from 1 to 32
      ch = i-1; //skip the header word, read the next 32
      //printf("%d\t%d\n",i,blockdata[i]);
      pdata[ch]=blockdata[i]&0xFFFFFFFF;
    }

    //moved up front
    //hard-coded 1 MHz clock
    clk_of_previous_trig = clk_of_this_trig;
    clk_of_this_trig = pdata[28];
    delta_trigger_time = clk_of_this_trig-clk_of_previous_trig;

    if(psis[0] == 2052) HF1( 4008,(float)clk_of_this_trig/1e6 , 1.f); // absolute time for Adam's switcher detector 


    //periodically update the pedestal values for the MWPC cathodes
    //every 60 seconds. Jianglai 07-30-2008
    //These new peds will be used in the next trigger.
    int this_minute = (clk_of_this_trig/1e6)/60;
    if(this_minute>=1&&clk_of_previous_trig<this_minute*60e6 
       && clk_of_this_trig>this_minute*60e6){
      if(runinfo.online_mode == 1){
	printf("%d Minutes after the run start. %d\t%d. "
	       "Updating MWPC cathode PEDs.\n", this_minute,
	       clk_of_previous_trig, clk_of_this_trig);
      }
      for(i=0;i<16;i++){
	pedE1[i] = get_ped(PDC2CALIB_ID_BASE+i);
	pedE2[i] = get_ped(PDC2CALIB_ID_BASE+i+16);
	pedW1[i] = get_ped(PDCCALIB_ID_BASE+i);
	pedW2[i] = get_ped(PDCCALIB_ID_BASE+i+16);      

	//printf("%f\t%f\t%f\t%f\n",pedE1[i],pedE2[i],pedW1[i],pedW2[i]);

      }
    }
    

    //now do bit operation to figure out if 2-4 trigger 
    //is active in this trigger
    if(psis[0]&0x1==1) trigger_counter[0]++;
    if((psis[0]>>1)&0x1==1) trigger_counter[1]++;    
    
    time_since_last_calc_in_sec = (clk_of_this_trig-clk_start)/1e6;

    const double time2update = 5.25; //s
    
    

    if(time_since_last_calc_in_sec>=time2update){//calculate rates every 5 s
      for(i=0;i<nch;i++){
	counts_since_last_calc[i] = pdata[i]-scaler_saved[i];
	rates_since_last_calc[i] = counts_since_last_calc[i]
	  /time_since_last_calc_in_sec;
	//save all counters values
	scaler_saved[i] = pdata[i];
      }
      //save 1 MHz clock 
      clk_start = clk_of_this_trig;

      //hard-coded: channel 22 is the proton charge counter
      double charge_in_nC = counts_since_last_calc[22];
      double current_in_uA = charge_in_nC/time2update;

      //calculate the trigger counter increment since last calculation
      long trigger_count_since_last[2];
      int ii;
      for(ii=0;ii<2;ii++){
	trigger_count_since_last[ii] 
	  = trigger_counter[ii]-trigger_counter_last[ii];
	//save the value
	trigger_counter_last[ii] = trigger_counter[ii];
      }


      //only do this for online analyzer
      if(runinfo.online_mode == 1){
	//print onto the screen
	printf("Graphite neutron monitor rate in the last %.2f s = %.2f Hz\n",
               time2update, counts_since_last_calc[22]/time2update);
	printf("east 2fold = %.1f Hz\n",rates_since_last_calc[8]);
	printf("west 2fold = %.1f Hz\n",rates_since_last_calc[9]);
	printf("ucn mon1 = %.1f Hz\n",rates_since_last_calc[10]);
	printf("ucn mon2 = %.1f Hz\n",rates_since_last_calc[11]);
	  // factor of 10 to account for pre-scaling
	printf("ucn mon3 = %.1f Hz\n",rates_since_last_calc[12]);
	printf("ucn mon4 = %.1f Hz\n",rates_since_last_calc[13]);
	printf("GMS east ref = %.1f Hz\n", rates_since_last_calc[5]);
	printf("GMS west ref = %.1f Hz\n", rates_since_last_calc[6]);
	printf("GMS global LED = %.1f Hz\n", rates_since_last_calc[7]);
	printf("In Hz, E1 = %.1f, E2 = %.1f, E3 = %.1f, E4 = %.1f\n",
	       rates_since_last_calc[16], rates_since_last_calc[17],
	       rates_since_last_calc[18], rates_since_last_calc[19]);
	printf("In Hz, W1 = %.1f, W2 = %.1f, W3 = %.1f, W4 = %.1f\n",
	       rates_since_last_calc[24], rates_since_last_calc[25],
	       //rates_since_last_calc[26], rates_since_last_calc[27]);
               rates_since_last_calc[30], rates_since_last_calc[27]);

	//deadtime couting for each channel
	printf("E2/4 triggers: %d (SIS), %d (SCLR), ",
	       trigger_count_since_last[0], counts_since_last_calc[8]);
	if(counts_since_last_calc[8]){
	  printf("Diff = %.2f%%\n",
		 (100.0*(counts_since_last_calc[8]-trigger_count_since_last[0]))/counts_since_last_calc[8]);
	} else printf("\n");

	printf("W2/4 triggers: %d (SIS), %d (SCLR), ",
	       trigger_count_since_last[1], counts_since_last_calc[9]);
	if(counts_since_last_calc[9]){
	  printf("Diff = %.2f%%\n",
		 (100.0*(counts_since_last_calc[9]-trigger_count_since_last[1]))/counts_since_last_calc[9]);
	} else printf("\n");
	
	float tdc_alarm_limit 
	  = 0.05*(trigger_count_since_last[0]+trigger_count_since_last[1]);
	//check TDC header/footer and event counter
	printf("Since the last update, there are %d events with "
	       "header/footer problems, and there are %d events "
	       "with event counter issuses! Alarm limit = %.0f\n", 
	       header_footer_error_sofar, nevent_counter_error_sofar,
	       tdc_alarm_limit); 
	
	//ignore tdc. The evnt no readout is strange for the TDC!

	int f_evnt_mis = ((qdc_en!=pdc1_en)||(qdc_en!=pdc2_en)||(qdc_en!=pdc3_en)||(pdc1_en!=pdc2_en)||(pdc1_en!=pdc3_en)||(pdc2_en!=pdc3_en));
	
	/////////// Old misalignment check: LED peak stability
	//float mean = HSTATI(7777, 1, "", 0);
	//float rms = HSTATI(7777, 2, "", 0);
	//printf("East Ref Tube LED trigger mean = %.1f rms = %.1f. ",mean,rms);
	//if((mean<1300||mean>1500) || f_evnt_mis) {
	
	if(f_evnt_mis) {
	  printf("EVENT BUFFER IS PROBABLY MIS-ALIGNED. STOP THIS RUN\n");
	  printf("QDC %d, TDC %d, PDC1 %d, PDC2 %d, PDC3 %d\n", qdc_en, tdc_en, pdc1_en, pdc2_en, pdc3_en);
	  char mycmd[200];
	  sprintf(mycmd,"ssh -f pcucn16 /home/daq/daq_scripts/playalarm");
	  system(mycmd);
	}
	else printf("Event Alignment Okay.\n");

	//HPRINT(7777);
	
	printf("-----------------------------------------------------------\n");

	int eastOK, westOK;
	HNDLE _hDB;
	cm_get_experiment_database(&_hDB, NULL);
	//set the alarm
	if(rates_since_last_calc[0]==0 &&
	   rates_since_last_calc[1]>0 ){
	  eastOK=0;
	  westOK=1;
	  db_set_value(_hDB, 0, "/System/Tmp/ETrigOK",
		       &eastOK, sizeof(eastOK), 1, TID_INT);
	} else if (rates_since_last_calc[0]>0 &&
		   rates_since_last_calc[1]==0 ){
	  westOK=0;
	  eastOK=1;
	  db_set_value(_hDB, 0, "/System/Tmp/WTrigOK",
		       &westOK, sizeof(westOK), 1, TID_INT);
	} else {
	  eastOK=westOK=1;
	  db_set_value(_hDB, 0, "/System/Tmp/ETrigOK",
		       &eastOK, sizeof(eastOK), 1, TID_INT);
	  db_set_value(_hDB, 0, "/System/Tmp/WTrigOK",
		       &westOK, sizeof(westOK), 1, TID_INT);
	}

	
	//open the file that stores the rates
	file_rate = fopen("ucna_rates.dat", "w");
	time(&now);
	sprintf(str, ctime(&now));
	fprintf(file_rate,"Time stamp\t%s\n",str);
	fprintf(file_rate,
		"average current in the last %.2f s = %.2f uA\n",time2update,
		current_in_uA);
        fprintf(file_rate,"east 2fold = %.1f Hz\n",rates_since_last_calc[8]);
        fprintf(file_rate,"west 2fold = %.1f Hz\n",rates_since_last_calc[9]);
        fprintf(file_rate,"ucn mon1 = %.1f Hz\n",rates_since_last_calc[10]);
        fprintf(file_rate,"ucn mon2 = %.1f Hz\n",rates_since_last_calc[11]*10.);
	  // factor of 10 to account for pre-scaling
        fprintf(file_rate,"ucn mon3 = %.1f Hz\n",rates_since_last_calc[12]);
        fprintf(file_rate,"ucn mon4 = %.1f Hz\n",rates_since_last_calc[13]);
        fprintf(file_rate,"GMS east ref = %.1f Hz\n", rates_since_last_calc[5]);
        fprintf(file_rate,"GMS west ref = %.1f Hz\n", rates_since_last_calc[6]);
        fprintf(file_rate,"GMS global LED = %.1f Hz\n", rates_since_last_calc[7]);
        fprintf(file_rate,"In Hz, E1 = %.1f, E2 = %.1f, E3 = %.1f, E4 = %.1f\n",
               rates_since_last_calc[16], rates_since_last_calc[17],
               rates_since_last_calc[18], rates_since_last_calc[19]);
        fprintf(file_rate,"In Hz, W1 = %.1f, W2 = %.1f, W3 = %.1f, W4 = %.1f\n",
               rates_since_last_calc[24], rates_since_last_calc[25],
               //rates_since_last_calc[26], rates_since_last_calc[27]);
               rates_since_last_calc[30], rates_since_last_calc[27]);

        //deadtime couting for each channel
        fprintf(file_rate,"E2/4 triggers: %d (SIS), %d (SCLR), ",
               trigger_count_since_last[0], counts_since_last_calc[8]);
        if(counts_since_last_calc[8]){
          fprintf(file_rate,"Diff = %.2f%%\n",
                 (100.0*(counts_since_last_calc[8]-trigger_count_since_last[0]))/counts_since_last_calc[8]);
        } else fprintf(file_rate,"\n");

        fprintf(file_rate,"W2/4 triggers: %d (SIS), %d (SCLR), ",
               trigger_count_since_last[1], counts_since_last_calc[9]);
        if(counts_since_last_calc[9]){
          fprintf(file_rate,"Diff = %.2f%%\n",
                 (100.0*(counts_since_last_calc[9]-trigger_count_since_last[1]))/counts_since_last_calc[9]);
        } else fprintf(file_rate,"\n");

	//fprintf(file_rate,"East Ref Tube LED trigger mean = %.1f rms = %.1f. ",mean,rms);

	//check TDC header/footer and event counter
	fprintf(file_rate,"Since the last update, there are %d events with "
		"header/footer problems, and there are %d events "
		"with event counter issuses! Alarm limit = %.0f\n", 
		header_footer_error_sofar, nevent_counter_error_sofar,
		tdc_alarm_limit); 
	
	if(header_footer_error_sofar>tdc_alarm_limit
	   ||nevent_counter_error_sofar>tdc_alarm_limit){
	  char mycmd[200];
	  sprintf(mycmd,
		  "ssh -f pcucn16 /home/daq/daq_scripts/playalarm");
	  system(mycmd);	  
	}

	//if(rms>200 || f_evnt_mis) { 
	if(f_evnt_mis) { 
	  fprintf(file_rate,"EVENT BUFFER IS PROBABLY MIS-ALIGNED. STOP THIS RUN\n");
	  fprintf(file_rate, "QDC %d, TDC %d, PDC1 %d, PDC2 %d, PDC3 %d\n", qdc_en, tdc_en, pdc1_en, pdc2_en, pdc3_en);	  
	  char mycmd[200];
	  sprintf(mycmd,
		  "ssh -f pcucn16 /home/daq/daq_scripts/playalarm");
	  system(mycmd);
	} else fprintf(file_rate,"Event Alignment Okay.\n");

        fprintf(file_rate,"-----------------------------------------------------------\n");

	
	fclose(file_rate);
      }//end of test of "online mode"
      HRESET(7777," ");
      HRESET(1077," ");
      header_footer_error_sofar=0;
      nevent_counter_error_sofar=0;
    }

  } else {
    printf("Can not find R_S830 bank from the raw data!!!\n");
  }
  bk_close(pevent,pdata+N_SCLR);


  //below is the code for ntupling the blinded clock
  //clk_of_this_trig is the true clock since run start
  //trigger_time is the true clock since the last h-gx
  bk_create(pevent, "CLK", TID_FLOAT, &clock);
  if(flipper_state==0) {
    clock[0]=((unsigned int)clk_of_this_trig)*blinding_factor1;
    clock[1]=((unsigned int)clk_of_this_trig)*blinding_factor2;
    clock[2]=((unsigned int)trigger_time)*blinding_factor1;
    clock[3]=((unsigned int)trigger_time)*blinding_factor2;
  } else {
    clock[0]=((unsigned int)clk_of_this_trig)*blinding_factor2;
    clock[1]=((unsigned int)clk_of_this_trig)*blinding_factor1;
    clock[2]=((unsigned int)trigger_time)*blinding_factor2;
    clock[3]=((unsigned int)trigger_time)*blinding_factor1;
  }    
  bk_close(pevent,clock+4);

  bk_create(pevent, "FLIP", TID_DWORD, &flipper);
  flipper[0] = upper16;
  bk_close(pevent,flipper+1);

  bk_create(pevent, "DELT", TID_DWORD, &deltaT);
  deltaT[0] = delta_trigger_time;
  bk_close(pevent,deltaT+1);
  
  //fdata[0]=pdata[1];
  //fdata[1]=pdata[2];
  //fdata[2]=pdata[3];
  //fdata[3]=pdata[4];
  //HFN(ADCCALIB_ID_BASE+99,fdata);
  //pmt1=pdata[5];
  //pmt2=pdata[2];
  //pmt3=pdata[3];
  //pmt4=pdata[4];
  //HFNT(ADCCALIB_ID_BASE+99);

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

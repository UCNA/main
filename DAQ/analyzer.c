/********************************************************************\

  Name:         analyzer.c
  Created by:   Stefan Ritt

  Contents:     System part of Analyzer code for sample experiment

  $Log: analyzer.c,v $
  Revision 1.1.1.1  2009/01/20 18:43:29  jliu
  Initial import of the online module.

  Revision 1.10  2002/05/10 05:22:34  pierre
  add MANA_LITE #ifdef

  Revision 1.9  2002/05/09 02:50:28  midas
  Removed initialization of 'Edit on start' by analyzer

  Revision 1.8  2002/05/08 19:54:40  midas
  Added extra parameter to function db_get_value()

  Revision 1.7  2000/11/20 12:29:37  midas
  Added use_tests flag in analyzer request

  Revision 1.6  2000/09/12 12:36:15  midas
  Removed test messages

  Revision 1.5  2000/08/11 11:43:50  midas
  Added cm_msg1 to produce messages which go to a differnt logging file

  Revision 1.4  2000/03/02 22:00:18  midas
  Changed events sent to double

  Revision 1.3  1998/10/29 14:18:19  midas
  Used hDB consistently

  Revision 1.2  1998/10/12 12:18:58  midas
  Added Log tag in header


\********************************************************************/
                                                        
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
#include <stdlib.h>
PAWC_DEFINE(20000000);
#endif


/*-- Globals -------------------------------------------------------*/

/* The analyzer name (client name) as seen by other MIDAS clients   */
char *analyzer_name = "Analyzer";

/* analyzer_loop is called with this interval in ms (0 to disable)  */
INT  analyzer_loop_period = 0; // Every 5 mins 

/* default ODB size */
INT  odb_size = DEFAULT_ODB_SIZE;

/* ODB structures */
RUNINFO          runinfo;
GLOBAL_PARAM     global_param;
EXP_PARAM        exp_param;
TRIGGER_SETTINGS trigger_settings;

/*-- Module declarations -------------------------------------------*/

extern ANA_MODULE scaler_accum_module;
extern ANA_MODULE adc_calib_module;
extern ANA_MODULE adc_summing_module;
extern ANA_MODULE tdc_calib_module;

ANA_MODULE *scaler_module[] = {
  &scaler_accum_module,
  NULL
};

ANA_MODULE *trigger_module[] = {
  &adc_calib_module,
  &adc_summing_module,
  &tdc_calib_module,
  NULL
};

/*-- Bank definitions ----------------------------------------------*/

ASUM_BANK_STR(asum_bank_str);

BANK_LIST trigger_bank_list[] = {

  /* online banks */
  { "R_PD", TID_DWORD, N_ADC+2, NULL },
  { "R_P2", TID_DWORD, N_ADC+2, NULL },
  { "R_P3", TID_DWORD, N_ADC+2, NULL },
  { "R_QD", TID_DWORD, N_ADC+2, NULL },
  { "R_TD", TID_DWORD, N_TDC+2, NULL },
  { "SIS0", TID_DWORD, 2, NULL },
  { "TDEN", TID_DWORD, 1, NULL },
  { "R_S830", TID_DWORD, N_SCLR+1, NULL }, /*33 words, including header*/
  { "R_ST0", TID_DWORD, N_SCLR+1, NULL }, /*33 words, including header*/
  
  
  /* calculated banks */
  { "PADC", TID_DWORD, N_ADC, NULL },
  { "PDC2", TID_DWORD, N_ADC, NULL },
  { "PDC3", TID_DWORD, N_ADC, NULL },
  { "QADC", TID_DWORD, N_ADC, NULL },
  { "TDC0", TID_DWORD, N_TDC, NULL },
  { "S830", TID_DWORD, N_SCLR, NULL },
  { "S820", TID_DWORD, N_SCLR, NULL },
  { "EVNB", TID_DWORD, 5, NULL },
  { "POSI", TID_FLOAT, 4, NULL },
  { "MULT", TID_WORD, 2, NULL },
  { "CADC", TID_FLOAT, N_ADC, NULL },
  { "ASUM", TID_STRUCT, sizeof(ASUM_BANK), asum_bank_str },
  { "CLK", TID_FLOAT, 4, NULL },
  { "FLIP", TID_DWORD, 1, NULL},
  { "DELT", TID_DWORD, 1, NULL},
  { "BKHF", TID_DWORD, 5, NULL },
  { "" },
};

BANK_LIST scaler_bank_list[] = {
  /* online banks */
  { "S820", TID_DWORD,  N_SCLR, NULL },
  { "S830", TID_DWORD,  N_SCLR, NULL },
  { "SIS0", TID_DWORD, 2, NULL },
  /* calculated banks */
  { "ACUM", TID_DOUBLE, N_ADC, NULL },
  { "" },
};

/*-- Event request list --------------------------------------------*/

ANALYZE_REQUEST analyze_request[] = {
  { "Trigger",            /* equipment name */
    1,                    /* event ID */
    TRIGGER_ALL,          /* trigger mask */
    GET_SOME,             /* get some events */
    "SYSTEM",             /* event buffer */
    TRUE,                 /* enabled */
    "", "", 
    NULL,                 /* analyzer routine */
    trigger_module,       /* module list */
    trigger_bank_list,    /* bank list */
    1000,                 /* RWNT buffer size */
    TRUE,                 /* Use tests for this event */
  },

  { "Scaler",             /* equipment name */
    2,                    /* event ID */
    TRIGGER_ALL,          /* trigger mask */
    GET_ALL,              /* get all events */
    "SYSTEM",             /* event buffer */
    TRUE,                 /* enabled */
    "", "", 
    NULL,                 /* analyzer routine */
    scaler_module,        /* module list */
    scaler_bank_list,     /* bank list */
    100,                  /* RWNT buffer size */
  },

  { "" }
};

/*-- Analyzer Init -------------------------------------------------*/
const double get_blinding_factor(const char* barestring);

INT analyzer_init()
{
HNDLE hDB, hKey;
char  str[80];

RUNINFO_STR(runinfo_str);
EXP_PARAM_STR(exp_param_str);
GLOBAL_PARAM_STR(global_param_str);
TRIGGER_SETTINGS_STR(trigger_settings_str);

  /* open ODB structures */
  cm_get_experiment_database(&hDB, NULL);
  db_create_record(hDB, 0, "/Runinfo", strcomb(runinfo_str));
  db_find_key(hDB, 0, "/Runinfo", &hKey);
  if (db_open_record(hDB, hKey, &runinfo, sizeof(runinfo), MODE_READ, NULL, NULL) != DB_SUCCESS)
    {
    cm_msg(MERROR, "analyzer_init", "Cannot open \"/Runinfo\" tree in ODB");
    return 0;
    }

  db_create_record(hDB, 0, "/Experiment/Run Parameters", strcomb(exp_param_str));
  db_find_key(hDB, 0, "/Experiment/Run Parameters", &hKey);
  if (db_open_record(hDB, hKey, &exp_param, sizeof(exp_param), MODE_READ, NULL, NULL) != DB_SUCCESS)
    {
    cm_msg(MERROR, "analyzer_init", "Cannot open \"/Experiment/Run Parameters\" tree in ODB");
    return 0;
    }

  sprintf(str, "/%s/Parameters/Global", analyzer_name);
  db_create_record(hDB, 0, str, strcomb(global_param_str));
  db_find_key(hDB, 0, str, &hKey);
  if (db_open_record(hDB, hKey, &global_param, sizeof(global_param), MODE_READ, NULL, NULL) != DB_SUCCESS)
    {
    cm_msg(MERROR, "analyzer_init", "Cannot open \"%s\" tree in ODB", str);
    return 0;
    }

  db_create_record(hDB, 0, "/Equipment/Trigger/Settings", strcomb(trigger_settings_str));
  db_find_key(hDB, 0, "/Equipment/Trigger/Settings", &hKey);

  if (db_open_record(hDB, hKey, &trigger_settings, sizeof(trigger_settings), MODE_READ, NULL, NULL) != DB_SUCCESS)
    {
    cm_msg(MERROR, "analyzer_init", "Cannot open \"/Equipment/Trigger/Settings\" tree in ODB");
    return 0;
    }

  /* create two blinding factors for the rates */
  const char* string1 = "Seed For UCNA 2011 East Flipper OFF or West Flipper ON";
  const char* string2 = "Seed For UCNA 2011 West Flipper ON or East Flipper OFF";
  blinding_factor1 = get_blinding_factor(string1);
  blinding_factor2 = get_blinding_factor(string2);

  if(blinding_factor1<0.998||blinding_factor1>1.002
     ||blinding_factor2<0.998||blinding_factor2>1.002){
    printf("Something wrong with the blinder!!!! Quit here!!!!\n");
    exit(1);
  }

  return SUCCESS;
}

/*-- Analyzer Exit -------------------------------------------------*/

INT analyzer_exit()
{
  return CM_SUCCESS;
}

/*-- Begin of Run --------------------------------------------------*/

INT ana_begin_of_run(INT run_number, char *error)
{
  return CM_SUCCESS;
}

/*-- End of Run ----------------------------------------------------*/

INT ana_end_of_run(INT run_number, char *error)
{
FILE   *f;
FILE   *fthis;
time_t now;
char   str[256];
int    size;
double n;
HNDLE  hDB;
BOOL   flag;
float  gms[8] = { 0.,0.,0.,0.,0.,0.,0.,0.}; 
float  nevent[2] = {0.,0.};

  cm_get_experiment_database(&hDB, NULL);

  /* update run log if run was written and running online */

  size = sizeof(flag);
  db_get_value(hDB, 0, "/Logger/Write data", &flag, &size, TID_BOOL, TRUE);
  
  size = sizeof(gms);
  db_set_value(hDB, 0, "/Analyzer/Parameters/Global/GMS Sum", &gms, size, 8, TID_FLOAT);
  size = sizeof(nevent);
  db_set_value(hDB, 0, "/Analyzer/Parameters/Global/GMS Events", &nevent, size, 2, TID_FLOAT);



  if (flag && runinfo.online_mode == 1)
    {
    /* update run log */
    size = sizeof(str);
    str[0] = 0;
    db_get_value(hDB, 0, "/Logger/Data Dir", str, &size, TID_STRING, TRUE);
    if (str[0] != 0)
      if (str[strlen(str)-1] != DIR_SEPARATOR)
        strcat(str, DIR_SEPARATOR_STR);
    strcat(str, "runlog.txt");
    
    f = fopen(str, "a");
    fthis = fopen("thisrun.info","w");
    
    time(&now);
    strcpy(str, ctime(&now));
    str[10] = 0;

    fprintf(f, "%s\t%3d\t", str, runinfo.run_number);
    
    fprintf(fthis, "%3d\n", runinfo.run_number); //run number
    fprintf(fthis, "%s\n", str); //start date
    
    strcpy(str, runinfo.start_time); //start time
    str[19] = 0;
    fprintf(f, "%s\t", str+11);
    fprintf(fthis, "%s\n", str+11);

    strcpy(str, ctime(&now));
    str[19] = 0;
    fprintf(f, "%s\t", str+11);
    fprintf(fthis, "%s\n", str+11); //end time

    size = sizeof(n);

    db_get_value(hDB, 0, "/Equipment/Trigger/Statistics/Events sent", &n, &size, TID_DOUBLE, TRUE);

    fprintf(f, "%5.1lfk\t", n/1000);
    fprintf(fthis, "%5.1lfk\n", n/1000);
    fprintf(f, "%s\t", exp_param.comment);
    fprintf(fthis, "%s\n", exp_param.comment);
    fprintf(f, "%i\t", (int)exp_param.gvo);
    fprintf(f, "%i\t", (int)exp_param.sfo);
    fprintf(f, "%i\t", (int)exp_param.calrun);
    fprintf(f, "%s\t", exp_param.flstate);
    fprintf(f, "%5.2fT\t", exp_param.scsf);
    fprintf(f, "%5.2fT\t", exp_param.afp);
    fprintf(f, "%5.2fT\n", exp_param.ppm);
    //  fprintf(fthis, "%s\n", exp_param.comment);
    
    
    fclose(f);
    fclose(fthis);

    system("/home/daq/daq_scripts/db_scripts/fill_run_info_DB.pl");
    }

  return CM_SUCCESS;
}

/*-- Pause Run -----------------------------------------------------*/

INT ana_pause_run(INT run_number, char *error)
{
  return CM_SUCCESS;
}

/*-- Resume Run ----------------------------------------------------*/

INT ana_resume_run(INT run_number, char *error)
{
  return CM_SUCCESS;
}

/*-- Analyzer Loop -------------------------------------------------*/

INT analyzer_loop()
{
  HNDLE  hDB;
  int i,size;
  cm_get_experiment_database(&hDB, NULL);

  // Calculate the average value of the PMT GMS peak during the 5min period.
  if(global_param.nevents[0]>0){
    global_param.ave_gms_0 = global_param.gms_value[0]/global_param.nevents[0];
    global_param.ave_gms_1 = global_param.gms_value[1]/global_param.nevents[0];
    global_param.ave_gms_2 = global_param.gms_value[2]/global_param.nevents[0];
    global_param.ave_gms_3 = global_param.gms_value[3]/global_param.nevents[0];
  } else {
    global_param.ave_gms_0 = 0.;
    global_param.ave_gms_1 = 0.;
    global_param.ave_gms_2 = 0.;
    global_param.ave_gms_3 = 0.;
  }
  if(global_param.nevents[1] >0 ){
    global_param.ave_gms_4 = global_param.gms_value[4]/global_param.nevents[1];
    global_param.ave_gms_5 = global_param.gms_value[5]/global_param.nevents[1];
    global_param.ave_gms_6 = global_param.gms_value[6]/global_param.nevents[1];
    global_param.ave_gms_7 = global_param.gms_value[7]/global_param.nevents[1];
  } else {
    global_param.ave_gms_4 = 0.;
    global_param.ave_gms_5 = 0.;
    global_param.ave_gms_6 = 0.;
    global_param.ave_gms_7 = 0.;
  }
    
  global_param.nevents[0] = 0.;  // Reset the number of events
  global_param.nevents[1] = 0.;

  for(i = 0;i<8;i++){
    global_param.gms_value[i] = 0.;
  }

  db_set_value(hDB, 0, "/Analyzer/Parameters/Global/GMS Sum", &global_param.gms_value, sizeof(global_param.gms_value), 8, TID_FLOAT);
  
  size = (int)sizeof(global_param.ave_gms_0);

  db_set_value(hDB, 0, "/Analyzer/Parameters/Global/GMS Average 0", &global_param.ave_gms_0, size, 1, TID_FLOAT);
  db_set_value(hDB, 0, "/Analyzer/Parameters/Global/GMS Average 1", &global_param.ave_gms_1, size, 1, TID_FLOAT);
  db_set_value(hDB, 0, "/Analyzer/Parameters/Global/GMS Average 2", &global_param.ave_gms_2, size, 1, TID_FLOAT);
  db_set_value(hDB, 0, "/Analyzer/Parameters/Global/GMS Average 3", &global_param.ave_gms_3, size, 1, TID_FLOAT);
  db_set_value(hDB, 0, "/Analyzer/Parameters/Global/GMS Average 4", &global_param.ave_gms_4, size, 1, TID_FLOAT);
  db_set_value(hDB, 0, "/Analyzer/Parameters/Global/GMS Average 5", &global_param.ave_gms_5, size, 1, TID_FLOAT);
  db_set_value(hDB, 0, "/Analyzer/Parameters/Global/GMS Average 6", &global_param.ave_gms_6, size, 1, TID_FLOAT);
  db_set_value(hDB, 0, "/Analyzer/Parameters/Global/GMS Average 7", &global_param.ave_gms_7, size, 1, TID_FLOAT);

  size = (int)sizeof(global_param.nevents);
  db_set_value(hDB, 0, "/Analyzer/Parameters/Global/GMS Events", &global_param.nevents, size, 2, TID_FLOAT);
  
  

  return CM_SUCCESS;

}

/*------------------------------------------------------------------*/
int getseed(const char* barestring)
{
  unsigned long long finalseed;
  int bitcount;
  unsigned int tempout = 0;

  //  This is an attempt to build a portable 64-bit constant
  unsigned long long longmask = (0x7FFFFFFF);
  int i;

  longmask<<=32;
  longmask|=0xFFFFFFFF;


  finalseed = 0;
  bitcount  = 0;
  for (i=0; i<strlen(barestring); i++){
    if ( ((barestring[i])&0xC0)!=0 && ((barestring[i])&0xC0)!=0xC0){
      finalseed = ((finalseed&longmask)<<1) | (((barestring[i])&0x40)>>6);
      bitcount++;
    }
    if ( ((barestring[i])&0x30)!=0 && ((barestring[i])&0x30)!=0x30){
      finalseed = ((finalseed&longmask)<<1) | (((barestring[i])&0x10)>>4);
      bitcount++;
    }
    if ( ((barestring[i])&0xC)!=0 && ((barestring[i])&0xC)!=0xC){
      finalseed = ((finalseed&longmask)<<1) | (((barestring[i])&0x4)>>2);
      bitcount++;
    }
    if ( ((barestring[i])&0x3)!=0 && ((barestring[i])&0x3)!=0x3){
      finalseed = ((finalseed&longmask)<<1) | ((barestring[i])&0x1);
      bitcount++;
    }
  }
  for (i=0; i<(192-bitcount); i++){
    if ((finalseed & 0x800000) == 0x800000){
      finalseed = ((finalseed^0x00000d)<<1) | 0x1;
    } else {
      finalseed<<=1;
    }
  }
  tempout = (finalseed&0xFFFFFFFF) ^ ((finalseed>>32)&0xFFFFFFFF);
  if ((tempout&0x80000000) == 0x80000000){
    tempout = -1 * (tempout&0x7FFFFFFF);
  } else {
    tempout =  1 * (tempout&0x7FFFFFFF);
  }
  return tempout;
}

const double get_blinding_factor(const char* seedstring)
{
  int seed = getseed(seedstring);
  srand(seed);
  double ranno = (double) random() / (double) 0x7fffffff;
  
  const double factor = 0.05;
  const double expectedA = 0.04;
  const double range = factor*expectedA; 

  return 1.0 + 2.0*(ranno-0.5)*range;
}

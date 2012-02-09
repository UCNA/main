/********************************************************************\

  Name:         experim.h
  Created by:   ODBedit program

  Contents:     This file contains C structures for the "Experiment"
                tree in the ODB and the "/Analyzer/Parameters" tree.

                Additionally, it contains the "Settings" subtree for
                all items listed under "/Equipment" as well as their
                event definition.

                It can be used by the frontend and analyzer to work
                with these information.

                All C structures are accompanied with a string represen-
                tation which can be used in the db_create_record function
                to setup an ODB structure which matches the C structure.

  Created on:   Sat May 31 18:19:12 2008

\********************************************************************/

#define EXP_PARAM_DEFINED

typedef struct {
  char      comment[80];
  BOOL      gvo;
  BOOL      sfo;
  BOOL      calrun;
  char      flstate[80];
  float     scsf;
  float     afp;
  float     ppm;
  //  char      octet[80];
} EXP_PARAM;

#define EXP_PARAM_STR(_name) char *_name[] = {\
"[.]",\
"Comment = STRING : [80] west pmt test",\
"Gate Valve Open = BOOL : n",\
"Spin Flipper On = BOOL : n",\
"Calibration Run = BOOL : n",\
"Flapper State = STRING : [80]",\
"SCS Field = FLOAT : 0",\
"AFP Field = FLOAT : 0",\
"PPM Field = FLOAT : 0",\
/*"Octet Type = STRING : [80] A1 Bkgd Off", */\
"",\
NULL }

#define EXP_EDIT_DEFINED

typedef struct {
  char      comment[80];
  BOOL      gvo;
  BOOL      sfo;
  BOOL      calrun;
  char      flstate[80];
  float     scsf;
  float     afp;
  float     ppm;
} EXP_EDIT;

#define EXP_EDIT_STR(_name) char *_name[] = {\
"[.]",\
"Comment = LINK : [35] /Experiment/Run Parameters/Comment",\
"Gate Valve Open = LINK : /Experiment/Run Parameters/Gate Valve Open",\
"Spin Flipper On =  LINK : /Experiment/Run Parameters/Spin Flipper On",\
"Calibration Run =  LINK : /Experiment/Run Parameters/Calibration Run",\
"Flapper State =  LINK : /Experiment/Run Parameters/Flapper State",\
"SCS Field =  LINK : /Experiment/Run Parameters/SCS Field",\
"AFP Field =  LINK : /Experiment/Run Parameters/AFP Field",\
"PPM Field =  LINK : /Experiment/Run Parameters/PPM Field",\
"",\
NULL }

#ifndef EXCL_ADC_CALIBRATION

#define ADC_CALIBRATION_PARAM_DEFINED

typedef struct {
  INT       pedestal[8];
  float     software_gain[8];
  double    histo_threshold;
} ADC_CALIBRATION_PARAM;

#define ADC_CALIBRATION_PARAM_STR(_name) char *_name[] = {\
"[.]",\
"Pedestal = INT[8] :",\
"[0] 174",\
"[1] 194",\
"[2] 176",\
"[3] 182",\
"[4] 185",\
"[5] 215",\
"[6] 202",\
"[7] 202",\
"Software Gain = FLOAT[8] :",\
"[0] 1",\
"[1] 1",\
"[2] 1",\
"[3] 1",\
"[4] 1",\
"[5] 1",\
"[6] 1",\
"[7] 1",\
"Histo threshold = DOUBLE : 20",\
"",\
NULL }

#endif

#ifndef EXCL_ADC_SUMMING

#define ADC_SUMMING_PARAM_DEFINED

typedef struct {
  float     adc_threshold;
} ADC_SUMMING_PARAM;

#define ADC_SUMMING_PARAM_STR(_name) char *_name[] = {\
"[.]",\
"ADC threshold = FLOAT : 5",\
"",\
NULL }

#endif

#ifndef EXCL_TDC_CALIBRATION

#define TDC_CALIBRATION_PARAM_DEFINED

typedef struct {
} TDC_CALIBRATION_PARAM;

#define TDC_CALIBRATION_PARAM_STR(_name) char *_name[] = {\
"",\
NULL }

#endif

#ifndef EXCL_GLOBAL

#define GLOBAL_PARAM_DEFINED

typedef struct {
  float     adc_threshold;
  float     gms_value[8];
  float     nevents[2];
  float     ave_gms_0;
  float     ave_gms_1;
  float     ave_gms_2;
  float     ave_gms_3;
  float     ave_gms_4;
  float     ave_gms_5;
  float     ave_gms_6;
  float     ave_gms_7;
} GLOBAL_PARAM;

#define GLOBAL_PARAM_STR(_name) char *_name[] = {\
"[.]",\
"ADC Threshold = FLOAT : 5",\
"GMS Sum = FLOAT[8] :",\
"[0] = 0.",\
"[1] = 0.",\
"[2] = 0.",\
"[3] = 0.",\
"[4] = 0.",\
"[5] = 0.",\
"[6] = 0.",\
"[7] = 0.",\
"GMS Events = FLOAT[2] :",\
"[0] = 0.",\
"[1] = 0.",\
"GMS Average 0 = FLOAT : 0",\
"GMS Average 1 = FLOAT : 0",\
"GMS Average 2 = FLOAT : 0",\
"GMS Average 3 = FLOAT : 0",\
"GMS Average 4 = FLOAT : 0",\
"GMS Average 5 = FLOAT : 0",\
"GMS Average 6 = FLOAT : 0",\
"GMS Average 7 = FLOAT : 0",\
"",\
NULL }

#endif

#ifndef EXCL_TRIGGER

#define ASUM_BANK_DEFINED

typedef struct {
  float     sum;
} ASUM_BANK;

#define ASUM_BANK_STR(_name) char *_name[] = {\
"[.]",\
"Sum = FLOAT : -1526",\
"",\
NULL }

#define TRIGGER_SETTINGS_DEFINED

typedef struct {
  INT       disc_level;
  char      plu1_logica[8][80];
  char      plu1_logicb[8][80];
  char      plu2_logica[8][80];
  char      plu2_logicb[8][80];
  INT       disc_width;
  INT       gg_delay[8];
  INT       gg_width[8];
  INT       gg_mux;
  INT       longgg_width[2];
  INT       disc_level_array[16];
} TRIGGER_SETTINGS;

#define TRIGGER_SETTINGS_STR(_name) char *_name[] = {\
"[.]",\
"disc_level = INT : 20",\
"plu1_logicA = STRING[8] :",\
"[80] [0] | [1] | [2] | [3] | [4] | [5] | [6] |[7]",\
"[80] [0] & [1]",\
"[80] 0",\
"[80] 0",\
"[80] 0",\
"[80] 0",\
"[80] 0",\
"[80] 0",\
"plu1_logicB = STRING[8] :",\
"[80] [0]|[1]|[2]|[3]",\
"[80] [0]&[1] | [0]&[2] | [0]&[3] | [1]&[2] | [1]&[3] | [2]&[3]",\
"[80] [0]&[1]&[2] | [0]&[1]&[3] | [0]&[2]&[3] | [1]&[2]&[3]",\
"[80] [0]&[1]&[2]&[3]",\
"[80] [0]&[1] | [0]&[2] | [0]&[3] | [1]&[2] | [1]&[3] | [2]&[3]",\
"[80] [0]&[1] | [0]&[2] | [0]&[3] | [1]&[2] | [1]&[3] | [2]&[3]",\
"[80] [0]&[1] | [0]&[2] | [0]&[3] | [1]&[2] | [1]&[3] | [2]",\
"[80] 0",\
"plu2_logicA = STRING[8] :",\
"[80] [0]|[1]|[2]|[3]",\
"[80] [0]&[1] | [0]&[2] | [0]&[3] | [1]&[2] | [1]&[3] | [2]&[3]",\
"[80] [0]&[1]&[2] | [0]&[1]&[3] | [0]&[2]&[3] | [1]&[2]&[3]",\
"[80] [0]&[1]&[2]&[3]",\
"[80] [0]&[1] | [0]&[2] | [0]&[3] | [1]&[2] | [1]&[3] | [2]&[3]",\
"[80] [0]&[1] | [0]&[2] | [0]&[3] | [1]&[2] | [1]&[3] | [2]&[3]",\
"[80] [0]&[1] | [0]&[2] | [0]&[3] | [1]&[2] | [1]&[3] | [2]",\
"[80] 0",\
"plu2_logicB = STRING[8] :",\
"[80] [0]|[1]|[2]|[3]",\
"[80] [0]&[1] | [0]&[2] | [0]&[3] | [1]&[2] | [1]&[3] | [2]&[3]",\
"[80] [0]&[1]&[2] | [0]&[1]&[3] | [0]&[2]&[3] | [1]&[2]&[3]",\
"[80] [0]&[1]&[2]&[3]",\
"[80] [0]&[1] | [0]&[2] | [0]&[3] | [1]&[2] | [1]&[3] | [2]&[3]",\
"[80] [0]&[1] | [0]&[2] | [0]&[3] | [1]&[2] | [1]&[3] | [2]&[3]",\
"[80] [4]",\
"[80] [0]",\
"disc_width = INT : 255",\
"gg_delay = INT[8] :",\
"[0] 0",\
"[1] 0",\
"[2] 0",\
"[3] 0",\
"[4] 0",\
"[5] 0",\
"[6] 0",\
"[7] 0",\
"gg_width = INT[8] :",\
"[0] 60",\
"[1] 60",\
"[2] 60",\
"[3] 0",\
"[4] 0",\
"[5] 0",\
"[6] 0",\
"[7] 0",\
"gg_mux = INT : 1",\
"longgg_width = INT[2] :",\
"[0] 96",\
"[1] 96",\
"disc_level_array = INT[16] :",\
"[0] 20",\
"[1] 20",\
"[2] 20",\
"[3] 20",\
"[4] 20",\
"[5] 20",\
"[6] 20",\
"[7] 20",\
"[8] 20",\
"[9] 20",\
"[10] 20",\
"[11] 20",\
"[12] 20",\
"[13] 20",\
"[14] 20",\
"[15] 20",\
"",\
NULL }

#define TRIGGER_COMMON_DEFINED

typedef struct {
  WORD      event_id;
  WORD      trigger_mask;
  char      buffer[32];
  INT       type;
  INT       source;
  char      format[8];
  BOOL      enabled;
  INT       read_on;
  INT       period;
  double    event_limit;
  DWORD     num_subevents;
  INT       log_history;
  char      frontend_host[32];
  char      frontend_name[32];
  char      frontend_file_name[256];
} TRIGGER_COMMON;

#define TRIGGER_COMMON_STR(_name) char *_name[] = {\
"[.]",\
"Event ID = WORD : 1",\
"Trigger mask = WORD : 0",\
"Buffer = STRING : [32] SYSTEM",\
"Type = INT : 2",\
"Source = INT : 16777215",\
"Format = STRING : [8] MIDAS",\
"Enabled = BOOL : y",\
"Read on = INT : 257",\
"Period = INT : 500",\
"Event limit = DWORD : 0",\
"Num subevents = DWORD : 0",\
"Log history = INT : 0",\
"Frontend host = STRING : [32] pcucn17",\
"Frontend name = STRING : [32] UCNA Frontend",\
"Frontend file name = STRING : [256] frontend.c",\
"",\
NULL }

#endif

#ifndef EXCL_SCALER

#define SCALER_COMMON_DEFINED

typedef struct {
  WORD      event_id;
  WORD      trigger_mask;
  char      buffer[32];
  INT       type;
  INT       source;
  char      format[8];
  BOOL      enabled;
  INT       read_on;
  INT       period;
  double    event_limit;
  DWORD     num_subevents;
  INT       log_history;
  char      frontend_host[32];
  char      frontend_name[32];
  char      frontend_file_name[256];
} SCALER_COMMON;

#define SCALER_COMMON_STR(_name) char *_name[] = {\
"[.]",\
"Event ID = WORD : 2",\
"Trigger mask = WORD : 0",\
"Buffer = STRING : [32] SYSTEM",\
"Type = INT : 1",\
"Source = INT : 0",\
"Format = STRING : [8] MIDAS",\
"Enabled = BOOL : n",\
"Read on = INT : 377",\
"Period = INT : 5000",\
"Event limit = DOUBLE : 0",\
"Num subevents = DWORD : 0",\
"Log history = INT : 0",\
"Frontend host = STRING : [32] pcucn17",\
"Frontend name = STRING : [32] UCNA Frontend",\
"Frontend file name = STRING : [256] frontend.c",\
"",\
NULL }

#endif


#include <stdio.h>
#include <stdlib.h>
#include "midas.h"
#include "mcstd.h"
#include "experim.h"
#include <fcntl.h>
#include "time.h"
#include <sys/time.h>

/* make frontend functions callable from the C framework */
#ifdef __cplusplus
extern "C" {
#endif


/*-- Globals -------------------------------------------------------*/

/* The frontend name (client name) as seen by other MIDAS clients   */
char *frontend_name = "UCNA Frontend";
/* The frontend file name, don't change it */
char *frontend_file_name = __FILE__;

/* frontend_loop is called periodically if this variable is TRUE    */
BOOL frontend_call_loop = TRUE;

/* a frontend status page is displayed with this frequency in ms */
INT display_period = 3000;

/* maximum event size produced by this frontend */
INT max_event_size = 10000;

/* maximum event size for fragmented events (EQ_FRAGMENTED) */
INT max_event_size_frag = 5*1024*1024;

/* buffer size to hold events */
INT event_buffer_size = 10*10000;

INT crate1;

/* Define the base addresses of our modules. */
u_int32_t PADC_base =0x0001;  /* CAEN V785 32 channel PADC */
u_int32_t PDC2_base =0x0003; /* The 2nd PADC, on the right*/
u_int32_t PDC3_base =0x0005; /* The 3rd PADC */
u_int32_t QADC_base =0x0021;  /* CAEN V792 32 channel QADC */
u_int32_t TDC_base  =0x00aa;    /* CAEN V775 32 channel TDC */
u_int32_t SCA820_base=0x0031;  /* SCA: CAEN V820 scaler */
u_int32_t SCA830_base=0x0041;  /* SCA: CAEN V830 scaler */
u_int32_t GG_base=0x0061; /* CAEN V486 gate and delay generator*/
u_int32_t LongGG_base=0x0073; /*CAEN V462 long gate generator */
u_int32_t PLU1_base=0x0051; /* CAEN V495 PLU #1 */
u_int32_t PLU2_base=0x0053; /* CAEN V495 PLU #2 */  
u_int32_t DISC_base=0xeeee; /* CAEN V895 discriminator*/
u_int32_t IO_base=0x383838; /* SIS3600 input register*/


#define N_ADC  32
#define N_TDC  32  
#define N_SCLR 32  /* SCA */

/* define an event counter */
INT evt_ctr;
DWORD start_poll;
DWORD end_read;


WORD dready_adc, dready_tdc,dready_pdc, dready_pdc2, dready_pdc3;
WORD dready_820, dready_830;

TRIGGER_SETTINGS trigger_settings;

/*-- Function declarations -----------------------------------------*/

INT frontend_init();
INT frontend_exit();
INT begin_of_run(INT run_number, char *error);
INT end_of_run(INT run_number, char *error);
INT pause_run(INT run_number, char *error);
INT resume_run(INT run_number, char *error);
INT frontend_loop();

INT read_trigger_event(char *pevent, INT off);
INT read_scaler_event(char *pevent, INT off);

void trigger_update();
u_int16_t eval_string (char* function, u_int16_t replacements[]);
BOOL wait_end_cycle(int transition, BOOL first);

/*-- Equipment list ------------------------------------------------*/

#undef USE_INT

EQUIPMENT equipment[] = {

  { "Trigger",            /* equipment name */
   1, 0,                 /* event ID, trigger mask */
  "SYSTEM",             /* event buffer */
  #ifdef USE_INT
  EQ_INTERRUPT,         /* equipment type */
  #else
  EQ_POLLED,            /* equipment type */
  #endif
  LAM_SOURCE(0,0xFFFFFF),/* event source crate 0, all stations */
  "MIDAS",              /* format */
  TRUE,                 /* enabled */
  RO_RUNNING |          /* read only when running */
  RO_ODB,               /* and update ODB */ 
  500,                  /* poll for 500ms */
  0,                    /* stop run after this event limit */
  0,                    /* number of sub events */
  0,                    /* don't log history */
  "", "", "",
  read_trigger_event,   /* readout routine */
  },

  { "Scaler",             /* equipment name */
    2, 0,                 /* event ID, trigger mask */
    "SYSTEM",             /* event buffer */
        EQ_PERIODIC  
	  /*    EQ_MANUAL_TRIG, */      /* equipment type */
    /*EQ_POLLED*/,            /* SCA: change to polled equipment */
    0,                    /* event source */
    "MIDAS",              /* format */
    FALSE,                 /* enabled */
    RO_RUNNING |
    RO_TRANSITIONS |      /* read when running and on transitions */
    RO_ODB,               /* and update ODB */ 
    5000,                /* read every 5 sec */
    0,                    /* stop run after this event limit */
    0,                    /* number of sub events */
    0,                    /* log history */
    "", "", "",
    read_scaler_event,    /* readout routine */
  },

  { "" }
};

#ifdef __cplusplus
}
#endif

/********************************************************************\
              Callback routines for system transitions

  These routines are called whenever a system transition like start/
  stop of a run occurs. The routines are called on the following
  occations:

  frontend_init:  When the frontend program is started. This routine
                  should initialize the hardware.
  
  frontend_exit:  When the frontend program is shut down. Can be used
                  to releas any locked resources like memory, commu-
                  nications ports etc.

  begin_of_run:   When a new run is started. Clear scalers, open
                  rungates, etc.

  end_of_run:     Called on a request to stop a run. Can send 
                  end-of-run event and close run gates.

  pause_run:      When a run is paused. Should disable trigger events.

  resume_run:     When a run is resumed. Should enable trigger events.

\********************************************************************/

/*-- Frontend Init -------------------------------------------------*/

INT frontend_init()
{
  /* hardware initialization done when frontend starts up*/
  int d1, returncode,i,j/*,crap*/;
  short int threshold,pedestal,bitset2,contr1;
  u_int16_t d3,d2,PLUin[8],PLUout[8];
  DWORD temp;
  HNDLE hDB, hkey;
  TRIGGER_SETTINGS_STR(trigger_settings_str);

  // register for deferred transition
  cm_register_deferred_transition(TR_STOP, wait_end_cycle);

  cm_get_experiment_database(&hDB, NULL);
  db_create_record(hDB,0,"/Equipment/Trigger/Settings",
		   strcomb(trigger_settings_str));
  db_find_key(hDB,0,"Equipment/Trigger/Settings",&hkey);
  if(db_open_record(hDB,hkey,&trigger_settings,sizeof(trigger_settings),
		    MODE_READ,trigger_update,NULL) != DB_SUCCESS) {
    cm_msg(MERROR, "frontend_init", "Cannot open Trigger Settings in ODB");
    return -1;
  }

  /* Open the connection to the VME environment */
  crate1 = open("/tmp/sis1100", O_RDWR, 0 );
  printf("crate1=%x\n",crate1);

  /* Perform a VME system reset */
  returncode = vmesysreset(crate1);
  printf("vmesysreset returncode=%d\n",returncode);
  sleep(10);

  // Note: after a system reset, the SIS 3600 takes a second or two to recover
  // This is currently made up in the routine trigger_update(), which takes
  // several seconds to perform due to repeated shell calls.
  // If trigger_update() is ever optimized, it may be necessary to insert
  // delay, or a polling routine on the SIS 3600 here.  Jeff

  // initialize trigger settings
  trigger_update();


  /*CAEN V486 Gate and Delay generator.*/
  //moved to trigger_update()
 
  /*CAEN V495 PLU*/
  //moved to trigger_update()

  /*CAEN V895 discriminator*/
  //moved to trigger_update()


/*QADC:  initialize the ADC */
/*QADC: software reset*/
  returncode = vme_A32D16_write(crate1, QADC_base << 16 | 0x1006 , 1 << 7);
  printf("vme_A32D16_write returncode=%d\n",returncode);
  vme_A32D16_write(crate1, QADC_base << 16 | 0x1008 , 1 << 7);

/*QADC: set Threshold*/
  threshold=1;  /*threshold times 16 is the threshold in channel*/
  for(i=0;i<32;i++) {
    //if(i==0 || i>4)
    //vme_A32D16_write(crate1 , QADC_base << 16 | 0x1080 + i*2 , threshold*20);
    //else
     vme_A32D16_write(crate1 , QADC_base << 16 | 0x1080 + i*2 , threshold);
     //printf("Setting threshold at %x to %x\n",QADC_base << 16 | 0x1080 + i*2,threshold);
  }

  /*QADC: set Pedestal*/
  pedestal=180;  /*It's recomended to set Pedestal>=60*/
  vme_A32D16_write(crate1, QADC_base << 16 | 0x1060, pedestal);
  /*QADC: disable zero suppression, Bit Set 2 Register*/
  vme_A32D16_read(crate1, QADC_base << 16 | 0x1032, &bitset2);
  vme_A32D16_write(crate1, QADC_base << 16 | 0x1032, bitset2 | 0x10 );
  /*QADC: disable overflow suppression. */
  vme_A32D16_read(crate1, QADC_base << 16 | 0x1032, &bitset2);
  vme_A32D16_write(crate1, QADC_base << 16 | 0x1032, bitset2 | 0x08 );
  /*QADC: BLOCK transfer setting, Control Register 1*/
  vme_A32D16_read(crate1, QADC_base << 16 | 0x1010, &contr1);
  vme_A32D16_write(crate1, QADC_base << 16 | 0x1010, contr1 | 0x04 );
  /*QADC: EMPTY PROG*/
  vme_A32D16_read(crate1, QADC_base << 16 | 0x1032, &bitset2);
  vme_A32D16_write(crate1, QADC_base << 16 | 0x1032, bitset2 | 1<<12 );

  
/*TDC: software reset*/
  vme_A32D16_write(crate1, TDC_base << 16 | 0x1006 , 1 << 7);
  vme_A32D16_write(crate1, TDC_base << 16 | 0x1008 , 1 << 7);

/*TDC: Set the VME TDC to Common Start/Stop Mode */
  //vme_A32D16_write(crate1 , TDC_base << 16 | 0x1034 , 1<<10); /* COM START */
  vme_A32D16_write(crate1 , TDC_base << 16 | 0x1032 , 1<<10); /* COM STOP */

/*TDC: set Threshold*/
  vme_A32D16_write(crate1 , TDC_base << 16 | 0x1080, 0);
  threshold=0x0;  /*threshold times 16 is the threshold in channel*/
  for(i=1;i<32;i++) {
     vme_A32D16_write( crate1 , TDC_base << 16 | 0x1080 + i*2 , threshold);
  }

/*TDC: Set Full Scale Range*/
/*digi=0xff+(225/1060)*(140-t),0xff--140ns, 0x1e--1200ns, 0xc8~180ns */
  //vme_A32D16_write(crate1 , TDC_base << 16 | 0x1060 , 0x1E);
  //vme_A32D16_write(crate1 , TDC_base << 16 | 0x1060 , 0xff);
  vme_A32D16_write(crate1 , TDC_base << 16 | 0x1060 , 0xc8);
  

  /*TDC: disable zero suppression, Bit Set 2 Register*/
  vme_A32D16_read(crate1, TDC_base << 16 | 0x1032, &bitset2);
  vme_A32D16_write(crate1, TDC_base << 16 | 0x1032, bitset2 | 0x10 );
  /*TDC: disable overflow suppression. */
  vme_A32D16_read(crate1, TDC_base << 16 | 0x1032, &bitset2);
  vme_A32D16_write(crate1, TDC_base << 16 | 0x1032, bitset2 | 0x08 );
  /*TDC: BLOCK transfer setting, Control Register 1*/
  vme_A32D16_read(crate1, TDC_base << 16 | 0x1010, &contr1);
  vme_A32D16_write(crate1, TDC_base << 16 | 0x1010, contr1 | 0x04 );
  /*TDC: EMPTY PROG*/
  vme_A32D16_read(crate1, TDC_base << 16 | 0x1032, &bitset2);
  vme_A32D16_write(crate1, TDC_base << 16 | 0x1032, bitset2 | 1<<12 );

  


  vme_A32D16_read(crate1, TDC_base << 16 | 0x1010, &contr1);
  printf("TDC control register 1: %x!!!!!!\n",contr1);
  vme_A32D16_read(crate1, TDC_base << 16 | 0x1032, &bitset2);
  printf("TDC bit set 2 register: %x!!!!!!\n",bitset2);
  
  char message[200];
  sprintf(message, "TDC bit set 2 register: %x!!!!!!", bitset2);
  cm_msg(MLOG, "my program", message);


  /* SCA830: Start with a "software reset": page 32 of manual */
  vme_A32D16_write(crate1 , SCA830_base << 16 | 0x1120 , 0x0);

  /* SCA830: Set scaler control register to preferred acquisition state: */
  /* page 24 of manual; choose ACQ mode=10 periodic; take defaults   */
  /* for all other bits */
  /* set ACQ mode=01 trigger random. */
  ///*set data format=26bit*/

  //Use "trigger random mode"
  //enable the header, set the 6th bit to 1
  vme_A32D16_write(crate1 , SCA830_base << 16 | 0x1108 , 0x21 /*| 1<<2*/);
  
  //block transfer # of events = "1"
  vme_A32D16_write(crate1 , SCA830_base << 16 | 0x1130, 0x1);

  /*SCA830: Channel Enable, P23 of the manual.*/
  /*Enable all channels 16.*/
  vme_A32D32_write(crate1 , SCA830_base << 16 | 0x1100 , 0xffffffff);

  //"V820" is now actually a spare V830, so use 830 code for now
/*   /\* Set up SCA820 *\/ */
/*   /\* SCA820: Start with a "software reset": page 32 of manual *\/ */
/*   vme_A32D16_write(crate1 , SCA820_base << 16 | 0x1120 , 0x0); */

/*   /\* SCA820: Set scaler control register to preferred acquisition state: *\/ */
/*   /\* page 24 of manual; choose ACQ mode=10 periodic; take defaults   *\/ */
/*   /\* for all other bits *\/ */
/*   /\* set ACQ mode=01 trigger random. *\/ */
/*   vme_A32D16_write(crate1 , SCA820_base << 16 | 0x1108 , 0x00 /\*| 1<<2*\/);

 */

  /* SCA820: Start with a "software reset": page 32 of manual */
  vme_A32D16_write(crate1 , SCA820_base << 16 | 0x1120 , 0x0);

  /* SCA830: Set scaler control register to preferred acquisition state: */
  /* page 24 of manual; choose ACQ mode=10 periodic; take defaults   */
  /* for all other bits */
  /* set ACQ mode=01 trigger random. */
  ///*set data format=26bit*/

  //Use "trigger random mode"
  //enable the header, set the 6th bit to 1
  vme_A32D16_write(crate1 , SCA820_base << 16 | 0x1108 , 0x21 /*| 1<<2*/);
  
  //block transfer # of events = "1"
  vme_A32D16_write(crate1 , SCA820_base << 16 | 0x1130, 0x1);

  /*SCA830: Channel Enable, P23 of the manual.*/
  /*Enable all channels 16.*/
  vme_A32D32_write(crate1 , SCA820_base << 16 | 0x1100 , 0xffffffff);  
  

/*PADC: software reset*/
  vme_A32D16_write(crate1 ,  PADC_base << 16 | 0x1006 , 1 << 7);
  vme_A32D16_write(crate1 ,  PADC_base << 16 | 0x1008 , 1 << 7);

/*PADC: set Threshold*/
  threshold=0;//10;  /*threshold times 16 is the threshold in channel*/
  for(i=0;i<32;i++) {
    vme_A32D16_write(crate1 , PADC_base << 16 | 0x1080 + i*2 , threshold);
    printf("Setting threshold at %x to %x\n",PADC_base << 16 | 0x1080 + i*2,threshold);
    //vme_A32D16_read(crate1 , ADC_base << 16 | 0x1080 + i*2 , &d);
    //printf("Read back %x\n",d);
  }

  /*PADC: BLOCK transfer setting, Control Register 1*/
  vme_A32D16_read(crate1, PADC_base << 16 | 0x1010, &contr1);
  vme_A32D16_write(crate1, PADC_base << 16 | 0x1010, contr1 | 0x04 );
  /*PADC: disable zero suppression, Bit Set 2 Register*/
  vme_A32D16_read(crate1, PADC_base << 16 | 0x1032, &bitset2);
  vme_A32D16_write(crate1, PADC_base << 16 | 0x1032, bitset2 | 0x10 );
  /*PADC: disable overflow suppression. */
  vme_A32D16_read(crate1, PADC_base << 16 | 0x1032, &bitset2);
  vme_A32D16_write(crate1, PADC_base << 16 | 0x1032, bitset2 | 0x08 );
  /*PADC: disable sliding scale.*/
  //vme_A32D16_read(crate1, PADC_base << 16 | 0x1032, &bitset2);
  //vme_A32D16_write(crate1, PADC_base << 16 | 0x1032, bitset2 & 0xff7f );
  /*PADC: EMPTY PROG*/
  vme_A32D16_read(crate1, PADC_base << 16 | 0x1032, &bitset2);
  vme_A32D16_write(crate1, PADC_base << 16 | 0x1032, bitset2 | 1<<12 );


/*PDC2: software reset*/
  vme_A32D16_write(crate1 ,  PDC2_base << 16 | 0x1006 , 1 << 7);
  vme_A32D16_write(crate1 ,  PDC2_base << 16 | 0x1008 , 1 << 7);

/*PDC2: set Threshold*/
  threshold=0;//10;  /*threshold times 16 is the threshold in channel*/
  for(i=0;i<32;i++) {
    vme_A32D16_write(crate1 , PDC2_base << 16 | 0x1080 + i*2 , threshold);
    printf("Setting threshold at %x to %x\n",PDC2_base << 16 | 0x1080 + i*2,threshold);
    //vme_A32D16_read(crate1 , ADC_base << 16 | 0x1080 + i*2 , &d);
    //printf("Read back %x\n",d);
  }

  /*PDC2: BLOCK transfer setting, Control Register 1*/
  vme_A32D16_read(crate1, PDC2_base << 16 | 0x1010, &contr1);
  vme_A32D16_write(crate1, PDC2_base << 16 | 0x1010, contr1 | 0x04 );
  /*PDC2: disable zero suppression, Bit Set 2 Register*/
  vme_A32D16_read(crate1, PDC2_base << 16 | 0x1032, &bitset2);
  vme_A32D16_write(crate1, PDC2_base << 16 | 0x1032, bitset2 | 0x10 );
  /*PDC2: disable overflow suppression. */
  vme_A32D16_read(crate1, PDC2_base << 16 | 0x1032, &bitset2);
  vme_A32D16_write(crate1, PDC2_base << 16 | 0x1032, bitset2 | 0x08 );
  /*PDC2: disable sliding scale.*/
  //vme_A32D16_read(crate1, PDC2_base << 16 | 0x1032, &bitset2);
  //vme_A32D16_write(crate1, PDC2_base << 16 | 0x1032, bitset2 & 0xff7f );
  /*PDC2: EMPTY PROG*/
  vme_A32D16_read(crate1, PDC2_base << 16 | 0x1032, &bitset2);
  vme_A32D16_write(crate1, PDC2_base << 16 | 0x1032, bitset2 | 1<<12 );


  /*implementing the new PADC module form TAMU. copy & paste from above*/
  /*PDC3: software reset*/
  vme_A32D16_write(crate1 ,  PDC3_base << 16 | 0x1006 , 1 << 7);
  vme_A32D16_write(crate1 ,  PDC3_base << 16 | 0x1008 , 1 << 7);

  /*PDC3: set Threshold*/
  threshold=0;//10;  /*threshold times 16 is the threshold in channel*/
  for(i=0;i<32;i++) {
    vme_A32D16_write(crate1 , PDC3_base << 16 | 0x1080 + i*2 , threshold);
    printf("Setting threshold at %x to %x\n",PDC3_base << 16 | 0x1080 + i*2,threshold);
    //vme_A32D16_read(crate1 , ADC_base << 16 | 0x1080 + i*2 , &d);
    //printf("Read back %x\n",d);
  }

  /*PDC3: BLOCK transfer setting, Control Register 1*/
  vme_A32D16_read(crate1, PDC3_base << 16 | 0x1010, &contr1);
  vme_A32D16_write(crate1, PDC3_base << 16 | 0x1010, contr1 | 0x04 );
  /*PDC3: disable zero suppression, Bit Set 2 Register*/
  vme_A32D16_read(crate1, PDC3_base << 16 | 0x1032, &bitset2);
  vme_A32D16_write(crate1, PDC3_base << 16 | 0x1032, bitset2 | 0x10 );
  /*PDC3: disable overflow suppression. */
  vme_A32D16_read(crate1, PDC3_base << 16 | 0x1032, &bitset2);
  vme_A32D16_write(crate1, PDC3_base << 16 | 0x1032, bitset2 | 0x08 );
  /*PDC3: disable sliding scale.*/
  //vme_A32D16_read(crate1, PDC3_base << 16 | 0x1032, &bitset2);
  //vme_A32D16_write(crate1, PDC3_base << 16 | 0x1032, bitset2 & 0xff7f );
  /*PDC3: EMPTY PROG*/
  vme_A32D16_read(crate1, PDC3_base << 16 | 0x1032, &bitset2);
  vme_A32D16_write(crate1, PDC3_base << 16 | 0x1032, bitset2 | 1<<12 );


  /* SCA: Set the time between internal triggers; see page 24 of the */
  /* manual.  Some useful values:  */
  /*    1 sec = 0x2625a0 (I set it to this, for now) */
  /*    10 sec = 0x17d7840 */
  /*    0.1 sec = 0x3d090 */
  //  vme_A24D32_write(crate1 , SCA_base << 16 | 0x1104 , 0x2625a0);
  //vme_A24D32_write(crate1 , SCA_base << 16 | 0x1104 , 0x2625a0);
  //vme_A24D32_write(crate1 , SCA_base << 16 | 0x1104 , 0x3d090);


  /*SIS3600 Input register*/
  /*global reset*/
  vme_A32D16_write(crate1, IO_base << 8 | 0x60, 1);
  
  /*configure:enable external next, strobe mode*/  
  vme_A32D32_write(crate1, IO_base << 8 | 0x0, 0xfe01ff00);
  //vme_A32D32_write(crate1, IO_base <<8 | 0x0, 0x10000);

  /*configure:enable external next, coincidence mode mode*/  
  //vme_A32D32_write(crate1, IO_base << 8 | 0x0, 0xf609ff00);
  
  temp=0;
  vme_A32D32_read(crate1, IO_base << 8 | 0x0, &temp);
  printf("temp=%x!!!!!!!!!!\n",temp);

  return SUCCESS;
}

void trigger_update() {
  
  int i,j,ii,jj;
  int er;
  u_int16_t d2,d3,PLUin[8],PLUout[8];

  /*CAEN V486 Gate and Delay generator.*/
  er = vme_A32D16_write(crate1, GG_base << 16 | 0x0040, 0x20); /*normal mode*/
  printf("return code vme_A32D16_write = %d\n",er);
  //vme_A32D16_write(crate1, GG_base << 16 | 0x0040, 0xe0|1<<3); /*test mode*/ 
  for (i=0;i<8;i++) {
    er = vme_A32D16_write(crate1, GG_base << 16 | (0x0010 + i*2),
		     trigger_settings.gg_delay[i]);
    printf("return code vme_A32D16_write = %d\n",er);
    er = vme_A32D16_write(crate1, GG_base << 16 | (0x0020 + i*2),
		     trigger_settings.gg_width[i]);
    printf("return code vme_A32D16_write = %d\n",er);
  }
  //select mux output
  vme_A32D16_write(crate1, GG_base << 16 | 0x0040,
		   trigger_settings.gg_mux << 5);

  /* CAEN V462 Gate and Delay Generator */
  /* Note: 24-bit addressing only! */
  // set gate length on channel 0
  er = vme_A24D16_write(crate1, LongGG_base << 16 | 0x0002,
			trigger_settings.longgg_width[0]>>16); //MSB
  
  vme_A24D16_write(crate1, LongGG_base << 16 | 0x0004,
		   trigger_settings.longgg_width[0]&0xffff); //LSB
  // set gate length on channel 1
  vme_A24D16_write(crate1, LongGG_base << 16 | 0x0006,
		   trigger_settings.longgg_width[1]>>16); //MSB
  vme_A24D16_write(crate1, LongGG_base << 16 | 0x0008,
		   trigger_settings.longgg_width[1]&0xffff); //LSB



  /*CAEN V895 discriminator*/
  /*enable*/
  /* jianglai implemented array variable so that different channel */
  /* can have different threshold */
  vme_A32D16_write(crate1, DISC_base << 16 | 0x4A, 0xFFFF);
  for(i=0; i<16; i++) {
    vme_A32D16_write(crate1, DISC_base << 16 | (0x00+i*2) ,
		     trigger_settings.disc_level_array[i]);
    printf("disc level is set to %d.\n",trigger_settings.disc_level_array[i]);
  }

  /*output width*/
  vme_A32D16_write(crate1, DISC_base << 16 | 0x40,
		   trigger_settings.disc_width);
  vme_A32D16_write(crate1, DISC_base << 16 | 0x42,
		   trigger_settings.disc_width);
  printf("disc width is set to %d.\n",trigger_settings.disc_width);
  /*Majority threshold*/
  //vme_A32D16_write(crate1, DISC_base << 16 | 0x48, 19); /*level 2*/

#define PLU1
  //#undef PLU1            //Master trigger
#ifdef PLU1
  //PLU #1
  /*CAEN V495 PLU*/
  /*prepare to write RAM*/
  vme_A32D16_write(crate1,PLU1_base << 16 | 0x002,0xf70f); 
  /*default value of the configuration register is f800.*/
  for(i=0; i<256; i++) {
    /*RAM of section A*/
    for(j=0; j<8; j++)
      PLUin[j]=(i>>j) & 0x01;
    for(jj=0; jj<8; jj++)
      PLUout[jj]=eval_string(trigger_settings.plu1_logica[jj],PLUin);
    d2=0;
    for(j=0; j<8; j++) 
      d2 = d2 | (PLUout[j]<<j);
    d2=d2 & 0x0ff;
    vme_A32D16_write(crate1, PLU1_base << 16 | (0x200+i*2), d2);
    vme_A32D16_read(crate1,PLU1_base << 16 | (0x200+i*2),&d3);
    //    printf("Setting RAM at %x to %x\n",PLU1_base << 16 | (0x200+i*2),d3);
  }
  printf("PLU1 section A logic is reset!\n");
  for(i=0; i<256; i++) {
    /*RAM of section B*/
    for(j=0; j<8; j++)
      PLUin[j]=(i>>j) & 0x01;
    for(jj=0; jj<8; jj++)
      PLUout[jj]=eval_string(trigger_settings.plu1_logicb[jj],PLUin);
    d2=0;
    for(j=0; j<8; j++) 
      d2 = d2 | (PLUout[j]<<j);
    d2=d2 & 0x0ff;
    vme_A32D16_write(crate1, PLU1_base << 16 | (0x400+i*2), d2);
  }
  printf("PLU1 section B logic is reset!\n");
  
  /*setting PLU to LOOK-UP TABLE mode*/
  vme_A32D16_write(crate1,PLU1_base << 16 | 0x002,0xf40c); 
  //end of PLU #1
#endif

#define PLU2
#ifdef PLU2          //East(0-7), West(8-15) trigger
  //PLU #2
  /*CAEN V495 PLU*/
  /*prepare to write RAM*/
  vme_A32D16_write(crate1,PLU2_base << 16 | 0x002,0xf70f); 
  /*default value of the configuration register is f800.*/
  for(i=0; i<256; i++) {
    /*RAM of section A*/
    for(j=0; j<8; j++)
      PLUin[j]=(i>>j) & 0x01;
    for(jj=0; jj<8; jj++)
      PLUout[jj]=eval_string(trigger_settings.plu2_logica[jj],PLUin);
    d2=0;
    for(j=0; j<8; j++) 
      d2 = d2 | (PLUout[j]<<j);
    d2=d2 & 0x0ff;
    vme_A32D16_write(crate1, PLU2_base << 16 | (0x200+i*2), d2);
  }
  printf("PLU2 section A logic is reset!\n");
  for(i=0; i<256; i++) {
    /*RAM of section B*/
    for(j=0; j<8; j++)
      PLUin[j]=(i>>j) & 0x01;
    for(jj=0; jj<8; jj++)
      PLUout[jj]=eval_string(trigger_settings.plu2_logicb[jj],PLUin);
    d2=0;
    for(j=0; j<8; j++) 
      d2 = d2 | (PLUout[j]<<j);
    d2=d2 & 0x0ff;
    vme_A32D16_write(crate1, PLU2_base << 16 | (0x400+i*2), d2);
  }
  printf("PLU2 section B logic is reset!\n");
  
  /*setting PLU to LOOK-UP TABLE mode*/
  vme_A32D16_write(crate1,PLU2_base << 16 | 0x002,0xf40c); 
  //end of PLU #2
#endif
}


u_int16_t eval_string (char* function, u_int16_t replacements[]) {
  // returns numerical evaluation of function, replacing "[0]" in string
  // with value of replacements[0], etc.
  char cmd[80];
  char* ptr;
  char temp_char[3];
  int ii=0;
  int jj;
  FILE *fp;

  int value;

  strcpy(cmd,"echo $((");
  ptr=cmd+strlen(cmd);

  // substitute things of form [#] to val of replacements in cmd string
  while (*function) {
    if(*function!='[') {
      *ptr=*function;
    } else {
      function++;
      sprintf(ptr,"%x",replacements[atoi(function)]);
      function++;
    }
    function++;
    ptr++;
  }
  *ptr='\0';

  // finish off and evaluate cmd using bash
  strcat(cmd,"))");

  if((fp = popen(cmd,"r")) != NULL) { 
    fscanf(fp,"%d",&value);
    if(value != 0 && value != 1) {
      printf("err in PLU settings!!! %s %s\n",cmd);
      value = 0;
    }
  } else printf("err in popen!!!\n");
  pclose(fp);

  //  printf("eval_string %s %d\n",cmd,value);

  return value;

}

/*-- Frontend Exit -------------------------------------------------*/

INT frontend_exit()
{

/* Close the connection to the VME environment */
  close(crate1);

  return SUCCESS;
}

/*-- Begin of Run --------------------------------------------------*/

INT begin_of_run(INT run_number, char *error)
{
  int d1,a;
  DWORD d32;  


/* reset event counter */
  evt_ctr = 0;

  /*QADC: DATA reset*/
  vme_A32D16_write(crate1 ,  QADC_base << 16 | 0x1032 , 1 << 2);
  vme_A32D16_write(crate1 ,  QADC_base << 16 | 0x1034 , 1 << 2);
  /*TDC: DATA reset*/
  vme_A32D16_write(crate1 ,  TDC_base << 16 | 0x1032 , 1 << 2);
  vme_A32D16_write(crate1 ,  TDC_base << 16 | 0x1034 , 1 << 2);
  /*PDC: DATA reset*/
  vme_A32D16_write(crate1 ,  PADC_base << 16 | 0x1032 , 1 << 2);
  vme_A32D16_write(crate1 ,  PADC_base << 16 | 0x1034 , 1 << 2);
  /*PDC2: DATA reset*/
  vme_A32D16_write(crate1 ,  PDC2_base << 16 | 0x1032 , 1 << 2);
  vme_A32D16_write(crate1 ,  PDC2_base << 16 | 0x1034 , 1 << 2);
  /*PDC3: DATA reset*/
  vme_A32D16_write(crate1 ,  PDC3_base << 16 | 0x1032 , 1 << 2);
  vme_A32D16_write(crate1 ,  PDC3_base << 16 | 0x1034 , 1 << 2);
  



  /*QADC & TDC &PADC1,2: Event counter reset*/
  vme_A32D16_write(crate1 ,  QADC_base << 16 | 0x1040 , 1);
  vme_A32D16_write(crate1 ,  TDC_base << 16 | 0x1040 , 1);
  vme_A32D16_write(crate1 ,  PADC_base << 16 | 0x1040 , 1);
  vme_A32D16_write(crate1 ,  PDC2_base << 16 | 0x1040 , 1);
  vme_A32D16_write(crate1 ,  PDC3_base << 16 | 0x1040 , 1);

  /* SCA830: software clear of the counter registers, page 33 of the manual*/
  vme_A32D16_write(crate1 , SCA830_base << 16 | 0x1122 , 0x0);

  //note "820" is physically a "830" now
  /* SCA820: software clear of the counter registers, page 33 of the manual*/
  vme_A32D16_write(crate1 , SCA820_base << 16 | 0x1122 , 0x0);

  /* SIS3600 input register, clear FIFO and logic*/
  vme_A32D16_write(crate1, IO_base << 8 | 0x20, 1);
  /*enable NEXT clock logic*/
  vme_A32D16_write(crate1, IO_base<<8 | 0x28, 1);
  /*issue first next clock pulse to start counting.*/
  //vme_A32D16_write(crate1, IO_base<<8 | 0x24, 1);

 /* 
  vme_A32D32_read(crate1, IO_base<<8 |0x0, &d32);
  if(d32 & 0x100) {
    ss_printf(0, 11, "SIS3600 FIFO empty!");
  }
  else {
    vme_A32D32_read(crate1, IO_base << 8 | 0x100, &d32);
    ss_printf(0, 11, "read out SIS3600!  ");
  }
 */

  ss_printf(0,8, "================================================================================");
  ss_printf(0,9, "V820 readings (0-3)  :");
  ss_printf(0,10,"V820 readings (4-7)  :");
  ss_printf(0,11,"V820 readings (8-11) :");
  ss_printf(0,12,"V820 readings (12-15):");
  ss_printf(0,13, "================================================================================");
  
  ss_printf(0,15,"end of run warnings:");
  ss_printf(0,16,"                                                             ");
  ss_printf(0,17,"                                                             ");
  ss_printf(0,18,"================================================================================");
  ss_printf(0,19,"Warnings during the run:");
  ss_printf(0,20,"                                                             ");
  ss_printf(0,19,"                                                             ");
    

/*enable main trigger*/
  s3100_control_write(crate1, 0x80, 0x01);
  return SUCCESS;
}


/*-- Defer End of Run ----------------------------------------------*/
BOOL wait_end_cycle(int transition, BOOL first) {
  static DWORD tfirst;
  if (first) {
    tfirst = ss_time();
    // disable main trigger
    s3100_control_write(crate1, 0x80, 0x010000); 
    ss_printf(40,18,"Stop requested");
  }
  // wait 3 seconds
  return ss_time() > tfirst + 3;
}


/*-- End of Run ----------------------------------------------------*/

INT end_of_run(INT run_number, char *error) {
  INT i=0;
  DWORD d32,temp;
  u_int16_t d1;
  
  /*disable main trigger*/
  s3100_control_write(crate1, 0x80, 0x010000); 
		  
  /*disable next logic for SIS3600*/
  vme_A32D16_write(crate1, IO_base<<8 | 0x2c, 1);

  /* check if SIS3600 has any events left over*/
  vme_A32D32_read(crate1, IO_base<<8 |0x0, &d32);
  while(!(d32 & 0x100)) {
    i++;
    vme_A32D32_read(crate1, IO_base << 8 | 0x100, &temp);
    vme_A32D32_read(crate1, IO_base<<8 |0x0, &d32);
  }
  if(i)
    ss_printf(0,13,"SIS3600 has %d events remaining at end of run!  ", i);
  else
    ss_printf(0,13,"SIS3600 has no events remaining at end of run!  ");

  /* check of TDC has any events still in buffer*/
  vme_A32D16_read(crate1 ,  TDC_base << 16 | 0x100E  ,&d1);
  if ((d1) & 0x01 ) /*DRADY is on bit 1*/
    ss_printf(0,14,"TDC still has events remaining at end of run!", i);
  else
    ss_printf(0,14,"TDC has no events remaining at end of run!   ", i);

  /* automatically launch replay command */
  char mycmd[200];
  sprintf(mycmd,
	  "ssh -f pcucn18 /home/daq/daq_scripts/start_ana %05d",run_number);
  system(mycmd);
  return SUCCESS;
}

/*-- Pause Run -----------------------------------------------------*/

INT pause_run(INT run_number, char *error)
{
  return SUCCESS;
}

/*-- Resume Run ----------------------------------------------------*/

INT resume_run(INT run_number, char *error)
{
  return SUCCESS;
}

/*-- Frontend Loop -------------------------------------------------*/

INT frontend_loop()
{
  /* if frontend_call_loop is true, this routine gets called when
     the frontend is idle or once between every event */
  if(end_read - start_poll !=0){
    ss_printf(0,20,"Polling time    : %15.10f \n",(float)(end_read - start_poll));///CLOCKS_PER_SEC);
    ss_printf(0,21,"Start Poll Time : %15.10f \n",(float) start_poll);///CLOCKS_PER_SEC);
    ss_printf(0,22,"End Read Timing : %15.10f \n",(float) end_read);//CLOCKS_PER_SEC);
    ss_printf(0,23,"Clocks PER SEC  : %15.10f \n",(float) CLOCKS_PER_SEC); 
  }
  return SUCCESS;
}

/*------------------------------------------------------------------*/

/********************************************************************\
  
  Readout routines for different events

\********************************************************************/

/*-- Trigger event routines ----------------------------------------*/

INT poll_event(INT source, INT count, BOOL test) {
/* Polling routine for events. Returns TRUE if event
   is available. If test equals TRUE, don't return. The test
   flag is used to time the polling */
  int   i,p;
  u_int16_t d1,d2,d3,d4,d5,d6,d7;

/*   static DWORD my_previous_poll=0; */
/*   DWORD this_time = ss_millitime(); */
/*   printf("time diff %d\n",this_time-my_previous_poll); */
/*   my_previous_poll = this_time; */

  for (i=0 ; i<count ; i++) {
    /* Polling is done on DREADY of the ADC Status Register 1 */
    /* Then we will read out both ADC and TDC; error checking */
    /* is done in the event readout section to make sure we stay */
    /* in sync. */
    vme_A32D16_read(crate1 ,  QADC_base << 16 | 0x100E  ,&d1);
    vme_A32D16_read(crate1 ,  TDC_base  << 16 | 0x100E  ,&d2);
    vme_A32D16_read(crate1 ,  PADC_base << 16 | 0x100E  ,&d3);
    vme_A32D16_read(crate1 ,  PDC2_base << 16 | 0x100E  ,&d4);
    vme_A32D16_read(crate1 ,  PDC3_base << 16 | 0x100E  ,&d5);
    vme_A32D16_read(crate1 ,  SCA820_base << 16 | 0x110E  ,&d6);
    vme_A32D16_read(crate1 ,  SCA830_base << 16 | 0x110E  ,&d7);
   
    if ((d1 & d2 & d3 & d4 & d5 & d6 & d7) & 0x01 ) /*DREADY is on bit 1*/
      if (!test) {
	dready_adc=d1 & 0x0001;
	dready_tdc=d2 & 0x0001;
	dready_pdc=d3 & 0x0001;
	dready_pdc2=d4 & 0x0001;
	dready_pdc3=d5 & 0x0001;
	dready_820=d6 & 0x0001;
	dready_830=d7 & 0x0001;
	/*software trigger for Scaler830*/
	//vme_A32D16_write(crate1, SCA830_base << 16 | 0x1124, 0x01);
	return 1;
      }
  }
  
  return 0;
}

/*-- Interrupt configuration ---------------------------------------*/

INT interrupt_configure(INT cmd, INT source, PTYPE adr)
{
  switch(cmd)
    {
    case CMD_INTERRUPT_ENABLE:
      break;
    case CMD_INTERRUPT_DISABLE:
      break;
    case CMD_INTERRUPT_ATTACH:
      break;
    case CMD_INTERRUPT_DETACH:
      break;
    }
  return SUCCESS;
}

/*-- Event readout -------------------------------------------------*/

#define USE_BLT_MODE 1   /*whether to use Block Transfer*/

INT read_trigger_event(char *pevent, INT off)
{
  WORD i,ch,d1,valid;   /*in MIDAS: WORD is unsigned short int.*/
  DWORD d32,buf_addr,*pdata,*blockdata,*eventn,temp;
  u_int32_t nch, nchread;
  INT  q, returncode;
  u_int32_t sis_status_bit8, sis_status_bit9, sis_status_bit10, sis_status_bit11, sis_status_bit12;
  //char bkname[20];
  static INT err_number = 1;
 
  //start_poll = ss_millitime();

  struct timeval tv;
  struct timezone tz;
  struct tm *tm;
  gettimeofday(&tv, &tz);
  tm=localtime(&tv.tv_sec);
  start_poll = tv.tv_usec;
  
 /* init bank structure */
  bk_init(pevent);

  /* create SIS0 bank */
  bk_create(pevent, "SIS0", TID_DWORD, &pdata);
  vme_A32D32_read(crate1, IO_base<<8 | 0x0, &d32);
  /*
  i=0;
  while(!(d32 & 0x100)) {
    i++;
    vme_A32D32_read(crate1, IO_base << 8 | 0x100, &temp);
    vme_A32D32_read(crate1, IO_base<<8 |0x0, &d32);
  }
  *pdata++=temp;
  *pdata++=i;
  */
  


  sis_status_bit8=(d32>>8)&0x1;
  sis_status_bit9=(d32>>9)&0x1;
  sis_status_bit10=(d32>>10)&0x1;
  sis_status_bit11=(d32>>11)&0x1;
  sis_status_bit12=(d32>>12)&0x1;

  if(sis_status_bit8||sis_status_bit11||sis_status_bit12){
    printf("SIS FIFO empty %x almost empty %x half full %x almost full %x full %x \n", sis_status_bit8, sis_status_bit9, sis_status_bit10, sis_status_bit11, sis_status_bit12);
  }
  if(sis_status_bit8) {
    ss_printf(0,17,"SIS3600 FIFO empty!\n");
    *pdata++=-1;
  }
  else {
    vme_A32D32_read(crate1, IO_base << 8 | 0x100, &d32);
    *pdata++ = d32;
  }
  
  bk_close(pevent, pdata);

#ifdef USE_BLT_MODE
  if ((dready_tdc == 0) || (dready_adc == 0) || (dready_pdc == 0) || (dready_pdc2 == 0) || (dready_pdc3 == 0)) {
    err_number++;
    ss_printf(0,24,"%d Inconsitent ready conditions: tdc %x adc %x pdc %x pdc2 %x pdc3 %x \n",err_number,dready_tdc,dready_adc,
           dready_pdc, dready_pdc2,dready_pdc3);
  }

  if(dready_tdc) {
    /* create raw TDC banks */
    bk_create(pevent, "R_TD", TID_DWORD, &blockdata);
    buf_addr=TDC_base << 16 | 0x0000;
    returncode=vme_A32BLT32_read(crate1, buf_addr, blockdata, 34, &nchread);
    if(nchread !=34 || returncode <0) {
      ss_printf(0,18,"TDC block transfer error! %d   %d\n",  nchread, returncode );
      // return 1;
    }
    bk_close(pevent,blockdata+34);  

    /*
    bk_create(pevent, "TDEN", TID_DWORD, &eventn);
    vme_A32D16_read(crate1, TDC_base<<16 | 0x1024, &d32);
    *eventn++=d32 & 0x0ffff;
    bk_close(pevent,eventn);
    */
    dready_tdc=0;
  }

  if(dready_adc) {
    /* create raw Charge ADC banks */
    bk_create(pevent, "R_QD", TID_DWORD, &blockdata);
    buf_addr=QADC_base << 16 | 0x0000;
    returncode=vme_A32BLT32_read(crate1, buf_addr, blockdata, 34, &nchread);
    if(nchread !=34 || returncode <0) {
      ss_printf(0,19,"QADC block transfer error! %d   %d\n",  nchread, returncode );
      //return 1;
    }
    bk_close(pevent,blockdata+34);  
    dready_adc=0;
  }

  if(dready_pdc) {
    /* create raw peak ADC banks */
    bk_create(pevent, "R_PD", TID_DWORD, &blockdata);
    buf_addr=PADC_base << 16 | 0x0000;
    returncode=vme_A32BLT32_read(crate1, buf_addr, blockdata, 34, &nchread);
    if(nchread !=34 || returncode <0) {
      ss_printf(0,20,"PADC1 block transfer error! %d   %d\n",  nchread, returncode );
      //return 1;
    }
    bk_close(pevent,blockdata+34);  
    dready_pdc=0;
  }

  if(dready_pdc2) {
    /* create raw peak ADC banks */
    bk_create(pevent, "R_P2", TID_DWORD, &blockdata);
    buf_addr=PDC2_base << 16 | 0x0000;
    returncode=vme_A32BLT32_read(crate1, buf_addr, blockdata, 34, &nchread);
    if(nchread !=34 || returncode <0) {
      ss_printf(0,20,"PADC2 block transfer error! %d   %d\n",  nchread, returncode );
      //return 1;
    }
    bk_close(pevent,blockdata+34);  
    dready_pdc2=0;
  }


  if(dready_pdc3) {
    /* create raw peak ADC banks */
    bk_create(pevent, "R_P3", TID_DWORD, &blockdata);
    buf_addr=PDC3_base << 16 | 0x0000;
    returncode=vme_A32BLT32_read(crate1, buf_addr, blockdata, 34, &nchread);
    if(nchread !=34 || returncode <0) {
      ss_printf(0,20,"PADC3 block transfer error! %d   %d\n",  nchread, returncode );
      //return 1;
    }
    bk_close(pevent,blockdata+34);  
    dready_pdc3=0;
  }

  if(dready_830) {
    /* create RAW SCLR830 bank */
    bk_create(pevent, "R_S830", TID_DWORD, &blockdata);
    buf_addr=SCA830_base << 16 | 0x0000;
    returncode=vme_A32BLT32_read(crate1, buf_addr, blockdata, N_SCLR+1, &nchread);
    if(nchread != N_SCLR+1|| returncode <0) {
      ss_printf(0,20,"SCA830 block transfer error! %d   %d\n",  nchread, returncode );
    }
    bk_close(pevent, blockdata+N_SCLR+1);
    dready_830=0;
  }
  
  if(dready_820) {
    /* create RAW SCLR820 bank */
    bk_create(pevent, "R_ST0", TID_DWORD, &blockdata);
    buf_addr=SCA820_base << 16 | 0x0000;
    returncode=vme_A32BLT32_read(crate1, buf_addr, blockdata, N_SCLR+1, &nchread);
    if(nchread != N_SCLR+1|| returncode <0) {
      ss_printf(0,20,"SCA820 block transfer error! %d   %d\n",  nchread, returncode );
    }
    bk_close(pevent, blockdata+N_SCLR+1);
    dready_820=0;
  }

  /*
  bk_create(pevent, "QADC", TID_DWORD, &pdata);

  for(i=0; i<N_ADC; i++) 
    pdata[i]=-1;

  nch=(blockdata[0] & 0x03F00)>>8;

  for(i=1; i<=nch; i++) {
    ch=(blockdata[i] & 0x3F0000)>>16;
    pdata[ch]=blockdata[i] & 0x0FFF;
  }    
  bk_close(pevent,pdata+N_ADC);  
  */

  /*
  bk_create(pevent, "TDC0", TID_DWORD, &pdata);

  for(i=0; i<N_TDC; i++) 
    pdata[i]=-1;

  nch=(blockdata[0] & 0x03F00)>>8;

  for(i=1; i<=nch; i++) {
    ch=(blockdata[i] & 0x3F0000)>>16;
    pdata[ch]=blockdata[i] & 0x0FFF;
  }    
  bk_close(pevent,pdata+N_TDC);  
  */

#endif

#ifndef USE_BLT_MODE
 /* create Charge ADC banks */
  if(dready_adc) {
    bk_create(pevent, "R_QD", TID_DWORD, &pdata);
    buf_addr=QADC_base << 16 | 0x0000;

    vme_A32D32_read(crate1, buf_addr, &d32);
    pdata[0]=d32;
    nch=(d32 & 0x03F00)>>8;
    for(i=0;i<nch+1;i++) {  //including EOB.
      vme_A32D32_read(crate1, buf_addr, &d32);
      pdata[i+1]=d32;
    }
    bk_close(pevent,pdata+34);  

    dready_adc=0;
  }
    //bk_create(pevent, "QADC", TID_DWORD, &pdata);

    //for(i=0; i<N_ADC; i++) 
    //pdata[i]=0;

    //vme_A32D32_read(crate1, buf_addr, &d32);
    //nch=(d32 & 0x03F00)>>8;
      
    //for(i=0;i<nch;i++) {
    //vme_A32D32_read(crate1, buf_addr, &d32);
    //ch=(d32 & 0x3F0000)>>16;
    //pdata[ch]=d32 & 0x0FFF;
    //}

    //vme_A32D32_read(crate1, buf_addr, &d32); /*read EOB.*/
    //bk_close(pevent,pdata+N_ADC);  

 /* create Peak ADC banks */
  if(dready_pdc) {
    bk_create(pevent, "R_PD", TID_DWORD, &pdata);

    buf_addr=PADC_base << 16 | 0x0000;

    vme_A32D32_read(crate1, buf_addr, &d32);
    pdata[0]=d32;
    nch=(d32 & 0x03F00)>>8;
    for(i=0;i<nch+1;i++) {
      vme_A32D32_read(crate1, buf_addr, &d32);
      pdata[i+1]=d32;
    }
    bk_close(pevent,pdata+34);  


    dready_pdc=0;
  }

 /* create Peak ADC banks */
  if(dready_pdc2) {
    bk_create(pevent, "R_P2", TID_DWORD, &pdata);

    buf_addr=PDC2_base << 16 | 0x0000;

    vme_A32D32_read(crate1, buf_addr, &d32);
    pdata[0]=d32;
    nch=(d32 & 0x03F00)>>8;
    for(i=0;i<nch+1;i++) {
      vme_A32D32_read(crate1, buf_addr, &d32);
      pdata[i+1]=d32;
    }
    bk_close(pevent,pdata+34);  


    dready_pdc2=0;
  }

 /* create Peak ADC3 banks */
  if(dready_pdc3) {
    bk_create(pevent, "R_P3", TID_DWORD, &pdata);

    buf_addr=PDC3_base << 16 | 0x0000;

    vme_A32D32_read(crate1, buf_addr, &d32);
    pdata[0]=d32;
    nch=(d32 & 0x03F00)>>8;
    for(i=0;i<nch+1;i++) {
      vme_A32D32_read(crate1, buf_addr, &d32);
      pdata[i+1]=d32;
    }
    bk_close(pevent,pdata+34);  


    dready_pdc3=0;
  }

 /* create TDC banks */
  if(dready_tdc) {
    bk_create(pevent, "R_TD", TID_DWORD, &pdata);

    buf_addr=TDC_base << 16 | 0x0000;

    vme_A32D32_read(crate1, buf_addr, &d32);
    pdata[0]=d32;
    nch=(d32 & 0x03F00)>>8;
    for(i=0;i<nch+1;i++) {
      vme_A32D32_read(crate1, buf_addr, &d32);
      printf("%d\n",(d32>>24)&0x3);
      pdata[i+1]=d32;
    }
    bk_close(pevent,pdata+34);  

    dready_tdc=0;
  }
  //buf_addr=TDC_base << 16 | 0x0000;

  //bk_create(pevent, "TDC0", TID_DWORD, &pdata);

  //for(i=0; i<N_TDC; i++) 
    //pdata[i]=4096;//0;

  //vme_A32D32_read(crate1, buf_addr, &d32);
  //nch=(d32 & 0x03F00)>>8;
  ////valid=(d32 & 0x7000000)>>24;
  ////if(valid==6) return;
      
  //for(i=0;i<nch;i++) {
    //vme_A32D32_read(crate1, buf_addr, &d32);
    //ch=(d32 & 0x3F0000)>>16;
    //pdata[ch]=d32 & 0x0FFF;
  //}

  //vme_A32D32_read(crate1, buf_addr, &d32); /*read EOB.*/
  //temp=d32 & 0xFFFFFF;
  //bk_close(pevent,pdata+N_TDC);  

    ////bk_create(pevent, "TDEN", TID_DWORD, &eventn);
    ////    vme_A32D16_read(crate1, TDC_base<<16 | 0x1024, &d32);
    //*eventn++=temp;
    //bk_close(pevent,eventn);
#endif

/*   /\* create SIS0 bank *\/ */
/*   bk_create(pevent, "SIS0", TID_DWORD, &pdata); */
/*   vme_A32D32_read(crate1, IO_base<<8 | 0x0, &d32); */
/*   /\* */
/*   i=0; */
/*   while(!(d32 & 0x100)) { */
/*     i++; */
/*     vme_A32D32_read(crate1, IO_base << 8 | 0x100, &temp); */
/*     vme_A32D32_read(crate1, IO_base<<8 |0x0, &d32); */
/*   } */
/*   *pdata++=temp; */
/*   *pdata++=i; */
/*   *\/ */
  
/*   if(d32 & 0x100) { */
/*     ss_printf(0,17,"SIS3600 FIFO empty!\n"); */
/*     *pdata++=-1; */
/*   } */
/*   else { */
/*     vme_A32D32_read(crate1, IO_base << 8 | 0x100, &d32); */
/*     *pdata++ = d32; */
/*   } */
  
/*   bk_close(pevent, pdata); */

  /*software trigger for V830 scaler*/
  //vme_A32D16_write(crate1, SCA830_base << 16 | 0x1124, 0x01);

/*   /\* create SCLR830 bank *\/ */
/*   bk_create(pevent, "S830", TID_DWORD, &pdata); */

/* /\*   for(i=0; i<32; i++) { *\/ */
/* /\*     vme_A32D32_read(crate1, SCA830_base << 16 | (0x1000+i*4), &d32); *\/ */
/* /\*     *pdata++ = d32; *\/ */
/* /\*   } *\/ */
/*   for(i=0; i<32; i++) { */
/*     vme_A32D32_read(crate1, SCA830_base << 16 | (0x0000), &d32); */
/*     *pdata++ = d32; */
/*   } */
		  
/*   bk_close(pevent, pdata); */

  /*software trigger*/
  //vme_A32D16_write(crate1, SCA820_base << 16 | 0x1124, 0x01);
  
/*   /\* create SCLR 820 bank *\/ */
/*   bk_create(pevent, "S820", TID_DWORD, &pdata); */
/*   for(i=0; i<32; i++) { */
/*     //jianglai's hack for now. only read out the first channel in 820 */
/*     //to minimize deadtime */
/*     if(i==0) vme_A32D32_read(crate1, SCA820_base << 16 | (0x1000+i*4), &d32); */
/*     else d32=0; */
/*     *pdata++ = d32; */
/*   } */
/*   bk_close(pevent, pdata); */

  //end_read = ss_millitime();
  gettimeofday(&tv, &tz);
  tm=localtime(&tv.tv_sec);
  end_read = tv.tv_usec;
  return bk_size(pevent);
}

/*-- Scaler event --------------------------------------------------*/

INT read_scaler_event(char *pevent, INT off)
{
  WORD d2,b,i;
  DWORD *pdata, a;

  /* init bank structure */
  bk_init(pevent);

  /*software trigger*/
  //vme_A32D16_write(crate1, SCA820_base << 16 | 0x1124, 0x01);

  /* create SCLR 820 bank */
  bk_create(pevent, "S820", TID_DWORD, &pdata);

  for(i=0; i<32; i++) {
    vme_A32D32_read(crate1, SCA820_base << 16 | (0x1000+i*4), &a);
    *pdata++ = a;
    if (i < 4) {
      ss_printf(25+i*12,9,"%10d",a);
    }
    else if (i < 8) {
      ss_printf(25+(i-4)*12,10,"%10d",a);
    }
     else if (i < 12) {
      ss_printf(25+(i-8)*12,11,"%10d",a);
     }
     else if (i < 16) {
      ss_printf(25+(i-12)*12,12,"%10d",a);
     }
  }

  bk_close(pevent, pdata);

  /*software trigger for V830 scaler*/
  //vme_A32D16_write(crate1, SCA830_base << 16 | 0x1124, 0x01);

/*   /\* create SCLR 830 bank *\/ */
/*   bk_create(pevent, "S830", TID_DWORD, &pdata); */
/*   for(i=0; i<32; i++) { */
/*     vme_A32D32_read(crate1, SCA830_base << 16 | (0x1000+i*4), &a); */
/*     *pdata++ = a; */
/*   } */
/*   bk_close(pevent, pdata); */


  /* create SIS0 bank */
  bk_create(pevent, "SIS0", TID_DWORD, &pdata);
  vme_A32D32_read(crate1, IO_base<<8 |0x0, &a);

  if(a & 0x100) {
    ss_printf(0,17,"SIS3600 FIFO empty!\n");
    *pdata++=-1;
  }
  else {
    vme_A32D32_read(crate1, IO_base << 8 | 0x100, &a);
    *pdata++ = a;
  }

  bk_close(pevent, pdata);

  return bk_size(pevent);
}

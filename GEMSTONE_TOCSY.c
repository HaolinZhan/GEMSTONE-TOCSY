/* GEMSTONE: Gradient-Enhanced Multiplet-Selective Targeted-Observation NMR Experiment

// Peter Kiraly, Nicolas Kern, Mateusz P. Plesniak, Mathias Nilsson, David J. Procter, Gareth A. Morris, and Ralph William Adams 
// Angew.Chem. (2020) https://doi.org/10.1002/anie.202011642

_r1 1st release, same as used for publication  
*/


#include <standard.h>
#include <chempack.h>

static int	ph1 [8] = {0,0, 0,0, 0,0,0,0},	//v1: excitation
		ph2 [8] = {0,0, 0,0, 1,1,1,1},	//v2: 1st 180ad
		ph3 [8] = {0,1, 0,1, 0,1,0,1},	//v3: sel180
		ph4 [8] = {0,0, 1,1, 0,0,1,1},	//v4: 2nd 180ad
		ph5 [8] = {2,0, 0,2, 0,2,2,0},	//oph
		ph6 [8] = {0,0, 0,0, 0,0,0,0};	//v5: non-sel90 during ZQF-TOCSY  
		

pulsesequence()
{

double	tpwr = getval("tpwr"),
	pw = getval("pw"),

	pw180_e = getval("pw180_e"),
	pwr180_e = getval("pwr180_e"),


//sel180 refocusing of the active spin
	pw180_a = getval("pw180_a"),
	pwr180_a = getval("pwr180_a"),
// wurst pulses
	pw180_bb1 = getval("pw180_bb1"),
	pwr180_bb1 = getval("pwr180_bb1"),
	pw180_bb2 = getval("pw180_bb2"),
	pwr180_bb2 = getval("pwr180_bb2"),
//PFG parameters
	hsgt = getval("hsgt"),
	hsglvl = getval("hsglvl"),

	gt1 = getval("gt1"),  
	gzlvl1 = getval("gzlvl1"),
	gt2 = getval("gt2"),  
	gzlvl2 = getval("gzlvl2"),
	gt3 = getval("gt3"),  
	gzlvl3 = getval("gzlvl3"),
	gt4 = getval("gt4"),  
	gzlvl4 = getval("gzlvl4"),

	gstab = getval("gstab"),
	gzlvl7 = getval("gzlvl7"),
	gzlvl8 = -1.0*gzlvl7, //getval("gzlvl8"), //gzlvl7*(pw180_bb1/pw180_bb2),
	gstab0 = getval("gstab0"),

//new add
  	slpwrT = getval("slpwrT"),
        slpwT = getval("slpwT"),
  	mixT = getval("mixT"),
        zqfpw1 = getval("zqfpw1"),
        zqfpwr1 = getval("zqfpwr1"),
        zqfpw2 = getval("zqfpw2"),
        zqfpwr2 = getval("zqfpwr2"),
        gzlvlzq1 = getval("gzlvlzq1"),
        gzlvlzq2 = getval("gzlvlzq2");


int kpph = getval("kpph");
char	sspul[MAXSTR],	shp_a[MAXSTR], shp_bb1[MAXSTR], shp_bb2[MAXSTR], shp_e[MAXSTR], slpatT[MAXSTR], zqfpat1[MAXSTR], zqfpat2[MAXSTR];
char kp_cssf[MAXSTR];

	getstr("kp_cssf",kp_cssf);       
	getstr("shp_e",shp_e);		
  
	getstr("sspul",sspul);
	getstr("shp_a",shp_a);		//created by macro kp_makePS9(2) using bw_a/kp_phincr_a
	getstr("shp_bb1",shp_bb1);
	getstr("shp_bb2",shp_bb2);   

   getstr("slpatT",slpatT);
   getstr("zqfpat1",zqfpat1);
   getstr("zqfpat2",zqfpat2);

   if (strcmp(slpatT,"mlev17c") &&
        strcmp(slpatT,"dipsi2") &&
        strcmp(slpatT,"dipsi3") &&
        strcmp(slpatT,"mlev17") &&
        strcmp(slpatT,"mlev16"))
        abort_message("SpinLock pattern %s not supported!.\n", slpatT);     


if (kpph==0)
{
	settable(t1,8,ph1);
	settable(t2,8,ph2);
	settable(t3,8,ph3);
	settable(t4,8,ph4);
	settable(t5,8,ph5);
	settable(t6,8,ph6);
}
else
{
	settable(t1,kpph,ph1);
	settable(t2,kpph,ph2);
	settable(t3,kpph,ph3);
	settable(t4,kpph,ph4);
	settable(t5,kpph,ph5);
	settable(t6,kpph,ph6);
}
	getelem(t1,ct,v1);
	getelem(t2,ct,v2);
	getelem(t3,ct,v3);
	getelem(t4,ct,v4);
	getelem(t6,ct,v6);
	getelem(t5,ct,oph);



status(A);
	txphase(zero);
        obsoffset(tof);
	obspower(tpwr);
	delay(0.001);

	if (sspul[0]=='y') 
	{
	zgradpulse(hsglvl,hsgt); 
	delay(gstab);
	rgpulse(pw,zero,rof1,rof1);
	zgradpulse(hsglvl,hsgt); 
	delay(gstab);
	}
lk_sample();
	delay(d1);
lk_hold();

delay(0.001);

status(B);

rgpulse(pw, v1, rof1, 0.0);



	obspower(pwr180_bb1);
	delay(gstab0);
	delay(rof2);
	rgradient('z',gzlvl7);
	shaped_pulse(shp_bb1,pw180_bb1,v2,rof1,rof1);	
	rgradient('z',0.0);
	delay(gstab0);

	obspower(pwr180_a);
        zgradpulse(gzlvl2,gt2);
/*	rgradient('x',gzlvl2);
	delay(gt2);
	rgradient('x',0.0);	
*/
        delay(gstab);
	shaped_pulse(shp_a,pw180_a,v3,rof1,rof1);	
        zgradpulse(gzlvl2,gt2);
/*
	rgradient('x',gzlvl2);
	delay(gt2);
	rgradient('x',0.0);	
*/
        delay(gstab);


	obspower(pwr180_bb2);
	delay(gstab0);
	rgradient('z',1.0*gzlvl8);
	shaped_pulse(shp_bb2,pw180_bb2,v4,rof1,rof1);	
	rgradient('z',0.0);
	delay(gstab0);

status(C);
	obspower(tpwr);
	rgpulse(pw,v6,rof1,rof1);

        if (mixT > 0.0)
        {
           if (getflag("Gzqfilt"))
           {
            obspower(zqfpwr1);
            rgradient('z',gzlvlzq1);
            delay(100.0e-6);
            shaped_pulse(zqfpat1,zqfpw1,zero,rof1,rof1);
            delay(100.0e-6);
            rgradient('z',0.0);
            delay(100e-6);
           }
           obspower(slpwrT);
           zgradpulse(gzlvl3,gt3);
           delay(gt3);
        
           if (dps_flag)
             rgpulse(mixT,v4,0.0,0.0);
           else
             SpinLock(slpatT,mixT,slpwT,v4);

           if (getflag("Gzqfilt"))
           {
            obspower(zqfpwr2);
            rgradient('z',gzlvlzq2);
            delay(100.0e-6);
            shaped_pulse(zqfpat2,zqfpw2,zero,rof1,rof1);
            delay(100.0e-6);
            rgradient('z',0.0);
            delay(100e-6);
           }
           obspower(tpwr);
           zgradpulse(gzlvl4,gt4);
           delay(gt4);
         }

        rgpulse(pw,v6,rof1,rof2);

status(D);

	startacq(getval("alfa"));
	sample(getval("at"));
	recoff();	
	endacq();


delay(0.05);	
lk_sample();
delay(0.05);
}


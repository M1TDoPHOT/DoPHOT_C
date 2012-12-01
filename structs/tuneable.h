#ifndef TUNEABLE_H_
#define TUNEABLE_H_

#define NSMAX    20000 // max number of allowable stars
#define NPMAX    11    // max number of allowable star parameters
#define MMAX     14    // needed by chisq, improve, isearch, shape, twofit
                       // = max npars+3 for max number of pars of two star fit
#define NPSKYMAX 7     // max number of allowable sky  parametera
#define NPSKY    7     // number of parameters in plane sky model
#define NSKYFIT  3     
#define NPHUB    7     // number of parameters in hubble sky model
#define NHUBFIT  7
#define NFF      20    // number of flags or files
#define NEP      7     // number of emperical parameters
#define IHSIDE   32    // side length of emperical box
#define NAPPLE   5     // number of APPLE parameters (apperature fit params)
#define MAXFIL   1024  //needed by subraster struct to set size of arrays 

#define nphead   400   // max number of header lines in parameter files

#define pi 3.14159265f
#define TINY 1.0e-20
#define OBLITVAL -60000 
#define MAGIC 2147483647

struct{
     float gxwid, gywid, skyguess;
}tune1_;
// gxwid is the psf fwhm width^2 times a const.

struct{
//     short int irect[2], krect[2];
//     float arect[2];
     short int* irect; //length 2
     short int* krect; //length 2
     float* arect; //length 2
}tune2_;
// irect holds x and y fixbox values
// arect holds x and y aperature box values
// krect NFITBOXFIRST_X, Y

struct{
    float tmin, tmax, tfac;
    int ibot, itop;
    int dofinalfit;
}tune3_;
/*ITOP- Saturated pixel value */
/*IBOT- Lowest allowed pixel value */
/*TMIN- is the minimum threshold for the sky THRESHMIN */
/*TMAX- THRESHMAX */
/*TFAC- THRESHDEC */
/* THRESHDEC The separation between successive thresholds in units of 2.0 raised to the power THRESHDEC. Thus if this parameter is 1.0, and THRESHMIN = 20, and THRESHMAX = 1000, the thresholds (relative to the sky value) would be 640, 320, 160, 80, 40, and 20. */
/*DOFINALFIT Do a last fitting iteration after all objects found.  This extra loop will ensure that objects of type 3 found ont he final threshold pass are fully fit, but it will slow DoPHOT */

struct{
    int nit, npar, nfit1, nfit2;
}tune4_;
/*NIT- NFITITER */
/*NPAR- NPARAM */
/*NFIT1- NFITMAG */
/*NFIT2- NFITSHAPE */

struct{
    float fac, xpnd, gxpnd;
}tune5_;
/*FAC- RESIDNOISE */
/*XPND- FOOTPRINT_NOISE */

struct{
    int nphsub, nphob, icrit;
    float ctpersat, widobl, cmax;
}tune6_;
/*NPHSUB – When subtracting an object from an image, DoPHOT goes out to where the PSF fitted to an object reaches this surface brightness. */
/*NPHOB – When obliterating objects, DoPHOT determines where the PSF that has been fitted to an ‘obliteratable’ object (see the discussion of ICRIT, CENTINTMAX and CTPER- SAT below) reaches a surface brightness equal to NPHOB, in DN. The size of the side of an obliteration box for a given object is then set equal to twice this linear distance. */
/*ICRIT – If there are this many contiguous saturated pixels (as determined by testing the pixel values with ITOP), then DoPHOT obliterates the object.*/
/*CTPERSAT – If an object is to be obliterated then its central intensity is assumed to be CTPERSAT times the number of saturated pixels in the object. This helps DoPHOT calculate reasonably-sized obliteration boxes.*/
/*CMAX - CENTINTMAX*/
/*WIDOBL - COSOBLSIZE*/

struct{
    int n0left, n0right;
}tune7_;

struct{
    float stograt, discrim;
//    float sig[3];
    float* sig; //length 3
}tune8_;
/*stograt- STARGALKNOB is the STar Over Galaxy chisq RATio which sets the threshold for when large objects are double stars or galaxies in the function SHAPE */
/*discrim- STARCOSKNOB */
/*sig[] = SIMGA1, SIGMA2, SIGMA3 */

struct{
    float chicrit, xtra, crit7, snlim, bumpcrit, sn2cos;
}tune9_;
/*CHICRIT- CHI2MINBIG */
/*XTRA- XTRA */
/*CRIT7- SNLIM7 */
/*SNLIM- SNLIM */
/*BUMPCRIT- SNLIMMASK */
/*sn2cos- SNLIMCOS is cosmic ray signal to noise */

struct{
    float enuff4, enuff7;
}tune10_;
/*ENUFF4, ENUFF7 – At least ENUFF4 of the pixels in a fit subraster must be present in or- der to perform a four-parameter fit (see NFITMAG above). ENUFF7 specifies this limit for seven-parameter fits*/

struct{
    float eperdn, rnoise;
}tune11_;
/*EPERDN – Electrons per DN for the detector.*/
/*RDNOISE – The readout noise of the detector in electrons.*/

struct{
    int ixby2, iyby2;
}tune12_;
// maskbox x and y radii

struct{
    int nbadleft, nbadright, nbadtop, nbadbot;
}tune13_;
/*NBADLEFT, NBADRIGHT, NBADTOP, NBADBOT – Pixels that are located closer to the edges of the image than these values are ignored by DoPHOT – sort of like ‘software trim- ming’. Left, right, top, and bottom correspond to small x, large x, large y, and small y, respectively.*/

struct{
    int lverb;
}tune14_;
/*lverb- Specifies the logfile's verbosity.  Higher equates to more reporting in the logfile */

struct{
//    float acc[NPMAX], alim[NPMAX], ava[NPMAX];
    float* acc; //length NPMAX
    float* alim; //length NPMAX
    float* ava; //length NPMAX
    int max_perf; 
}tune15_;
/*ACC- RELACC */
/*ALIM- ABSLIM */
/*AVA- guess at average values? */
/*max_perf- maximum number of perfect stars to be used for computing average shape parameters.  useful as dimmer stars are often problematic beyond what their noise might indicate. */ 

struct{
//    char flags[NFF][?], files[NFF][?];
    char** flags;
    char** files;
}tune16_;

struct{
    float beta4, beta6;
}tune17_;
/*BETA4, BETA6  These parameters modify the z4 and z6 terms in the pseudo-gaussian PSF function*/

struct{
    float pixthresh;
    int nthpix;
}tune18_;
/*PIXTHRESH – This parameter allows DoPHOT’s finding algorithm to trigger on pixels that are above PIXTHRESH times the local noise. This parameter represents one of the masking tape/baling wire type fixes that has been incorporated into the code in order to allow the finding algorithm to find faint stars in well-sampled data. By setting PIXTHRESH too low, however, you force the program waste time trying to find objects where only a single noisy pixel exists. Binning well-sampled data and leaving PIXTHRESH at 1.0 is usually more effective at finding faint stars.*/
/*NTHPIX – This is really not a filtering parameter, but one used in the search for objects. It controls how often the search routine updates the value of the background from the model (be it PLANE, HUBBLE or MEDIAN) when initially raster searching along X for peaks in the data. NTHPIX is the interval in pixels for this update. Setting NTHPIX 0 or negative (the default) sets NTHPIX = SQRT(X-size of image in pixels) to the near- est integer.*/

struct{
    float apmagmaxerr;
}tune19_;
/*APMAGMAXERR – Only aperture photometry results with errors smaller than the value specified by this parameter are reported in the output photometry file.*/

struct{
    int jhxwid, jhywid, mprec;
}tune20_;
/* for SKYTYPE = 'MEDIAN' only */
/*JHXWID, JHYWID – These specify the HALF size of the median filtering box along X and Y respectively. They must be integer values. If you specify a non-positive value, the size will be scaled automatically from the FWHM parameter. The default for both is 0, which forces auto-scaling.*/
/*MPREC – This controls the precision to which the median model is calculated. If the noise in the image background is say 100 DN, you do not need a model computed to 1 DN precision. You can ask instead for precision to n DN by setting the value of the MPREC parameter to n. This must be an integer.*/

struct{
    int fixpos;
}tune21_;
/*FIXPOS – This flag can only have the values YES or NO. Any illegal or unitelligible value for this keyword in the parameter file will cause this parameter to default to NO. If FIX- POS is YES, then the program assumes you have an input template file with the posi- tions of objects as described in detail in Section 11.1 and 11.2 of the MANUAL*/

struct{
    float eminsig;
    int ihtab;
    int empok;
    float emthrsh;
    int emenab;
    float efac, gfac;
}tune22_;
/*ihtab- NEMP_PSF_BOX*/
/*eminsig- EMP_REJ_RAD_SIG*/
/*efac- EMP_RESIDNOISE*/
/*emthrsh- THRESHEMP*/

struct{
    int nempski;
    float targx, targy, targz;
}tune23_;
/*nempski- N_EMP_SKIP is the number of bright template stars to skip in empirical model*/ 
/*TARGX Y Z = EMP_STAR_X Y X */

#endif 

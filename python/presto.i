%module presto

%{
#include "presto.h"
#include "ransomfft.h"
%}

%include numpy.i

%{

/* Note the following is copied and slightly modified from Nick 
   Patavalis' PPGPLOT python module.                            */

static PyObject *
arrtofloatvector (PyObject *o, float **v, int *vsz){
/* arrtofloatvector(array, pointer to data, number of points):
      Convert a 1-D Numeric array into single-precision floats.
      NOTE:  You could lose some accuracy!!!  */

    PyArrayObject *a1, *af1, *af2;
    int ownedaf1=0;
    
    /* Check if args are arrays. */
    if (!PyArray_Check(o)) {
	PyErr_SetString(PyExc_TypeError,\
			"tofloatvector() object is not an array");
	return(NULL);
    }
    a1 = (PyArrayObject *)o;

    /* Check if args are vectors. */
    if (a1->nd != 1) {
	PyErr_SetString(PyExc_TypeError,\
			"tofloatvector() object is not a vector");
	return(NULL);
    }
    
    switch (a1->descr->type_num) {
    case PyArray_FLOAT:
	af1 = a1;
	break;
    case PyArray_CHAR: 
    case PyArray_UBYTE: 
    case PyArray_SBYTE:
    case PyArray_SHORT: 
    case PyArray_INT: 
    case PyArray_LONG:
    case PyArray_DOUBLE:
	if (!(af1 = (PyArrayObject *)PyArray_Cast(a1,PyArray_FLOAT))) {
	    PyErr_SetString(PyExc_TypeError,\
			    "tofloatvector() cannot cast vector to floats");
	    return(NULL);
	}
	ownedaf1 = 1;
	break;
    default:
	PyErr_SetString(PyExc_TypeError,\
			"tofloatvector() cannot cast vector to floats");
	return(NULL);
	break;
    }
    
    af2 = af1;
    if (PyArray_As1D((PyObject **)&af2, (char **)v, vsz, 
		     PyArray_FLOAT) == -1) {
	af2 = NULL;
    }
    
    if (ownedaf1)
	Py_DECREF(af1);
    return((PyObject *)af2);
}


static PyObject* 
tofloatvector(PyObject* self, PyObject* args) {
/* tofloatvector(array):
      Convert a 1-D Numeric array into single-precision floats.
      NOTE:  You could lose some accuracy!!!  */
  PyObject *arr;
  float *data=NULL;
  int numdata=0;

  if (!PyArg_ParseTuple(args, "O", &arr))
    return NULL;
  return (PyObject *)arrtofloatvector(arr, &data, &numdata);
}

static PyObject* 
float_to_complex(PyObject* self, PyObject* args) {
/* float_to_complex(array):
      Convert a 1-D Numeric float array into a 1-D
      single-precision complex array.  (i.e. assume that
      the float array really held the real and complex values
      of a complex array)                                      */
  PyArrayObject *arr;
  float *data=NULL;
  int numdata=0;

  if (!PyArg_ParseTuple(args, "O", &arr))
    return NULL;
  /* Check that obj is really an array */
  if (!PyArray_Check(arr)) {
    PyErr_SetString(PyExc_TypeError,\
		    "float_to_complex() Argument is not an array");
    return NULL;
  }
  /* Check the type */
  if (arr->descr->type_num != PyArray_FLOAT) {
    PyErr_SetString(PyExc_TypeError,\
		    "float_to_complex() Array is not made of floats.");
    return NULL;
  }
  /* Insure that the array is really 1-D */
  if (arr->nd != 1) {
    PyErr_SetString(PyExc_TypeError,\
		    "float_to_complex() Array is not 1-D.");
    return NULL;
  }
  /* Insure that the array has an even number of points */
  if (arr->dimensions[0] & 1) {
    PyErr_SetString(PyExc_TypeError,\
		    "float_to_complex() Array needs an even number of points.");
    return NULL;
  }
  /* Perform the conversion  (NOTE:  This is probably _way_ unsafe...)*/
  arr->descr = PyArray_DescrFromType(PyArray_CFLOAT);
  arr->dimensions[0] >>= 1;
  arr->strides[0] = arr->descr->elsize;
  Py_INCREF(Py_None);
  return Py_None;
}


static PyObject* 
complex_to_float(PyObject* self, PyObject* args) {
/* complex_to_float(array):
      Convert a 1-D Numeric complex array into a 1-D
      single-precision float array holding the real and imaginary
      values in contiguous memory locations. */
  PyArrayObject *arr;
  float *data=NULL;
  int numdata=0;

  if (!PyArg_ParseTuple(args, "O", &arr))
    return NULL;
  /* Check that obj is really an array */
  if (!PyArray_Check(arr)) {
    PyErr_SetString(PyExc_TypeError,\
		    "complex_to_float() Argument is not an array");
    return NULL;
  }
  /* Check the type */
  if (arr->descr->type_num != PyArray_CFLOAT) {
    PyErr_SetString(PyExc_TypeError,\
		    "complex_to_float() Array is not made of complex #'s.");
    return NULL;
  }
  /* Insure that the array is really 1-D */
  if (arr->nd != 1) {
    PyErr_SetString(PyExc_TypeError,\
		    "complex_to_float() Array is not 1-D.");
    return NULL;
  }
  /* Perform the conversion  (NOTE:  This is probably _way_ unsafe...)*/
  arr->descr = PyArray_DescrFromType(PyArray_FLOAT);
  arr->dimensions[0] <<= 1;
  arr->strides[0] = arr->descr->elsize;
  Py_INCREF(Py_None);
  return Py_None;
}

float *power_arr(fcomplex *dft, long numfreqs){
  /* Determine the spectral powers of the Fourier amplitudes 'dft'*/
  float powargr, powargi, *powers;
  long i;

  powers = gen_fvect(numfreqs);
  for (i=0; i<numfreqs; i++)
    powers[i] = POWER(dft[i].r, dft[i].i);
  return powers;
}

float *phase_arr(fcomplex *dft, long numfreqs){
  /* Determine the spectral phases of the Fourier amplitudes 'dft'*/
  float phsargr, phsargi, phstmp, *phases;
  long i;

  phases = gen_fvect(numfreqs);
  for (i=0; i<numfreqs; i++)
    phases[i] = PHASE(dft[i].r, dft[i].i);
  return phases;
}

double *dpower_arr(dcomplex *dft,  long numfreqs){
  /* Determine the spectral powers of the Fourier amplitudes 'dft'*/
  double powargr, powargi, *powers;
  long i;

  powers = gen_dvect(numfreqs);
  for (i=0; i<numfreqs; i++)
    powers[i] = POWER(dft[i].r, dft[i].i);
  return powers;
}

double *dphase_arr(dcomplex *dft, long numfreqs){
  /* Determine the spectral phases of the Fourier amplitudes 'dft'*/
  double phsargr, phsargi, phstmp, *phases;
  long i;

  phases = gen_dvect(numfreqs);
  for (i=0; i<numfreqs; i++)
    phases[i] = PHASE(dft[i].r, dft[i].i);
  return phases;
}

%}

%native(tofloatvector) PyObject *tofloatvector(PyObject *, PyObject *);
/* tofloatvector(arr):
      Convert the vector 'arr' into a vector of single-precision floats. */
%native(float_to_complex) PyObject *float_to_complex(PyObject *, PyObject *);
/* float_to_complex(array):
      Convert a 1-D Numeric float array into a 1-D
      single-precision complex array.  (i.e. assume that
      the float array really held the real and complex values
      of a complex array)                                      */
%native(complex_to_float) PyObject *complex_to_float(PyObject *, PyObject *);
/* complex_to_float(array):
      Convert a 1-D Numeric complex array into a 1-D
      single-precision float array holding the real and imaginary
      values in contiguous memory locations. */


#ifndef SQRT2
#define SQRT2         1.4142135623730950488016887242096980785696718753769
#endif
#ifndef PI
#define PI            3.1415926535897932384626433832795028841971693993751
#endif
#ifndef TWOPI
#define TWOPI         6.2831853071795864769252867665590057683943387987502
#endif
#ifndef DEGTORAD
#define DEGTORAD      0.017453292519943295769236907684886127134428718885417
#endif
#ifndef RADTODEG
#define RADTODEG      57.29577951308232087679815481410517033240547246656
#endif
#ifndef PIBYTWO
#define PIBYTWO       1.5707963267948966192313216916397514420985846996876
#endif
#ifndef SOL
#define SOL           299792458.0
#endif
#ifndef SECPERJULYR
#define SECPERJULYR   31557600.0
#endif
#ifndef SECPERDAY
#define SECPERDAY     86400.0
#endif
#ifndef ARCSEC2RAD
#define ARCSEC2RAD    4.8481368110953599358991410235794797595635330237270e-6
#endif
#ifndef SEC2RAD
#define SEC2RAD       7.2722052166430399038487115353692196393452995355905e-5
#endif

typedef enum {
     LOWACC = 0, HIGHACC = 1
} presto_interp_acc;

typedef enum {
     CONV = 0, CORR = 1, INPLACE_CONV = 2, INPLACE_CORR = 3
} presto_optype;

typedef enum {
     FFTDK = 0, FFTD = 1, FFTK = 2, NOFFTS = 3
} presto_ffts;

typedef enum {
     RAW = 0, PREPPED = 1, FFT = 2, SAME = 3
} presto_datainf;

%apply fcomplex* IN_1D_CFLOAT { fcomplex *dft };
%apply long ARRAYLEN { long numfreqs };
float *power_arr(fcomplex *dft, long numfreqs);
/* power_arr(dft, powers, numfreqs):
      Return the power spectrum of a complex FFT 'dft' in 'powers'. */

%apply fcomplex* IN_1D_CFLOAT { fcomplex *dft };
%apply long ARRAYLEN { long numfreqs };
float *phase_arr(fcomplex *dft, long numfreqs);
/* phase_arr(dft, phases, numfreqs):
      Return the spectral phases (deg) of a complex FFT 'dft' in 'phases'. */

%apply dcomplex* IN_1D_CDOUBLE { dcomplex *dft };
%apply long ARRAYLEN { long numfreqs };
double *dpower_arr(dcomplex *dft, long numfreqs);
/* dpower_arr(dft, powers, numfreqs):
      Return the power spectrum of a complex FFT 'dft' in 'powers'. */

%apply dcomplex* IN_1D_CDOUBLE { dcomplex *dft };
%apply long ARRAYLEN { long numfreqs };
double *dphase_arr(dcomplex *dft, long numfreqs);
/* dphase_arr(dft, phases, numfreqs):
      Return the spectral phases (deg) of a complex FFT 'dft' in 'phases'. */


typedef struct orbitparams {
  double p;	    /* Orbital period (s)                            */
  double e;	    /* Orbital eccentricity                          */
  double x;	    /* Projected semi-major axis (lt-sec)            */
  double w;	    /* Longitude of periapsis (deg)                  */
  double t;	    /* Time since last periastron passage (s)        */
  double pd;	    /* Orbital period derivative (s/yr)              */
  double wd;	    /* Advance of longitude of periapsis (deg/yr)    */
  %addmethods {
      orbitparams(){
	orbitparams *p = (orbitparams *)malloc(sizeof(orbitparams));
	return p;
      }
      ~orbitparams(){
	free(self);
      }
  }
} orbitparams;

/* These files describe makedata, makeinf, and PSR database structs  */
%include presto_psrdata.i
%include presto_structs.i

typedef struct RDERIVS {
  float pow;	   /* Power normalized with local power             */
  float phs;	   /* Signal phase                                  */
  float dpow;	   /* 1st deriv of power wrt fourier freq           */
  float dphs;	   /* 1st deriv of phase wrt fourier freq           */
  float d2pow;     /* 2nd deriv of power wrt fourier freq           */
  float d2phs;     /* 2nd deriv of power wrt fourier freq           */
  float locpow;    /* Local mean power level                        */
  %addmethods {
      rderivs(){
	rderivs *p = (rderivs *)malloc(sizeof(rderivs));
	return p;
      }
      ~rderivs(){
	free(self);
      }
  }
} rderivs;

typedef struct FOURIERPROPS {
  double r;        /* Fourier freq at start of run                  */
  float rerr;	   /* Error in fourier freq                         */
  double z;        /* Fourier freq deriv at start of run = T^2*fdot */
  float zerr;      /* Error in fourier freq deriv                   */
  double w;        /* Fourier 2nd freq deriv = T^3*fdotdot          */
  float werr;      /* Error in 2nd fourier freq deriv               */
  float pow;	   /* Power normalized with local power             */
  float powerr;    /* Error in power                                */
  float sig;	   /* Approx. Sigma level                           */
  float rawpow;    /* Raw signal power (unnormalized)               */
  float phs;	   /* Signal phase                                  */
  float phserr;    /* Error in phase                                */
  float cen;	   /* Centroid of signal pulsations over run        */
  float cenerr;    /* Error in centroid                             */
  float pur;	   /* Purity:  rms duration of signal over run      */
  float purerr;    /* Error in purity                               */
  float locpow;    /* Local mean power level                        */
  %addmethods {
    fourierprops(){
      fourierprops *p = (fourierprops *)malloc(sizeof(fourierprops));
      return p;
    }
    ~fourierprops(){
      free(self);
    }
  }
} fourierprops;

typedef struct BINARYPROPS {
  double ppsr;       /* Pulsar period (s)                            */
  double fpsr;       /* Pulsar freq (hz)                             */
  double rpsr;       /* Pulsar Fourier freq (bins)                   */
  double pbin;       /* Binary period (s)                            */
  double rbin;       /* Binary Fourier freq (bins)                   */
  double z;          /* Amplitude of frequency modulation            */
  double asinic;     /* Projected semi-major axis of orbit (lt-sec)  */
  double rdetect;    /* Mini-FFT candidate detected freq (bins)      */
  long nfftbins;     /* Number of FFT bins in mini-fft               */
  long lowbin;       /* Start bin of mini-FFT in original FFT        */
  float ppsrerr;     /* Error in pulsar period (s)                   */
  float fpsrerr;     /* Error in pulsar freq (hz)                    */
  float rpsrerr;     /* Error in pulsar Fourier freq (bins)          */
  float pbinerr;     /* Error in binary period (s)                   */
  float rbinerr;     /* Error in pulsar Fourier freq (bins)          */
  float zerr;        /* Error in freq modulation ampltude            */
  float asinicerr;   /* Error in a*sin(i)/c in (lt-sec)              */
  float rdetecterr;  /* Error in Mini-FFT detected freq (bins)       */
  float sig;	     /* Approx. Sigma level                          */
  float phs;	     /* Signal phase                                 */
  float phserr;      /* Error in phase                               */
  float cen;	     /* Centroid of signal pulsations over run       */
  float cenerr;      /* Error in centroid                            */
  float pur;	     /* Purity:  rms duration of signal over run     */
  float purerr;      /* Error in purity                              */
  float pow;         /* Power normalized with Nph equivalent         */
  float powerr;      /* Error in signal power                        */
  %addmethods {
    binaryprops(){
      binaryprops *p = (binaryprops *)malloc(sizeof(binaryprops));
      return p;
    }
    ~binaryprops(){
      free(self);
    }
  }
} binaryprops;


%apply float* IN_1D_FLOAT { float *data };
void frotate(float *data, long numbins, float bins_to_left);
/* Rotates a vector by bins_to_left places to the left.    */
/* numbins is the number of FLOATING points to move.       */

%apply double* IN_1D_DOUBLE { double *data };
void drotate(double *data, long numbins, double bins_to_left);
/* Rotates a vector by bins_to_left places to the left.    */
/* numbins is the number of DOUBLE points to move.         */


/*
 *   For all the following:
 *   E           = Eccentric anomaly (radians)
 *   Eo          = initial condition of E
 *   Eacc        = accuracy (1e-14 should be plenty good)
 *   t           = Time (sec)
 *   to          = orbital integration start time
 *   tend        = orbital integration ending time
 *   dt          = sampling interval for the integrator
 *   f_orb       = orbital frequency in hertz
 *   p_orb       = orbital period in seconds
 *   p_psr       = pulsar period in units of choice
 *   e or e_orb  = orbital eccentricity
 *   x or x_orb  = projected semi-major axis of orbit in lt-sec
 *   w or w_orb  = longitude of periastron (degrees)
 */

%apply long ARRAYLEN { long numpts };
double *dorbint(double Eo, long numpts, double dt, orbitparams *orb);
/* This routine integrates Keplar's Equation and returns a double      */
/* vector of the eccentric anomalys (E) for each point.  The initial   */
/* value for eccentric anomaly (usually determined by using            */
/* keplars_equation()) goes in Eo.  The time increment to use is dt,   */
/* total number of pts goes in 'numpts' and all of the various orbital */
/* parameters are found in *orb.  The routine uses 4th order Runge-    */
/* Kutta in a dumb mode (no adaptive step-size) since all we want is   */
/* tabulated results with even intervals.                              */


double keplars_eqn(double t, double p_orb, double e, double Eacc);
/* This routine solves Keplar's Equation at a single time t (sec) and  */
/* returns the value of the eccentric anomaly.  The orbital period (s) */
/* is in p_orb and the orbital eccentricity is in e.  Eacc is the      */
/* absolute accuracy in E that we want to achieve.  t is the time in   */
/* seconds since the last periapsis.  Uses Newton-Raphson.             */

%apply double* IN_1D_DOUBLE { double *E };
double lin_interp_E(double *E, double time, double to,
		    double dt, double maxt);
/* Return a linearly interpolated value of E at time 'time'.  */
/* to is the starting time and dt is the time interval of the */
/* evenly tabulated dat vector *E.                            */

%apply double* IN_1D_DOUBLE { double *E };
void E_to_phib(double *E, long numpoints, orbitparams *orb);
/* Convert eccentric anomalies returned by dorbint into */
/* Phase delays. */

%apply double* IN_1D_DOUBLE { double *E };
void E_to_v(double *E, long numpoints, orbitparams *orb);
/* Convert eccentric anomalies returned by dorbint into */
/* Pulsar line-of-sight velocity (km/s). */

%apply double* IN_1D_DOUBLE { double *E };
void E_to_p(double *E, long numpoints, double p_psr, orbitparams *orb);
/* Convert eccentric anomalies returned by dorbint into */
/* Pulse period. */

%apply double* IN_1D_DOUBLE { double *E };
void E_to_z(double *E, long numpoints, double p_psr, double T,
	    orbitparams *orb);
/* Convert eccentric anomalies returned by dorbint into */
/* Fourier f-dot. */

%apply double* IN_1D_DOUBLE { double *E };
void E_to_phib_BT(double *E, long numpoints, orbitparams *orb);
/* Convert eccentric anomalys (*E) to time delays */
/* using Blanford and Teukolsky Equations         */
/* This model is NOT currently in use.            */

int r_resp_halfwidth(presto_interp_acc accuracy);
  /*  Return the approximate kernel half width in FFT bins required    */
  /*  to achieve a fairly high accuracy correlation based correction   */
  /*  or interpolation for a standard Fourier signal.                  */
  /*  Arguments:                                                       */
  /*    'accuracy' is either LOWACC or HIGHACC.                        */
  /*  Notes:                                                           */
  /*    The result must be multiplied by 2*'numbetween' to get the     */
  /*    length of the array required to hold such a kernel.            */

int z_resp_halfwidth(double z, presto_interp_acc accuracy);
  /*  Return the approximate kernel half width in FFT bins required    */
  /*  to achieve a fairly high accuracy correlation based correction   */
  /*  or interpolation for a Fourier signal with constant f-dot. (i.e  */
  /*  a constant frequency derivative)                                 */
  /*  Arguments:                                                       */
  /*    'z' is the Fourier Frequency derivative (# of bins the signal  */
  /*       smears over during the observation).                        */
  /*    'accuracy' is either LOWACC or HIGHACC.                        */
  /*  Notes:                                                           */
  /*    The result must be multiplied by 2*'numbetween' to get the     */
  /*    length of the array required to hold such a kernel.            */

int w_resp_halfwidth(double z, double w, presto_interp_acc accuracy);
  /*  Return the approximate kernel half width in FFT bins required    */
  /*  to achieve a fairly high accuracy correlation based correction   */
  /*  or interpolation for a Fourier signal with an f-dot that (i.e    */
  /*  varies linearly in time -- a constant f-dotdot)                  */
  /*  Arguments:                                                       */
  /*    'z' is the average Fourier Frequency derivative (# of bins     */
  /*       the signal smears over during the observation).             */
  /*    'w' is the Fourier Frequency 2nd derivative (change in the     */
  /*       Fourier f-dot during the observation).                      */
  /*    'accuracy' is either LOWACC or HIGHACC.                        */
  /*  Notes:                                                           */
  /*    The result must be multiplied by 2*'numbetween' to get the     */
  /*    length of the array required to hold such a kernel.            */

int bin_resp_halfwidth(double ppsr, double T, orbitparams * orbit);
  /*  Return the approximate kernel half width in FFT bins required    */
  /*  to achieve a fairly high accuracy correlation based correction   */
  /*  or interpolation for a pulsar in a binary orbit.                 */
  /*  Arguments:                                                       */
  /*    'ppsr' is the period of the pusar in seconds.                  */
  /*    'T' is the length of the observation in seconds.               */
  /*    'orbit' is a ptr to a orbitparams structure containing the     */
  /*       Keplarian orbital parameters of the binary system.          */
  /*  Notes:                                                           */
  /*    The result must be multiplied by 2 * 'numbetween' to get the   */
  /*    length of the array required to hold such a kernel.            */

%apply int ARRAYLEN { int numkern };
fcomplex *gen_r_response(double roffset, int numbetween, int numkern);
  /*  Generate a complex response function for Fourier interpolation.  */
  /*  Arguments:                                                       */
  /*    'roffset' is the offset in Fourier bins for the full response  */
  /*       (i.e. At this point, the response would equal 1.0)          */
  /*    'numbetween' is the number of points to interpolate between    */
  /*       each standard FFT bin.  (i.e. 'numbetween' = 1 = interbins) */
  /*    'numkern' is the number of complex points that the kernel will */
  /*       contain.                                                    */

%apply int ARRAYLEN { int numkern };
fcomplex *gen_z_response(double roffset, int numbetween, double z,
			 int numkern);
  /*  Generate the response function for Fourier f-dot interpolation.  */
  /*  Arguments:                                                       */
  /*    'roffset' is the offset in Fourier bins for the full response  */
  /*       (i.e. At this point, the response would equal 1.0)          */
  /*    'numbetween' is the number of points to interpolate between    */
  /*       each standard FFT bin.  (i.e. 'numbetween' = 1 = interbins) */
  /*    'z' is the Fourier Frequency derivative (# of bins the signal  */
  /*       smears over during the observation).                        */
  /*    'numkern' is the number of complex points that the kernel will */
  /*       contain.                                                    */

%apply int ARRAYLEN { int numkern };
fcomplex *gen_w_response(double roffset, int numbetween, double z,
			 double w, int numkern);
  /*  Generate the response function for Fourier f-dot interpolation.  */
  /*  Arguments:                                                       */
  /*    'roffset' is the offset in Fourier bins for the full response  */
  /*       (i.e. At this point, the response would equal 1.0)          */
  /*    'numbetween' is the number of points to interpolate between    */
  /*       each standard FFT bin.  (i.e. 'numbetween' = 1 = interbins) */
  /*    'z' is the average Fourier Frequency derivative (# of bins     */
  /*       the signal smears over during the observation).             */
  /*    'w' is the Fourier Frequency 2nd derivative (change in the     */
  /*       Fourier f-dot during the observation).                      */
  /*    'numkern' is the number of complex points that the kernel will */
  /*       contain.                                                    */

%apply int ARRAYLEN { int numkern };
fcomplex *gen_bin_response(double roffset, int numbetween, double ppsr,
			   double T, orbitparams * orbit, int numkern);
  /*  Generate the Fourier response function for a sinusoidal PSR      */
  /*  signal from a binary orbit.                                      */
  /*  Arguments:                                                       */
  /*    'roffset' is the offset in Fourier bins for the full response  */
  /*       (i.e. At this point, the response would equal 1.0)          */
  /*    'numbetween' is the number of points to interpolate between    */
  /*       each standard FFT bin.  (i.e. 'numbetween' = 1 = interbins) */
  /*    'ppsr' is the period of the pusar in seconds.                  */
  /*    'T' is the length of the observation in seconds.               */
  /*    'orbit' is a ptr to a orbitparams structure containing the     */
  /*       Keplarian orbital parameters of the binary system.          */
  /*    'numkern' is the number of complex points that the kernel will */
  /*       contain.                                                    */

float get_numphotons(FILE * file);
  /* Return the total number of photons in the FFT file      */
  /* i.e.  it returns the value of the 0th frequency bin.    */
  /* Arguments:                                              */
  /*   'file' is a pointer to the file you want to access.   */

%apply fcomplex* IN_1D_CFLOAT { fcomplex *data };
float get_localpower(fcomplex *data, int numdata, double r);
  /* Return the local power level at specific FFT frequency.  */
  /* Arguments:                                               */
  /*   'data' is a pointer to a complex FFT.                  */
  /*   'numdata' is the number of complex points in 'data'.   */
  /*   'r' is the Fourier frequency in data that we want to   */
  /*      interpolate.                                        */

%apply fcomplex* IN_1D_CFLOAT { fcomplex *data };
float get_localpower3d(fcomplex *data, int numdata, double r,
		       double z, double w);
  /* Return the local power level around a specific FFT           */
  /* frequency, f-dot, and f-dotdot.                              */
  /* Arguments:                                                   */
  /*   'data' is a pointer to a complex FFT.                      */
  /*   'numdata' is the number of complex points in 'data'.       */
  /*   'r' is the Fourier frequency in data that we want to       */
  /*      interpolate.                                            */
  /*   'z' is the Fourier Frequency derivative (# of bins the     */
  /*       signal smears over during the observation).            */
  /*   'w' is the Fourier Frequency 2nd derivative (change in the */
  /*       Fourier f-dot during the observation).                 */

%apply fcomplex* IN_1D_CFLOAT { fcomplex *data };
void get_derivs3d(fcomplex *data, int numdata, double r,
                  double z, double w, float localpower,
                  rderivs *result);
  /* Return an rderives structure that contains the power,      */
  /* phase, and their first and second derivatives at a point   */
  /* in the F/F-dot/F-dortdot volume.                           */  
  /* Arguments:                                                 */
  /*   'data' is a pointer to a complex FFT.                    */
  /*   'numdata' is the number of complex points in 'data'.     */
  /*   'r' is the Fourier frequency in data that we want to     */
  /*      interpolate.                                          */
  /*   'z' is the Fourier Frequency derivative (# of bins the   */
  /*       signal smears over during the observation).          */
  /*   'w' is the Fourier Frequency 2nd derivative (change in   */
  /*       the Fourier f-dot during the observation).           */
  /*   'localpower' is the local power level around the signal. */
  /*   'result' is a pointer to an rderivs structure that will  */
  /*       contain the results.                                 */

void calc_props(rderivs data, double r, double z, double w, 
		fourierprops * result);
  /* Return a fourierprops structure that contains the various  */
  /* properties of a signal described by Middleditch, Deich,    */ 
  /* and Kulkarni in _Isolated_Pulsars_, 1993, p372.            */  
  /* Arguments:                                                 */
  /*   'data' is a pointer to an rderivs structure containing   */
  /*       derivative information about the peak in question.   */
  /*   'r' is the Fourier frequency in data that we want to     */
  /*      interpolate.                                          */
  /*   'z' is the Fourier Frequency derivative (# of bins the   */
  /*       signal smears over during the observation).          */
  /*   'w' is the Fourier Frequency second derivative.          */
  /*   'result' is a pointer to an fourierprops structure that  */
  /*       will contain the results.                            */

void calc_binprops(fourierprops * props, double T, int lowbin, 
		   int nfftbins, binaryprops * result);
  /* Return a binaryprops structure that contains the various     */
  /* estimates of the binary pulsar system from a mini-FFT.       */
  /* Arguments:                                                   */
  /*   'props' is a pointer to the candidate's fourierprops.      */
  /*   'T' is the total length (sec) of the original time series. */
  /*   'lowbin' is the Fourier bin number from the original FFT   */
  /*      the lowest bin in the mini-FFT.                         */
  /*   'nfftbins' is the number of bins in the mini-FFT.          */
  /*   'absnorm' is the value of the power normalization          */
  /*      constant for this mini-FFT.                             */
  /*   'result' is the returned binaryprops structure.            */

void calc_rzwerrs(fourierprops *props, double T, rzwerrs *result);
  /* Calculate periods, frequencies, their derivatives        */
  /* and their errors.                                        */
  /* Arguments:                                               */
  /*   'props' is a pointer to a fourierprops structure.      */
  /*   'T' is the length of the data set in sec (i.e. N*dt).  */
  /*   'result' is a pointer to the returned rzwerrs struct.  */

void nice_output_1(char *output, double val, double err, int len);
/* Generates a string in "output" of length len with "val" rounded   */
/*   to the appropriate decimal place and the error in parenthesis  */
/*   as in scientific journals.  The error has 1 decimal place.  */

void nice_output_2(char *output, double val, double err, int len);
/* Generates a string in "output" of length len with "val" rounded   */
/*   to the appropriate decimal place and the error in parenthesis  */
/*   as in scientific journals.  The error has 2 decimal places.  */

void print_candidate(fourierprops * cand, double dt, unsigned long N, 
		     float nph, int numerrdigits);
/* Outputs a 2 column summary of all the properties or a fourier peak  */

void print_bin_candidate(binaryprops * cand, int numerrdigits);
/* Outputs a 2 column summary of all the properties or a fourier peak  */

void get_rzw_cand(char *filenm, int candnum, fourierprops *cand);
/*  Read the rzw candidate file 'filenm' and return a        */
/*  pointer to the fourierprops that describes it.           */

void get_bin_cand(char *filenm, int candnum, binaryprops *cand);
/*  Read the bin candidate file 'filenm' and return a        */
/*  pointer to the binaryprops that describes it.            */

unsigned long chkfilelen(FILE *file, int size);
/* Return the length of a file (in blocks of 'size').  */

%apply int ARRAYLEN { int numpts };
fcomplex *read_fcomplex_file(FILE *file, int firstpt, int numpts);
/* Return an fcomplex vector with complex data taken from a file. */
/* Argumants:                                                     */
/*   'file' is a pointer to the file you want to access.          */
/*   'firstpt' is the number of the first point to get. (0 = 1st  */
/*       point in the file).  If < 0, the resulting array will    */
/*       be zero padded.                                          */
/*   'numpts' is the number of points to get from the file.       */
/*       If the number of bins to read takes us past the end of   */
/*       file, the returned vector will be zero padded.           */

%apply float* IN_1D_FLOAT { float *arr };
int prune_powers(float *arr, long n, long numsumpows);
/* Sets powers that are more than PRUNELEV standard devs above */
/* the median value to the median value.  Thereby "cleaning"   */
/* the spectrum of high power signals.                         */

%apply float* IN_1D_FLOAT { float *arr };
float selectkth(long k, long n, float *arr);
/* Selects the kth largest value from the array arr */

double dms2rad(int deg, int min, double sec);
/* Convert degrees, minutes and seconds into radians */

double hms2rad(int hour, int min, double sec);
/* Convert hours, minutes and seconds into radians */

double sphere_ang_diff(double ra1, double dec1, double ra2, double dec2);
/* Return the spherical angle (radians) between two RA and DECS */

%apply fcomplex* IN_1D_CFLOAT { fcomplex *data, fcomplex *result };
void spread_with_pad(fcomplex *data, int numdata, 
		     fcomplex *result, int numresult, 
		     int numbetween, int numpad);
  /* Prepare the data array for correlation by spreading         */
  /*      the input data array and padding it.                   */
  /* Arguments:                                                  */
  /*   'data' is the FFT array to be prepared                    */
  /*   'numdata' is the number of complex points in 'data'       */
  /*   'result' is the prepped data array                        */
  /*   'numresult' is the number of complex points in 'result'   */
  /*   'numbetween' is the number of interpolated pts per bin    */
  /*   'numpad' is the number of bins to use as zero padding     */

%apply fcomplex* IN_1D_CFLOAT { fcomplex *data, fcomplex *result };
void spread_no_pad(fcomplex *data, int numdata, 
		   fcomplex *result, int numresult, 
		   int numbetween);
  /* Prepare the data array for correlation by spreading         */
  /*      the input data array.                                  */
  /* Arguments:                                                  */
  /*   'data' is the FFT array to be prepared                    */
  /*   'numdata' is the number of complex points in 'data'       */
  /*   'result' is the prepped data array                        */
  /*   'numresult' is the number of complex points in 'result'   */
  /*   'numbetween' is the number of interpolated pts per bin    */


%apply fcomplex* IN_1D_CFLOAT { fcomplex *data };
void paddata(fcomplex *data, int numdata, int numpad);
  /* Pad the last 'numpad' bins of 'data' with zeros.         */
  /* Arguments:                                               */
  /*   'data' is the FFT array to be padded                   */
  /*   'numdata' is the number of complex points in 'data'    */
  /*   'numpad' is the number of bins to use as zero padding  */

%apply fcomplex* IN_1D_CFLOAT { fcomplex *kernel, fcomplex *result };
void place_complex_kernel(fcomplex *kernel, int numkernel, 
			  fcomplex *result, int numresult);
  /* This routine places the kernel in a zero filled array */
  /* with half of the response at the beginning and half   */
  /* of the response at the end of the result array.  See  */
  /* Numerical Recipes in C 2ed, p 541 for more info.      */
  /* Arguments:                                            */
  /*   'kernel' is a complex response function.  Bin zero  */
  /*      response is in bin numkernel/2.                  */
  /*   'numkernel' is the number of points in the kernel.  */
  /*      This should be an even number.                   */
  /*   'result' is the result array.                       */
  /*   'numresult' is the number of points in the result.  */

%apply float* IN_1D_FLOAT { float *kernel, float *result };
void place_real_kernel(float *kernel, int numkernel, 
		       float *result, int numresult);
  /* This routine places the kernel in a zero filled array */
  /* with half of the response at the beginning and half   */
  /* of the response at the end of the result array.  See  */
  /* Numerical Recipes in C 2ed, p 541 for more info.      */
  /* Arguments:                                            */
  /*   'kernel' is a real-valued response function.  Bin   */
  /*      zero response is in bin numkernel/2.             */
  /*   'numkernel' is the number of points in the kernel.  */
  /*      This should be an even number.                   */
  /*   'result' is the result array.                       */
  /*   'numresult' is the number of points in the result.  */

%apply fcomplex* IN_1D_CFLOAT { fcomplex *data, fcomplex *result };
void chop_complex_ends(fcomplex *data, int numdata, 
		       fcomplex *result, int numresult, 
		       int chopbins);
  /* Chop the contaminated ends off of an array that has  */
  /* been correlated/convolved.                           */
  /* Arguments:                                           */
  /*   'data' is the array to chop.                       */
  /*   'numdata' is the number of points in data.         */
  /*   'result' is the resultant array.                   */
  /*   'numresult' is the number of points in the result. */
  /*   'chopbins' is the number of bins to chop on each   */
  /*      end of the data array.                          */

%apply float* IN_1D_FLOAT { float *data, float *result };
void chop_real_ends(float *data, int numdata, 
		    float *result, int numresult, 
		    int chopbins);
  /* Chop the contaminated ends off of an array that has  */
  /* been correlated/convolved.                           */
  /* Arguments:                                           */
  /*   'data' is the array to chop.                       */
  /*   'numdata' is the number of points in data.         */
  /*   'result' is the resultant array.                   */
  /*   'numresult' is the number of points in the result. */
  /*   'chopbins' is the number of bins to chop on each   */
  /*      end of the data array.                          */


%apply fcomplex* IN_1D_CFLOAT { fcomplex *data, fcomplex *kernel };
%apply long ARRAYLEN { long numdata };
fcomplex *complex_corr_conv(fcomplex *data, fcomplex *kernel, 
			    int numdata, 
			    presto_ffts ffts, presto_optype type);
  /* Perform and return a complex correlation or convolution.       */
  /* Arguments:                                                     */
  /*   'data' is the complex array to correlate/convolve.           */
  /*   'kernel' is the correlation/convolution kernel.              */
  /*   'numdata' is the length of 'data', 'kernel' and the result.  */
  /*   'ffts' describes how to perform the convolution/correlation. */
  /*      'ffts' = FFTDK:  FFT both the 'data' and the 'kernel'.    */
  /*      'ffts' = FFTD:  FFT only the 'data' not the 'kernel'.     */
  /*      'ffts' = FFTK:  FFT only the 'kernel' not the 'data'.     */
  /*      'ffts' = NOFFTS:  Don't FFT the 'data' or the 'kernel'.   */
  /*   'type' is the type of operation to perform.                  */
  /*      'type' = CONV:  return a convolution in a new vector.     */
  /*      'type' = CORR:  return a correlation in a new vector.     */
  /*      'type' = INPLACE_CONV:  convolution over-writes 'data'.   */
  /*      'type' = INPLACE_CORR:  correlation over-writes 'data'.   */

%apply float* IN_1D_FLOAT { float *data, float *kernel };
%apply long ARRAYLEN { long numdata };
float *real_corr_conv(float *data, float *kernel, int numdata, 
		      presto_ffts ffts, presto_optype type);
  /* Perform and return a real-valued correlation or convolution.   */
  /* Arguments:                                                     */
  /*   'data' is the complex array to correlate/convolve.           */
  /*   'kernel' is the correlation/convolution kernel.              */
  /*   'numdata' is the length of 'data', 'kernel' and the result.  */
  /*   'ffts' describes how to perform the convolution/correlation. */
  /*      'ffts' = FFTDK:  FFT both the 'data' and the 'kernel'.    */
  /*      'ffts' = FFTD:  FFT only the 'data' not the 'kernel'.     */
  /*      'ffts' = FFTK:  FFT only the 'kernel' not the 'data'.     */
  /*      'ffts' = NOFFTS:  Don't FFT the 'data' or the 'kernel'.   */
  /*   'type' is the type of operation to perform.                  */
  /*      'type' = CONV:  return a convolution in a new vector.     */
  /*      'type' = CORR:  return a correlation in a new vector.     */
  /*      'type' = INPLACE_CONV:  convolution over-writes 'data'.   */
  /*      'type' = INPLACE_CORR:  correlation over-writes 'data'.   */

%apply fcomplex* IN_1D_CFLOAT { fcomplex *data, fcomplex *kernel, fcomplex *result };
int corr_complex(fcomplex *data, int numdata, presto_datainf datainf, 
		 fcomplex *kern, int numkern, presto_datainf kerninf, 
		 fcomplex *result, int numresult, int lobin, 
		 int numbetween, int kern_half_width, presto_optype optype);
  /* This routine is a general correlation or convolution routine    */
  /* for complex data.  It can perform convolutions or correlations  */
  /* on raw complex data, data that is prepared for a convolution/   */
  /* correlation but not FFTd, or already FFTd data.  The kernel     */
  /* that it uses can also be raw, prepped, or FFTd.  If you call    */
  /* the routine multiple times with either the same kernel or data  */
  /* array, it uses a saved version of the array from the previous   */
  /* call to cut down on many processing steps. The return value     */
  /* tells how many usable (i.e.  non-contaminated) points were      */
  /* returned in the result array (the first value will be that of   */
  /* 'lobin').  This routine will _not_ perform in-place             */
  /* correlations or convolutions (i.e. it ignores those choices     */
  /* for 'optype').                                                  */
  /* Arguments:                                                      */
  /*   'data' is a complex array of the data to be interpolated.     */
  /*   'numdata' is the number of complex points in 'data'.          */
  /*   'datainf' is one of the following that describes the data:    */
  /*              RAW = Normal un-altered complex data.              */
  /*              PREPPED = Data has been padded and spread based    */
  /*                        on 'kern_half_width' and 'numbetween'    */
  /*                        and is ready to be FFTd.                 */
  /*              FFT = Data has already been prepared and FFTd.     */
  /*              SAME = Data is the same as the previous call.      */
  /*                        The routine uses its saved data.         */
  /*   'kern' is the correlation kernel.                             */
  /*   'numkern' is the number of complex points in 'kern'.          */
  /*   'kerninf' is one of the same choices as 'datainf' above.      */
  /*   'result' is the resulting complex array (must already exist). */
  /*   'numresult' is the number of complex points in 'result'.      */
  /*   'lobin' is the lowest fourier bin to convolve/correlate.      */
  /*   'numbetween' is the number of bins to spread the data points. */
  /*   'kern_half_width' is half the width (bins) of the raw kernel. */
  /*   'optype' is either CORR or CONV (correlation or convolution). */
  /* Notes:                                                          */
  /*   If either 'datainf' or 'kerninf' are of type PREPPED or FFT,  */
  /*   then the length of the FFTs used in the correlation/          */
  /*   convolution calculations will be of length 'numdata' or       */
  /*   'numkern'.  If both 'datainf' and 'kerninf' are of type       */
  /*   PREPPED or FFT then 'numdata' and 'numkern' must have the     */
  /*   same value.  In order for SAME values of 'datainf' and        */
  /*   'kerninf' to help out, the routine must be called with the    */
  /*   same values for 'kern_half_width' and 'numbetween' as well.   */

%apply fcomplex* IN_1D_CFLOAT { fcomplex *data, fcomplex *result };
void stretch_fft(fcomplex *data, int numdata, 
		 fcomplex *result, int numresult);
  /* This routine stretches and/or interpolates an FFT of length    */
  /* numdata.  It zeros 'result' where end-effects have occurred.   */
  /* This routine is usually used to co-add stretched FFTs to       */
  /* increase the signal-to-noise ratios of a detection.            */
  /* Arguments:                                                     */
  /*   'data' is a pointer to a complex FFT.                        */
  /*   'numdata' is the number of complex points in 'data'.         */
  /*   'result' is a pointer to the complex stretched array.        */
  /*   'numresult' is the number of complex points in 'result'.     */
  /* Notes:                                                         */
  /*   The ratio of 'numresult' to 'numdata' determines the amount  */
  /*   of stretching that will take place.  For example, if         */
  /*   'numresult' is twice 'numdata', then the data will be        */
  /*   stretched by a factor of two (i.e. interbinned).             */

%apply float* IN_1D_FLOAT { float *powers };
float *corr_loc_pow(float *powers, int numpowers);
  /* This routine determines the local power levels for every         */
  /* frequency in an FFT containing 'numpowers' complex frequencies.  */
  /* It sets the areas where end-effects are a problem to the         */
  /* local power level of the closest bin without end effect          */
  /* problems.  It returns a vector with the local power levels.      */
  /* Arguments:                                                       */
  /*   'powers' is a pointer to a fcomplex vector containing the FFT. */
  /*   'numpowers' is the number of complex points in 'powers'.       */


%apply fcomplex* IN_1D_CFLOAT { fcomplex *data };
%apply fcomplex *OUTPUT { fcomplex *ans }; 
void rz_interp(fcomplex *data, int numdata, double r, double z,
	       int kern_half_width, fcomplex *ans);
  /* This routine uses the correlation method to do a Fourier        */
  /* complex interpolation at a single point in the f-fdot plane.    */
  /* It does the correlations manually. (i.e. no FFTs)               */
  /* Arguments:                                                      */
  /*   'data' is a complex array of the data to be interpolated.     */
  /*   'numdata' is the number of complex points (bins) in data.     */
  /*   'r' is the Fourier frequency in data that we want to          */
  /*      interpolate.  This can (and should) be fractional.         */
  /*   'z' is the fdot to use (z=f-dot/T^2 (T is integration time)). */
  /*   'kern_half_width' is the half-width of the kernel in bins.    */
  /*   'ans' is the complex answer.                                  */

%apply fcomplex* IN_1D_CFLOAT { fcomplex *data };
%apply double *OUTPUT { double *rout };
double max_r_arr(fcomplex *data, int numdata, double rin, 
		 double *rout, rderivs *derivs);
/* Return the Fourier frequency that maximizes the power.  */


%apply fcomplex* IN_1D_CFLOAT { fcomplex *data };
%apply double *OUTPUT { double *rout, double *zout };
double max_rz_arr(fcomplex *data, int numdata, double rin, double zin, 
		  double *rout, double *zout, rderivs * derivs);
/* Return the Fourier frequency and Fourier f-dot that      */ 
/* maximizes the power.                                     */

%apply double* IN_1D_DOUBLE { double *prof, 
			      double *delays, 
			      double *onoffpairs};
%apply float* IN_1D_FLOAT { float *data, float *chiarr };
%apply long *OUTPUT { long *totnumfolded };
%apply double *OUTPUT { double *avg, double *var };
void foldfile(FILE *datafile, double dt, double *prof, long proflen, 
	      double fo, double fdot, double fdotdot, int binary, 
	      double *delays, double orbto, double orbdt, long numdelays, 
	      double *avg, double *var, float *chiarr, 
	      double *onoffpairs, long *totnumfolded);

%apply double* IN_1D_DOUBLE { double *prof, 
			      double *delays, 
			      double *onoffpairs};
%apply float* IN_1D_FLOAT { float *data };
%apply long *OUTPUT { long *totnumfolded };
void fold(float *data, long N, double dt, double tb, double *prof, 
	  long proflen, double fo, double fdot, double fdotdot, 
	  int binary, double *delays, double orbto, double orbdt, 
	  long numdelays, double *onoffpairs, long *totnumfolded);

double delay_from_dm(double dm, double freq_emitted);
/* Return the delay in seconds caused by dispersion, given  */
/* a Dispersion Measure (dm) in cm-3 pc, and the emitted    */
/* frequency (freq_emitted) of the pulsar in MHz.           */


double dm_from_delay(double delay, double freq_emitted);
/* Return the Dispersion Measure in cm-3 pc, that would     */
/* cause a pulse emitted at frequency 'freq_emitted' to be  */
/* delayed by 'delay' seconds.                              */


double doppler(double freq_observed, double voverc);
/* This routine returns the frequency emitted by a pulsar */
/* (in MHz) given that we observe the pulsar at frequency */
/* freq_observed (MHz) while moving with radial velocity  */
/* (in units of v/c) of voverc wrt the pulsar.            */


%apply fcomplex* IN_1D_CFLOAT { fcomplex *minifft };
%apply float* IN_1D_FLOAT { float *highpows, float *highfreqs };
void search_minifft(fcomplex *minifft, int numminifft,
		    float norm, int numcands, float *highpows,
		    float *highfreqs);
  /* This routine searches a short FFT (usually produced using the   */
  /* MiniFFT binary search method) and returns two vectors which     */
  /* contain the highest powers found and their Fourier frequencies. */
  /* The routine uses interbinning to help find the highest peaks.   */
  /* Arguments:                                                      */
  /*   'minifft' is the FFT to search (complex valued)               */
  /*   'numminifft' is the number of complex points in 'minifft'     */
  /*   'norm' is the value to multiply each pow power by to get      */
  /*      a normalized power spectrum.                               */
  /*   'numcands' is the length of the returned vectors.             */
  /*   'highpows' a vector containing the 'numcands' highest powers. */
  /*   'highfreqs' a vector containing the 'numcands' frequencies    */
  /*      where 'highpows' were found.                               */
  /* Notes:  The returned vectors must have already been allocated.  */
  /*   The returned vectors will be sorted by decreasing power.      */


%apply double* IN_1D_DOUBLE { double *topotimes, 
			      double *barytimes, 
			      double *voverc };
void barycenter(double *topotimes, double *barytimes, 
		double *voverc, long N, char *ra, char *dec,
		char *obs, char *ephem);
/* This routine uses TEMPO to correct a vector of           */
/* topocentric times (in *topotimes) to barycentric times   */
/* (in *barytimes) assuming an infinite observation         */
/* frequency.  The routine also returns values for the      */
/* radial velocity of the observation site (in units of     */
/* v/c) at the barycentric times.  All three vectors must   */
/* be initialized prior to calling.  The vector length for  */
/* all the vectors is 'N' points.  The RA and DEC (J2000)   */
/* of the observed object are passed as strings in the      */
/* following format: "hh:mm:ss.ssss" for RA and             */
/* "dd:mm:s.ssss" for DEC.  The observatory site is passed  */
/* as a 2 letter ITOA code.  This observatory code must be  */
/* found in obsys.dat (in the TEMPO paths).  The ephemeris  */
/* is either "DE200" or "DE400".                            */



/* The FAST FFTs we will use                                     */
/* Note:  A standard forward transform (with negative exponent)  */
/*        is called by setting isign=-1                          */

%apply fcomplex* IN_1D_CFLOAT { fcomplex *indata };
void fftwcall(fcomplex *indata, long nn, int isign);

%apply fcomplex* IN_1D_CFLOAT { fcomplex *indata };
void tablesixstepfft(fcomplex *indata, long nn, int isign);

%apply float* IN_1D_FLOAT { float *data };
%apply long ARRAYLEN { long n };
void realfft(float *data, unsigned long n, int isign);

//
//Note:  The following must be the last declarations in the file
//

%typemap(python, out) fcomplex ** {
  PyArrayObject *arr;
  int n[2];

  n[0] = _output_matrixcols;
  n[1] = _output_matrixrows;
  _output_matrixrows = 0;
  _output_matrixcols = 0;
  arr = (PyArrayObject *) \
    PyArray_FromDimsAndData(2, n, PyArray_CFLOAT, (char *)$source[0]);
  free($source);
  if (arr == NULL) return NULL;
  arr->dimensions[1] = ((int) (*_arg9) - _arg3) * _arg2 ;
  arr->strides[0] = arr->dimensions[1] * sizeof(fcomplex);
  arr->flags |= OWN_DATA;
  PyArray_INCREF(arr);
  $target = (PyObject *)arr;
}

%apply fcomplex* IN_1D_CFLOAT { fcomplex *data };
%apply int MATRIXROWS { int fftlen };
%apply int MATRIXCOLS { int numz };
%apply int *OUTPUT { int *nextbin };
fcomplex **corr_rz_plane(fcomplex *data, int numdata, int numbetween,
			 int startbin, double zlo, double zhi,
			 int numz, int fftlen,
			 presto_interp_acc accuracy, int *nextbin);
  /* This routine uses the correlation method to do Fourier          */
  /* complex interpolations of the f-fdot plane.                     */
  /* Arguments:                                                      */
  /*   'data' is a complex array of the data to be interpolated.     */
  /*   'numdata' is the number of complex points (bins) in data.     */
  /*   'numbetween' is the number of points to interpolate per bin.  */
  /*   'startbin' is the first bin to use in data for interpolation. */
  /*   'zlo' is the lowest fdot to use (z=f-dot/T^2)                 */
  /*   'zhi' is the highest fdot to use (z=f-dot/T^2)                */
  /*   'numz' is the number of z values to use to make the plane     */
  /*   'fftlen' is the # of complex pts in kernel and result.        */
  /*   'accuracy' is either HIGHACC or LOWACC.                       */
  /*   'nextbin' will contain the bin number of the first bin not    */
  /*      interpolated in data.                                      */

%typemap(python, out) fcomplex * {
  PyArrayObject *arr;
  int n;
  
  n = _output_arraylen;
  _output_arraylen = 0;
  arr = (PyArrayObject *) \
    PyArray_FromDimsAndData(1, (int *)&n, PyArray_CFLOAT, (char *)$source);
  if (arr == NULL) return NULL;
  arr->flags |= OWN_DATA;
  arr->dimensions[0] = ((int) (*_arg7) - _arg3) * _arg2;
  arr->strides[0] = arr->dimensions[1] * sizeof(fcomplex);
  PyArray_INCREF(arr);
  $target = (PyObject *)arr;
}

%apply fcomplex* IN_1D_CFLOAT { fcomplex *data };
%apply int ARRAYLEN { int fftlen };
%apply int *OUTPUT { int *nextbin };
fcomplex *corr_rz_interp(fcomplex *data, int numdata, int numbetween,
			 int startbin, double z, int fftlen,
			 presto_interp_acc accuracy, int *nextbin);
  /* This routine uses the correlation method to do a Fourier        */
  /* complex interpolation of a slice of the f-fdot plane.           */
  /* Arguments:                                                      */
  /*   'data' is a complex array of the data to be interpolated.     */
  /*   'numdata' is the number of complex points (bins) in data.     */
  /*   'numbetween' is the number of points to interpolate per bin.  */
  /*   'startbin' is the first bin to use in data for interpolation. */
  /*   'z' is the fdot to use (z=f-dot/T^2).                         */
  /*   'fftlen' is the # of complex pts in kernel and result.        */
  /*   'accuracy' is either HIGHACC or LOWACC.                       */
  /*   'nextbin' will contain the bin number of the first bin not    */
  /*      interpolated in data.                                      */

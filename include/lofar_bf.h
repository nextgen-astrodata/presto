/* based on sigproc_fb.h  and psrfits.h */

/*  LOFAR data is represented in a hierarchical format, 
    stored in HDF5 format; for ease of internal representation
    of this hierarchy we define structs for each of the subgroups
    SUB_ARRAY_POINTING, BEAM, STOKESS which then can be used in LOFARBF
                                                                  */  
#define MAXPFITSFILES 1000

/* TODO: This file contains DAL specific code, so it must be compiled/handled with g++ */
//#include "/opt/cep/dal/current/include/dal/lofar/BF_File.h" /* DAL BF_File class */

typedef struct{
  /* SUB_ARRAY_POINTING group header attributes */
  char expTimeStartUTC[40]; /* Start time of obs (UTC time) */
  char expTimeEndUTC[40];   /* End time of obs (UTC time) */
  double expTimeStartMJD; /* Start time of obs (MJD time) */
  double expTimeEndMJD;   /* End time of obs (MJD time) */
  double totalIntegrationTime;  /* Total integration time of this SAP */
  double pointRA;         /* J2000 RA of SAP at observation start */
  double pointDEC;        /* J2000 DEC of SAP at observation start */
  int observationNofBeams;  /* Number of beams in the entire ob- servation */
  int nofBeams;           /* Number of beams stored in this file */
  int ibeam;             /* Beam number used for this data */
}sap_info;

typedef struct{
  /* BEAM group header attributes */  
  char targets[120];  /* List of observation targets/- sources observed within this Beam */
  int nofStations;    /* Number of stations used within this Beam */
  //char stationsList[300]; /* List of stations used for this Beam */
  unsigned long  nofSamples;  /* List of stations used for this Beam */
  double samplingRate;      /* List of stations used for this Beam */
  double samplingTime;      /* Sampling time */
  int channelsPerSubband;   /* Number of channels for each sub- band */
  double subbandWidth;      /* Subband width, typically 156250 or 195312.5 */
  double channelWidth;      /* */
  double tracking;        /* Tracking method */
  /* J2000 right ascension of the cen- ter of the beam at start of obser- vation (at LOFAR core)*/
  double pointRA;
  /* J2000 declination of the cen- ter of the beam at start of obser- vation (at LOFAR core)*/
  double pointDEC;
  /* RA offset of beam from sub- array pointing at start of obser- vation */
  double pointOffsetRA;
  /* Dec offset of beam from sub- array pointing at start of obser- vation */  
  double pointOffsetDEC;
  /* The diameter of the beam el- lipse in the major axis a (RA), at zenith at center frequency */
  double beamDiameterRA;
  char beamDiameterRAUnit[10];
  /* The diameter of the beam ellipse in the minor axis b (Dec), at zenith at center frequency */
  double beamDiameterDEC; 
  /* Beam center frequency - the mid- dle frequency of all the channel frequencies */
  double beamFrequencyCenter;
  char beamFrequencyCenterUnit[10]; 
 /* Are the data folded within this Beam? */
  int foldedData;
  /* The fold period if the data were folded */
  double foldPeriod;
  char foldPeriodUnit[10];
  /* Where the data dedispersed in- coherently	(‘INCOHERENT’), dedispersed	coherently 
    (‘COHERENT’), or not dedis- persed at all (‘NONE’) */
  double dedispersion;
  /* The dispersion measure applied to the data, if data were de- dispersed */
  double dedispersionMeasure;
  char dedispersionMeasureUnit[10];
}beam_info;

/* Coordinates */
typedef struct{
  double *refLocationValue;   /* Numerical value(s) of the reference location */
  char *refLocationUnit[10];  /* Physical unit(s) for the reference location */
  /* Identifier for the reference system of the location; see Tab. 10 for a list of recognized values. */
  char refLocationFrame[20];
  double refTimeValue;        /* Numerical value of the reference time */
  char refTimeUnit[10];       /* Physical unit of the reference time */
  char refTimeFrame[20];      /* Identifier for the reference time system used */
  int nofCoordinates;         /* Identifier for the reference time system used */
  int nofAxes;                /* Identifier for the reference time system used */
  char *coordinateTypes[20];  /* embedded coordinate object types */
}coordinates_info;

typedef struct{
  char storageType[10];   /* Coordinate storage type */
  int nofAxes;            /* Number of coordinate axes */
  char **axisNames;       /* World axis names */
  char **axisUnits;       /* World axis units */
  double referenceValue;  /* Reference value */
  double referencePixel;  /* Reference pixel */
  double increment;       /* Coordinate increment */
  /* The World Coordinate Reference scaling delta matrix */
// double *PC;
//  double *axisValuesPixel;  /* Reference pixels*/
//  double *axisValuesWorld;  /* Reference values */
}timeCoordinate_info;

typedef struct{
  char coordinateType[20];  /* Coordinate Type descriptor: 'Spectral' */
  char storageType[20];     /* Coordinate storage type: 'Tabular' */
  int nofAxes;              /* Number of coordinate axes */
  //char **axisNames;       /* World axis name */
  //char **axisUnits;       /* World axis units */
  double referenceValue;    /* Reference value (CRVAL) */
  double referencePixel;    /* Reference value (CRPIX) */
  double increment;         /* Coordinate increment (CDELT) */
  //double *PC;               /* The World Coordinate Reference scaling delta matrix */
  //double *axisValuesPixel;  /* Reference pixels */
  //double *axisValuesWorld;  /* Reference values */
}spectralCoordinate_info;

/* Stokes Dataet */
typedef struct{
  char dataType[10];          /* Data type used */
  char stokesComponent[10];   /* Stokes component stored within the dataset */
  int nofSamples;             /* Number of bins along the time axis of the dataset */
 /* Number of frequency sub-bands, into which the frequency domain is being divided */
  int nofSubbands;            
  int *nofChannels;       /* Number of channels (frequency bins) within each subband */
}stokes_info;

typedef struct LOFARBF {
  char filename[80];      /* LOFAR BFwriter h5 filename */
  char rawFilename[80];   /* LOFAR BFwriter raw filename */
  char fileDate[80];      /* File creation date */
  char fileType[80];      /* File type (see Table 5) in ICD003*/
  char telescope[80];     /* Name of the telescope */

  int telescope_id;      /* Telescope ID (see telescope_name() */

  char projectID[40];     /* Project identifier */
  char projectTitle[80];  /* Title	of	the	project */
  char projectPI[80];     /* Name of Principal Investigator */
  char projectCOI[160];     /* Name(s) of the Coinvestigator(s) */
  char projectContact[80];  /* Project contact	details */
  char observationID[20];   /* Observation unique identifier */
  char observationStartUTC[20]; /* Observation start date (UTC)*/
  char observationEndUTC[20];   /* Observation end date (UTC) */
  //float observationStartMJD[20];
  //float observationEndMJD[20];
  int observationNofStations;         /* nof. stations used during the ob- servation */
//  char observationStationList;      /* char string with comma separated list of stations */
  double observationFrequencyMin;     /* Observation minimum frequency */
  double observationFrequencyCenter;  /* Observation center frequency (weighted) */
  double observationFrequencyMax;     /* Observation maximum fre- quency*/
  /* Number of bits per sample in the incoming data stream from the stations to CEP/BlueGene */
  int observationNofBitsPerSample;
  double clockFrequency;            /* Clock frequency, valid values for LOFAR are 160.0 MHz and 200.0 MHz */
  char antennaSet[20];              /* Antenna set specification of observation */
  char filterSelection[20]; /* e.g. ‘LBA_10_70' */
  char targets[80];         /* observed targets  */
  char systemVersion[40];   /* Processing system name/version */
  char pipelineName[40];    /* Pipeline processing name */
  char pipelineVersion[40]; /* Pipeline processing version */
  char docName[40];       /* Interface Control Document name */
  char docVersion[20];    /* Interface Control Document version/issue number */
  //char notes[80]        /* Notes or comments */
  /* Whether the file was created ‘Online’ (stream to file) or ‘Of- fline’ (file to file) */
  char createOfflineOnline[10]; 
  char bfFormat[3];       /* " RAW" if the file is BeamFormed raw or "TAB" processed data */
  char bfVersion[20];     /* BeamFormed data format version number */
  double totalIntegrationTime;  /* Total integration time of the ob- servation */
  char observationDatatype[10]; /* Type of observation: Searching, timing, etc. */
  double subArrayPointingDiameter;  /* FWHM of the sub-array pointing at zenith at center frequency */
  double bandwidth;       /* Total bandwidth (excluding gaps) */
  /* Number of sub-array pointing synthesized beams in the entire observation */
  int observationNofSubArrayPointings;  /* Number of sub-array pointing synthesized beams in the entire observation */
  int nofSubArrayPointings; /* Number of sub-array pointing synthesized beams stored in this file */
  
  /* SUB_ARRAY_POINTING group header attributes */
  sap_info sap;

  /* BEAM group header attributes */  
  beam_info beam;

  int sumifs;            /* Whether the IFs are summed or not */
  int headerlen;         /* Header length in bytes */
} lofarbf_info;


/* For the time being try to reuse spectra_info from psrfits.h */
/*
struct spectra_info {
    char telescope[40];     // Telescope used
    char observer[40];      // Observer's name
    char source[40];        // Source name
    char frontend[40];      // Frontend used
    char backend[40];       // Backend or instrument used
    char project_id[40];    // Project identifier
    char date_obs[40];      // Start of observation (YYYY-MM-DDTHH:MM:SS.SSS)
    char ra_str[40];        // Right Ascension string (HH:MM:SS.SSSS)
    char dec_str[40];       // Declination string (DD:MM:SS.SSSS)
    char poln_type[40];     // Polarization recorded (LIN or CIRC)
    char poln_order[40];    // Order of polarizations (i.e. XXYYXYYX)
    long long N;            // Total number of spectra in observation
    double T;               // Total observation time (N * dt)
    double dt;              // Sample duration (s)
    double fctr;            // Center frequency of the observing band (MHz)
    double lo_freq;         // Center freq of the lowest channel (MHz)
    double hi_freq;         // Center freq of the highest channel (MHz)
    double orig_df;         // Original frequency spacing between the channels (MHz)
    double chan_dm;         // DM used to de-disperse the freq channels
    double df;              // Frequency spacing between the channels (MHz)
    double BW;              // Bandwidth of the observing band (MHz)
    double ra2000;          // RA  of observation (deg, J2000) at obs start
    double dec2000;         // Dec of observation (deg, J2000) at obs start
    double azimuth;         // Azimuth (commanded) at the start of the obs (deg)
    double zenith_ang;      // Zenith angle (commanded) at the start of the obs (deg)
    double beam_FWHM;       // Beam FWHM (deg)
    double time_per_subint; // Duration (in sec) of a full SUBINT entry
    int scan_number;        // Number of scan
    int tracking;           // Tracking (1) or drift scan (0)
    int orig_num_chan;      // Number of spectral channels per sample
    int num_channels;       // Number of spectral channels per sample
    int num_polns;          // 1, 2, or 4 (for Full Stokes)
    int summed_polns;       // Are polarizations summed? (1=Yes, 0=No)
    int FITS_typecode;      // FITS data typecode as per CFITSIO
    int bits_per_sample;    // Number of bits per sample
    int bytes_per_spectra;  // Number of bytes in a spectra (inluding all polns)
    int samples_per_spectra;// Number of samples in a spectra (inluding all polns)
    int bytes_per_subint;   // Number of bytes per SUBINT entry
    int spectra_per_subint; // Number of spectra per SUBINT entry
    int samples_per_subint; // Number of samples per SUBINT entry
    int num_files;          // Number of files in the observation
    int offs_sub_col;       // The number of the OFFS_SUB column in the SUBINT HDU
    int dat_wts_col;        // The number of the DAT_WTS column in the SUBINT HDU
    int dat_offs_col;       // The number of the DAT_OFFS column in the SUBINT HDU
    int dat_scl_col;        // The number of the DAT_SCL column in the SUBINT HDU
    int data_col;           // The number of the DATA column in the SUBINT HDU
    int apply_scale;        // Do we apply the scales to the data? (1=Yes, 0=No)
    int apply_offset;       // Do we apply the offsets to the data? (1=Yes, 0=No)
    int apply_weight;       // Do we apply the weights to the data? (1=Yes, 0=No)
    int apply_flipband;     // Do we invert the band? (1=Yes, 0=No)
    int flip_bytes;         // Hack to flip the order of the bits in a byte of data
    float zero_offset;      // A DC zero-offset value to apply to all the data
    float clip_sigma;       // Clipping value in standard deviations to use
    long double start_MJD[MAXPFITSFILES]; // Array of long double MJDs for the file starts
//    fitsfile *files[MAXPFITSFILES];       // Array of num_files FITS file pointers
    dal::BF_File files[MAXPFITSFILES];    // Array of num_files BF_File objects    
    int start_subint[MAXPFITSFILES];      // Starting ISUBINT in each file
    int num_subint[MAXPFITSFILES];        // Number of subints per file
    long long start_spec[MAXPFITSFILES];  // Starting spectra for each file (zero offset)
    long long num_spec[MAXPFITSFILES];    // Number of spectra per file
    long long num_pad[MAXPFITSFILES];     // Number of padding samples after each file
};
*/

int is_LOFARBF(const char *filename);

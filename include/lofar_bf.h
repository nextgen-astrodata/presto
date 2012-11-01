/* based on sigproc_fb.h  and psrfits.h */

/*  LOFAR data is represented in a hierarchical format, 
    stored in HDF5 format; for ease of internal representation
    of this hierarchy we define structs for each of the subgroups
    SUB_ARRAY_POINTING, BEAM, STOKESS which then can be used in LOFARBF
                                                                  */  
#define MAXPFITSFILES 1000

/* TODO: This file contains DAL specific code, so it must be compiled/handled with g++ */
//#include "/opt/cep/dal/current/include/dal/lofar/BF_File.h" /* DAL BF_File class */

/*
// For the time being try to reuse spectra_info from psrfits.h 
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
/* Functions according to PSRFITS file format handling */
//void read_LOFARBF_files(const char **filenames, int numfiles, struct spectra_info *s);
void read_LOFARBF_files(struct spectra_info *s);
int is_LOFARBF(const char *filename);
void print_LOFARBF_info(struct spectra_info *s);

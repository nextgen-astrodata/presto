/*
  Read information from LOFAR BFwriter HDF5 files  

  Author:       Sven Duscha (duscha@astron.nl)
  Date:         2012-09-13
  Last change:  2012-11-07
*/

#include "/home/duscha/DAL/include/dal/lofar/CLA_File.h" 	// DAL includes
#include "/home/duscha/DAL/include/dal/lofar/BF_File.h" 	// DAL includes
#include <iostream>
#include <vector>
#include <numeric>

// PRESTO C includes, with our own lofar_bf.h function declarations
extern "C"
{
	#include "presto.h"
	#include "mask.h"
	#include "backend_common.h"
	#include "misc_utils.h" 
	#include "lofar_bf.h"
}

using namespace std;
using namespace dal;

#ifndef CHARLEN
#define CHARLEN 40
#endif

#ifndef MAXLEN
#define MAXLEN 500
#endif

// Helper functions
std::string doubleToString(double num);
std::string intToString(int num);
std::string longToString(long num);
std::vector<int> findSaps(dal::BF_File &file);
std::vector<int> findBeams(dal::BF_SubArrayPointing &sap);


// Read and convert PSRFITS information from a group of files 
// and place the resulting info into a spectra_info structure.
//void read_LOFARBF_files(const char **filenames, int numfiles, struct spectra_info *s)
void read_LOFARBF_files(struct spectra_info *s)
{
    // TODO: make SAP, BEAM and STOKES selectable through function parameters
    int sapNr=0;
    int beamNr=0;
    int stokesNr=0;

		/*
		s->fitsfiles = (fitsfile **)malloc(sizeof(fitsfile *) * s->num_files);
		*/
    s->start_subint = gen_ivect(s->num_files);
    s->num_subint = gen_ivect(s->num_files);
    s->start_spec = (long long *)malloc(sizeof(long long) * s->num_files);
    s->num_spec = (long long *)malloc(sizeof(long long) * s->num_files);
    s->num_pad = (long long *)malloc(sizeof(long long) * s->num_files);
    s->start_MJD = (long double *)malloc(sizeof(long double) * s->num_files);
    s->N = 0;
    //s->num_beams = 1;

		/* TODO: Is this for LOFAR BF data relevant? */
    // By default, don't flip the band.  But don't change
    // the input value if it is aleady set to flip the band always
    //if (s->apply_flipband==-1) s->apply_flipband = 0;
    s->apply_flipband=0;    // LOFAR BF don't flip the band
    s->files=(FILE **)malloc(sizeof(FILE *) * s->num_files);

    // Step through the other files
    for (int ii = 0; ii < s->num_files; ii++)
		{
//#if DEBUG_OUT       
        printf("Reading '%s'\n", s->filenames[ii]);
//#endif
        // Is the file a LOFARBF file?
        if (!is_LOFARBF(s->filenames[ii])) {
            fprintf(stderr, 
                    "\nError!  File '%s' does not appear to be LOFARBF!\n", 
                    s->filenames[ii]);
            exit(1);
        }
        
				BF_File file(s->filenames[ii], BF_File::READ);			// Open the LOFARBF file
        s->files[ii]=fopen(s->filenames[ii], "r");          // get "r"  file pointer to file for PRESTO
 
				// Read all spectra_info equivalent information from LOFARBF file through DAL2.5
				strncpy(s->telescope, file.telescope().get().c_str() , CHARLEN);
				strncpy(s->observer, file.projectPI().get().c_str() , CHARLEN);

				BF_SubArrayPointing sap(file.subArrayPointing(findSaps(file)[0]));	// TODO: allow other SAPs
				BF_BeamGroup beam(sap.beam(findBeams(sap)[0]));				              // TODO: allow multiple Beams

				// Targets				
				vector<string> targets(beam.targets().get());	    // TODO: write all targets into ->source?
				strncpy(s->source, targets[0].c_str(), CHARLEN);	// spectra_info leaves only CHARLEN chars for (one) target 

				// Frontend and Backend
				// CLOCK_FREQUENCY, use this as "Frontend" entry, use FILTER_SELECTION as "Backend"
				if(file.clockFrequency().exists())
				{
	        if(file.clockFrequencyUnit().exists())
	        {
	          string backend=doubleToString(file.clockFrequency().value).append(" ").append(file.clockFrequencyUnit().value);
	          strncpy(s->backend, backend.c_str(), CHARLEN); 
	        }
	        else
	        {
		  		strncpy(s->backend, file.filterSelection().get().c_str(), CHARLEN);
				  }
				}
				if(file.antennaSet().exists())
				{
					strncpy(s->frontend, file.antennaSet().get().c_str(), CHARLEN);
				}
				// PROJECT_ID
				strncpy(s->project_id, file.projectID().get().c_str() , CHARLEN);

				// Date of Observation UTC
				strncpy(s->date_obs, file.observationStartUTC().get().c_str(),CHARLEN);

				// MJD Start time
				s->start_MJD=(long double*) malloc(sizeof(long double));
				s->start_MJD[0]=file.observationStartMJD().value;

				//strncpy(s->ra_str, doubleToString(sap.pointRA().value).c_str(), CHARLEN);
				//strncpy(s->ra_str, doubleToString(sap.pointDEC().value).c_str(), CHARLEN);

				// RA J2000 and Dec J2000 in degrees
				s->ra2000=sap.pointRA().value;
				s->dec2000=sap.pointDEC().value;
				// RA J2000 and Dec J2000 in HH:MM:SS
				int hh=0, mm=0;  // hour and minute for RA and DEC conversion
				double ss=0;     // seconds
			  deg2RAhms(s->ra2000, &hh, &mm, &ss);
			  char RA[40], DEC[40];
        ra_dec_to_string(RA, hh, mm, ss);
        strncpy(s->ra_str, RA, CHARLEN);
        deg2dms(s->dec2000, &hh, &mm, &ss);
        //string DEC=intToString(hh).append(":").append(intToString(mm)).append(":").append(intToString(ss));      				
        ra_dec_to_string(DEC, hh, mm, ss);
        strncpy(s->dec_str, DEC, CHARLEN);
				// Tracking
				if(beam.tracking().exists())
				{
					if(beam.tracking().get()=="J2000")				
						s->tracking=1; 
					else
						s->tracking=0;
				}
				// Azimuth and altitude, both are only optional
				/*
				if(beam.pointAzimuth().exists())
				{
					s->azimuth=beam.pointAzimuth().get();
				}	
				if(beam.pointAltitude().exists())
				{
					s->zenith_ang=90-beam.pointAltitude().get(); // Zenith_ang = 90 - altitude
				}	
        */
				s->azimuth=0;      // as long as these are not in the file, set them to -1
        s->zenith_ang=0;

				// Polarizations Type and order, Number of spectra, TODO: only I supported for now
				// Read polarization from the BEAM group, these are the polarizations in THIS file
    	 	strncpy(s->poln_type, "LIN", CHARLEN);   // Polarization recorded (LIN or CIRC), for LOFAR BF always "LIN"
        string stokesComponents;
        vector<string> stokesComponentsVec=beam.stokesComponents().get();
        for(unsigned int i=0; i<stokesComponentsVec.size(); i++)
        {
          stokesComponents.append(stokesComponentsVec[i]);
        }
        strncpy(s->poln_order, stokesComponents.c_str(), CHARLEN);
        s->num_polns=stokesComponents.size();
        /* Hack? if summed_polns is set to true, then s->num_pols is taken as 2, thought it should be 1 */
        s->summed_polns=0;            // summed polarizations are always summed for LOFAR, i.e. I= XX + YY

        // Number of beams and beam selection
        s->num_beams=sap.nofBeams().value;
        s->beamnum=0;                         // selected beam (TODO: allow beam selection)

        // Sample time, use unit to convert to "us" LOFAR BF default = "s"
        if(beam.samplingTimeUnit().get() == "s")
          s->dt=beam.samplingTime().value;        
        else if(beam.samplingTimeUnit().get() == "ms")        
          s->dt=beam.samplingTime().value*10e-3;
        else if(beam.samplingTimeUnit().get() == "us")        
          s->dt=beam.samplingTime().value*10e-6;        
        else if(beam.samplingTimeUnit().get() == "ns")        
          s->dt=beam.samplingTime().value*10e-9;                
        
        // Frequencies
        // Central frequency, default in LOFAR BF is "MHz"
        if(beam.beamFrequencyCenterUnit().get() == "MHz")
          s->fctr=beam.beamFrequencyCenter().value;
        if(beam.beamFrequencyCenterUnit().get() == "kHz")
          s->fctr=beam.beamFrequencyCenter().value*10e-3; 
        // Low channel / High channel
        if(beam.beamFrequencyCenterUnit().get() == "MHz")
        {
          s->lo_freq=file.observationFrequencyMin().value;
          s->hi_freq=file.observationFrequencyMax().value;
        }
        else if(beam.beamFrequencyCenterUnit().get() == "kHz")
        {
          s->lo_freq=file.observationFrequencyMin().value*10e-3; 
          s->hi_freq=file.observationFrequencyMax().value*10e-3; 
        }
        // Channel width (LOFAR BF default="Hz")
        if(beam.channelWidthUnit().get() == "kHz")
          s->orig_df=beam.channelWidth().value*10e-3;
        else if(beam.channelWidthUnit().get() == "Hz")
          s->df=beam.channelWidth().value*10e-6;
        // Bandwidth, check for Unit if it is MHz, kHz, or Hz
        if(file.bandwidthUnit().get()=="MHz")
          s->BW=file.bandwidth().value;
		    else if(file.bandwidthUnit().get()=="kHz")
		      s->BW=file.bandwidth().value * 1000;
		    else
		      s->BW=file.bandwidth().value * 1000000;
        // Beam FWHM: BEAM_DIAMETER in RA (default="arcmin")
        if(beam.beamDiameterRAUnit().get() == "deg")
          s->beam_FWHM=beam.beamDiameterRA().value;
        else if(beam.beamDiameterRAUnit().get() == "arcmin")
          s->beam_FWHM=beam.beamDiameterRA().value*60;
        else if(beam.beamDiameterRAUnit().get() == "arcsec")
          s->beam_FWHM=beam.beamDiameterRA().value*60;
        // number of channels
        s->orig_num_chan=beam.channelsPerSubband().value;    // Number of spectral channels per sample
        
        //BF_Stokes stokes=beam.stokes(stokesNr);         // Get STOKES_0 , TODO: make this selectable  
        vector<unsigned int> v=beam.stokes(stokesNr).nofChannels().get();    
        s->num_channels=std::accumulate(v.begin(), v.end(), 0); // use STL to sum over channels per subband
        // Total number of spectra samples in observation
     		s->N=beam.nofSamples().value*stokesComponents.size();

        if(beam.stokes(stokesNr).dataType().get()=="float")
          s->bits_per_sample=32;     // Bytes per sample default=4 (float)
        else if(beam.stokes(stokesNr).dataType().get()=="int")
          s->bits_per_sample=16;
        else if(beam.stokes(stokesNr).dataType().get()=="char")  
          s->bits_per_sample=8;
        
        s->T=sap.totalIntegrationTime().value;                // total time in observation         
        s->time_per_subint=sap.totalIntegrationTime().value;
        s->bytes_per_spectra=s->num_channels * s->bits_per_sample;  // Number of bytes in a spectra (inluding all polns)
        //s->samples_per_spectra=beam.nofSamples().value;  // TODO: is this correct?
        s->samples_per_spectra = s->num_polns * s->num_channels;  // <- this is from psrfits.c implementation
        s->zero_offset=0;         // LOFAR BF data has no zero offset
        s->header_offset=0;       // NO bytes to skip from header
        s->offset_to_spectra=0;   // bytes in file header (NOT needed for LOFAR)

        s->num_spec[ii] = beam.nofSamples().value;           // number of spectra per file
        s->num_pad[ii] = 0;       // there are no padding samples after each file
        
        s->spectra_per_subint=s->num_channels;
        s->samples_per_subint=beam.nofSamples().value;
        long long p=0;
        s->num_pad=&p;     // LOFAR BF H5 does not use padding
        
        //TODO: fill in remaining information into spectra_info
		}
}

/*
long long (*offset_to_LOFARBF_spectra)(long long specnum, struct spectra_info *s)
{
  // Open LOFARBF RAW file with actual data
  
  // set file pointer in s to new file
  
  // return offset=0 (LOFAR BF has only 1 spectrum per file)  
}
*/

// Helper functions to convert numbers to strings
//
// convert a double to a string
string doubleToString(double num)
{
	ostringstream osstream;
	osstream << num;
	return osstream.str();
}

// convert an int to a string
string intToString(int num)
{
	ostringstream osstream;
	osstream << num;
	return osstream.str();
}

// convert a long to a string
string longToString(long num)
{
	ostringstream osstream;
	osstream << num;
	return osstream.str();
}

// Find SAP Nr in file
vector<int> findSaps(BF_File &file)
{
	vector<int> saps;
	for(unsigned int i=0; i<file.nofSubArrayPointings().value; i++)
	{
		if(file.subArrayPointing(i).exists())
		{
			saps.push_back(i);
		}
	}

	return saps;
}

// Find BEAM Nr in SAP
vector<int> findBeams(BF_SubArrayPointing &sap)
{
	vector<int> beams;
	for(unsigned int i=0; i<sap.observationNofBeams().value; i++)
	{
		if(sap.beam(i).exists())
		{
			beams.push_back(i);
		}
	}

	return beams;
}

// Check if the file is a LOFAR BeamFormed HDF5 file 
int is_LOFARBF(const char *filename)
{
	string name(filename);
	if(name.find("_bf.h5"))	// Check file name convention first (see ICD005)
	{
		try								  // try to open file as a BF_File
		{		
			BF_File file(filename);
		}		
		catch(DALException e)
		{
			cout << "Error: " << filename << " could not be opened as a LOFAR BF_File." << endl;
			return(0);		// not a LOFAR BF_File
		}
	}
	else
	{
		return(0);			// not a LOFAR BF_File
	}
	return(1);				// filename is indeed a LOFAR BF_File
}

void get_LOFAR_BF_file_info(const char *filename)
{


}

void read_LOFAR_BF(FILE * infiles[], int numfiles, float *data,
              int numpts, double *dispdelays, int *padding,
              int *maskchans, int *nummasked, mask * obsmask)
{

}

void print_LOFARBF_info(struct spectra_info *s)
{
	printf("print_LOFARBF_info() s->num_files: %d\n", s->num_files);	/* DEBUG */
	//printf("From the LOFARBF HDF5 file '%s':\n", s->files[0]->_ptr->filename);
	// Loop over the files provided in spectra_info s
	int i=0;
	for(i=0; i<s->num_files; i++)
	{
 		cout << "Information for" << s->filenames[i] << endl;
    
	}
}

// Read a rawblock from a LOFAR BF raw file
int get_LOFARBF_rawblock(float *fdata, struct spectra_info *s, int *padding)
{
  return 0;
}

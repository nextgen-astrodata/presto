/*
  Read information from LOFAR BFwriter HDF5 files  

  Author:       Sven Duscha (duscha@astron.nl)
  Date:         2012-09-13
  Last change:  2012-11-01
*/

#include "/home/duscha/DAL/include/dal/lofar/CLA_File.h" 	// DAL includes
#include "/home/duscha/DAL/include/dal/lofar/BF_File.h" 	// DAL includes
#include <iostream>
#include <vector>

// PRESTO C includes, with our own lofar_bf.h function declarations
extern "C"
{
	#include "presto.h"
	#include "mask.h"
	#include "backend_common.h"
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

// Read and convert PSRFITS information from a group of files 
// and place the resulting info into a spectra_info structure.
//void read_LOFARBF_files(const char **filenames, int numfiles, struct spectra_info *s)
void read_LOFARBF_files(struct spectra_info *s)
{
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
    s->num_beams = 1;

		/* TODO: Is this for LOFAR BF data relevant? */
    // By default, don't flip the band.  But don't change
    // the input value if it is aleady set to flip the band always
    if (s->apply_flipband==-1) s->apply_flipband = 0;

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
        
				BF_File file(s->filenames[ii], BF_File::READ);						// Open the LOFARBF file

				// Read all spectra_info equivalent information from LOFARBF file through DAL2.5
				strncpy(s->telescope, file.telescope().get().c_str() , CHARLEN);
				strncpy(s->observer, file.projectPI().get().c_str() , CHARLEN);
        // TODO: when DAL copy constructor bug has been fixed, use functions again
				string sapNr=intToString(findSaps(file)[0]);  // DEBUG
				cout << "sapNr: " << sapNr << endl; 	// DEBUG	
				BF_SubArrayPointing sap(file.subArrayPointing(findSaps(file)[0]));	  // TODO: allow other SAPs
				string beamNr=intToString(findBeams(sap)[0]);
				BF_BeamGroup beam(sap.beam(findBeams(sap)[0]));				// TODO: allow multiple Beams
				cout << "beamNr: " << beamNr << endl; // DEBUG

/*
        // Get SAP Numbers from file
    	  vector<int> saps;
	      for(unsigned int i=0; i<file.nofSubArrayPointings().value; i++)
	      {
	        cout << "i: " << i << endl;
	      	if(file.subArrayPointing(i).exists())
		      {
			      saps.push_back(i);
		      }
	      }
        BF_SubArrayPointing sap(file, file.subArrayPointing(saps[0]).name());
        // Get BEAM numbers from file  
	      vector<int> beams;
        for(unsigned int i=0; i<sap.nofBeams().value; i++)
	      {
	        cout << "sap.nofBeams().value: " << sap.nofBeams().value << endl;   // DEBUG
		      cout << "sap.beam(i).exists(): " << sap.beam(i).exists() << endl;   // DEBUG
		      if(sap.beam(i).exists())
		      {
		        cout << "beams.push_back("<< i << ")" << endl;  // DEBUG
			      beams.push_back(i);
		      }
	      }
	      cout << "beams[0]: " << beams[0] << endl;
        if(beams.size()==0)
        {
          fprintf(stderr, "error: %s SAP %d does not contain any beams", s->filenames[ii], saps[0]);
          exit(1);
        }
        BF_BeamGroup beam(sap, sap.beam(beams[0]).name());
*/
				// Targets				
				vector<string> targets(beam.targets().get());	
				strncpy(s->source, targets[0].c_str(), CHARLEN);	// spectra_info leaves only CHARLEN chars for (one) target 

				// Frontend and Backend
				// SAMPLING_RATE, use this as "Frontend" entry, leave "Backend" empty
				if(file.antennaSet().exists())
				{
					strncpy(s->frontend, file.antennaSet().get().c_str(), CHARLEN);
				}
				if(beam.samplingRate().exists())
				{
		  		strncpy(s->backend, doubleToString(beam.samplingRate().value).c_str(), CHARLEN);
				}
				// PROJECT_ID
				strncpy(s->project_id, file.projectID().get().c_str() , CHARLEN);

				// Date of Observation UTC
				strncpy(s->date_obs, file.observationStartUTC().get().c_str(),CHARLEN);

				// MJD Start time
				s->start_MJD=(long double*) malloc(sizeof(long double));
				s->start_MJD[0]=file.observationStartMJD().value;

				// RA J2000 and Dec J2000 in HH:MM:SS
				//strncpy(s->ra_str, doubleToString(sap.pointRA().value).c_str(), CHARLEN);
				//strncpy(s->ra_str, doubleToString(sap.pointDEC().value).c_str(), CHARLEN);

				// RA J2000 and Dec J2000 in degrees
				s->ra2000=sap.pointRA().value;
				s->dec2000=sap.pointDEC().value;

				// Tracking
				if(beam.tracking().exists())
				{
					if(beam.tracking().get()=="J2000")				
						s->tracking=1; 
					else
						s->tracking=0;
				}
				// Azimuth and altitude
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

				// Polarizations Type and order, Number of spectra, TODO: only I supported for now
    	 	strncpy(s->poln_type, "LIN", CHARLEN);   // Polarization recorded (LIN or CIRC)
     		strncpy(s->poln_order, "I", CHARLEN);    // Order of polarizations (i.e. XXYYXYYX) TODO: get Stokes parameter

				//BF_Stokes stokes0(beam, "STOKES_0");	// open STOKES_0 group
     		//s->N=stokes0.nofSamples().value;    // Total number of spectra in observation
		}
}


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

void set_LOFAR_padvals(float *fpadvals, int good_padvals)
{  
   int ii;
   float sum_padvals = 0.0;

//   if (good_padvals) {
//      for (ii = 0; ii < numchan_st; ii++) {
//         padvals[ii] = newpadvals[ii] = (unsigned char) (fpadvals[ii] + 0.5);
//         sum_padvals += fpadvals[ii];
//      }
//      padval = (unsigned char) (sum_padvals / numchan_st + 0.5);
//   } else {
//      for (ii = 0; ii < numchan_st; ii++)
//         padvals[ii] = newpadvals[ii] = padval;
//   }

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


/*
  Read information from LOFAR BFwriter HDF5 files  

  Author:       Sven Duscha (duscha@astron.nl)
  Date:         2012-09-13
  Last change:  2012-10-17
*/

#include "/home/duscha/DAL/include/dal/lofar/CLA_File.h" 	// DAL includes
#include "/home/duscha/DAL/include/dal/lofar/BF_File.h" 	// DAL includes
#include <iostream>

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

// Read and convert PSRFITS information from a group of files 
// and place the resulting info into a spectra_info structure.
void read_LOFARBF_files(struct spectra_info *s)
{
/*
    s->datatype = PSRFITS;
    s->fitsfiles = (fitsfile **)malloc(sizeof(fitsfile *) * s->num_files);
    s->start_subint = gen_ivect(s->num_files);
    s->num_subint = gen_ivect(s->num_files);
    s->start_spec = (long long *)malloc(sizeof(long long) * s->num_files);
    s->num_spec = (long long *)malloc(sizeof(long long) * s->num_files);
    s->num_pad = (long long *)malloc(sizeof(long long) * s->num_files);
    s->start_MJD = (long double *)malloc(sizeof(long double) * s->num_files);
    s->N = 0;
    s->num_beams = 1;
*/
}

// Check if the file is a LOFAR BeamFormed HDF5 file 
int is_LOFAR_BF(const char *filename)
{
	string name=filename;
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
	//printf("From the LOFARBF HDF5 file '%s':\n", s->files[0]->_ptr->filename);
	// Loop over the files provided in spectra_info s
 
	int i=0;
	for(i=0; i<s->num_files; i++)
	{
 		cout << "Information for" << s->filenames[i] << endl;
	}
}


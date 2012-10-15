/*
  Read information from LOFAR BFwriter HDF5 files  

  Author:       Sven Duscha (duscha@astron.nl)
  Date:         2012-09-13
  Last change:  2012-10-15
*/

#include "presto.h"
#include "mask.h"
#include "backend_common.h"
#include "lofar_bf.h"

// Read and convert PSRFITS information from a group of files 
// and place the resulting info into a spectra_info structure.
void read_LOFARBF_files(struct spectra_info *s)
{

}

/* Check if the file is a LOFAR BeamFormed HDF5 file */
int is_LOFAR_BF(const char *filename)
{

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
   /*
   int ii;
   float sum_padvals = 0.0;

   if (good_padvals) {
      for (ii = 0; ii < numchan_st; ii++) {
         padvals[ii] = newpadvals[ii] = (unsigned char) (fpadvals[ii] + 0.5);
         sum_padvals += fpadvals[ii];
      }
      padval = (unsigned char) (sum_padvals / numchan_st + 0.5);
   } else {
      for (ii = 0; ii < numchan_st; ii++)
         padvals[ii] = newpadvals[ii] = padval;
   }
  */
}

void print_LOFARBF_info(struct spectra_info *s)
{
	//printf("From the LOFARBF HDF5 file '%s':\n", s->files[0]->_ptr->filename);
	/* Loop over the files provided in spectra_info s */
 
	unsigned int i=0;
	for(i=0; i<s->num_files; i++)
	{
 		printf("Information for %s:\n", s->filenames[i]);
	}
}

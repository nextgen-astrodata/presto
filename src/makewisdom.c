#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "meminfo.h"
#include "sfftw.h"

int main(void)
{
  FILE *wisdomfile;
  fftw_plan plan;
  int fftlen;

  fftlen = 2;

  /* Generate the wisdom... */

  printf("\nCreating Wisdom for FFTW.\n");
  printf("This may take a while...\n\n");
  printf("Generating plans for FFTs of length:\n");

  while (fftlen <= BIGFFTWSIZE) {
    printf("   %d\n", fftlen);
    plan = fftw_create_plan(fftlen, -1, \
	FFTW_MEASURE | FFTW_USE_WISDOM | FFTW_IN_PLACE);
    fftw_destroy_plan(plan);
    plan = fftw_create_plan(fftlen, 1, \
	FFTW_MEASURE | FFTW_USE_WISDOM | FFTW_IN_PLACE);
    fftw_destroy_plan(plan);
    fftlen <<= 1;
  }

  printf("Exporting wisdom to 'fftw_wisdom.txt'\n");

  /* Open wisdom file for writing... */

  wisdomfile = fopen("fftw_wisdom.txt", "w");

  /* Write the wisdom... */

  fftw_export_wisdom_to_file(wisdomfile);

  /* Cleanup... */

  fftw_forget_wisdom();
  fclose(wisdomfile);
  printf("Done.\n\n");

  return (0);

}
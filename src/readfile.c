#include <limits.h>
#include <ctype.h>
#include "presto.h"
#include "mask.h"
#include "multibeam.h"
#include "bpp.h"
#include "wapp_key.h"
#include "spigot.h"
#include "psrfits.h"
#include "lofar_bf.h" 
#include "sigproc_fb.h"
#include "readfile_cmd.h"

/* #define DEBUG */

#define PAGELEN 32              /* Set the page length to 32 lines */
#define NUMTYPES 18             /* number of file types to support */

void print_WAPP_hdr(struct HEADERP *hdr);

int BYTE_print(long count, char *obj_ptr);
int FLOAT_print(long count, char *obj_ptr);
int DOUBLE_print(long count, char *obj_ptr);
int FCPLEX_print(long count, char *obj_ptr);
int DCPLEX_print(long count, char *obj_ptr);
int SHORT_print(long count, char *obj_ptr);
int INT_print(long count, char *obj_ptr);
int LONG_print(long count, char *obj_ptr);
int RZWCAND_print(long count, char *obj_ptr);
int BINCAND_print(long count, char *obj_ptr);
int POSITION_print(long count, char *obj_ptr);
int PKMBHDR_print(long count, char *obj_ptr);
int BCPMHDR_print(long count, char *obj_ptr);
int WAPPHDR_print(long count, char *obj_ptr);
int SPIGOTHDR_print(long count, char *obj_ptr);
int SIGPROCHDR_print(long count, char *obj_ptr);
int PSRFITSHDR_print(long count, char *obj_ptr);
int LOFARHDR_print();
void print_rawbincand(rawbincand cand);
void set_WAPP_HEADER_version(struct HEADERP *hdr);

typedef enum {
    BYTE, FLOAT, 
    DOUBLE, FCPLEX, 
    DCPLEX, SHORT, INT, LONG,
    RZWCAND, BINCAND, POSITION, 
    PKMBHDR, BCPMHDR, WAPPHDR, SPIGOTHDR, SIGPROCHDR, PSRFITSHDR,
    LOFARHDR		/* BF writer HDF5 header of LOFAR */
} rawtypes;

typedef struct fcplex {
   float r;
   float i;
} fcplex;

typedef struct dcplex {
   double r;
   double i;
} dcplex;

int type_sizes[NUMTYPES] = {
   sizeof(unsigned char),
   sizeof(float),
   sizeof(double),
   sizeof(fcplex),
   sizeof(dcplex),
   sizeof(short),
   sizeof(int),
   sizeof(long),
   sizeof(fourierprops),
   sizeof(rawbincand),
   sizeof(position),
   49792, /* This is the length of a Parkes Multibeam record */
   32768, /* This is the length of a BCPM header */
   0,     /* "Special" length for a WAPP header */
   0,     /* "Special" length for a SPIGOT header */
   0,     /* "Special" length for a SIGPROC header */
   0      /* "Special" length for a LOFAR BF header */
};

int objs_at_a_time[NUMTYPES] = {
   PAGELEN, PAGELEN, PAGELEN, PAGELEN, PAGELEN, PAGELEN, PAGELEN,
   PAGELEN, 1, 1, PAGELEN, 1, 1, 1, 1, 1
};

/* You don't see this every day -- An array of pointers to functions: */
int (*print_funct_ptrs[NUMTYPES]) () = {BYTE_print,
                                        FLOAT_print,
                                        DOUBLE_print,
                                        FCPLEX_print,
                                        DCPLEX_print,
                                        SHORT_print,
                                        INT_print,
                                        LONG_print,
                                        RZWCAND_print,
                                        BINCAND_print,
                                        POSITION_print,
                                        PKMBHDR_print,
                                        BCPMHDR_print, 
                                        WAPPHDR_print, 
                                        SPIGOTHDR_print, 
                                        SIGPROCHDR_print,
                                        PSRFITSHDR_print,												       
                                        LOFARHDR_print};

/* A few global variables */

long N;
double dt, nph;

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

int main(int argc, char **argv)
{
   int index = -1, need_type = 0;
   int objs_read, objs_to_read, has_suffix;
   long i, j, ct;
   char *cptr, *data, *short_filenm, *extension, key = '\n';
   FILE *infile;
   Cmdline *cmd;
   infodata inf;

   /* Call usage() if we have no command line arguments */

   if (argc == 1) {
      Program = argv[0];
      usage();
      exit(0);
   }

   /* Parse the command line using the excellent program Clig */

   cmd = parseCmdline(argc, argv);

#ifdef DEBUG
   showOptionValues();
#endif

   //fprintf(stdout, "\n\n  PRESTO Binary File Reader\n");
   //fprintf(stdout, "     by Scott M. Ransom\n\n");

   /* Set our index value */

   if (cmd->bytP || cmd->sbytP)
      index = BYTE;
   else if (cmd->fltP || cmd->sfltP)
      index = FLOAT;
   else if (cmd->dblP || cmd->sdblP)
      index = DOUBLE;
   else if (cmd->fcxP || cmd->sfcxP)
      index = FCPLEX;
   else if (cmd->dcxP || cmd->sdcxP)
      index = DCPLEX;
   else if (cmd->shtP || cmd->sshtP)
      index = SHORT;
   else if (cmd->igrP || cmd->sigrP)
      index = INT;
   else if (cmd->lngP || cmd->slngP)
      index = LONG;
   else if (cmd->rzwP || cmd->srzwP)
      index = RZWCAND;
   else if (cmd->binP || cmd->sbinP)
      index = BINCAND;
   else if (cmd->posP || cmd->sposP)
      index = POSITION;
   else if (cmd->pkmbP)
      index = PKMBHDR;
   else if (cmd->bcpmP)
      index = BCPMHDR;
   else if (cmd->wappP)
      index = WAPPHDR;
   else if (cmd->spigotP)
      index = SPIGOTHDR;
   else if (cmd->filterbankP)
      index = SIGPROCHDR;
   else if (cmd->psrfitsP)
      index = PSRFITSHDR;
   else if (cmd->lofarP)
      index = LOFARHDR;

   /* Try to determine the data type from the file name */

   if (index == -1) {
      has_suffix = split_root_suffix(cmd->argv[0], &short_filenm, &extension);
      if (!has_suffix) {
         need_type = 1;
      } else {
         if (strlen(extension) < 2) { /* LOFAR BF writer uses ".h5" files, so we need to support .h5 */
            need_type = 1;
         } else {
            if (0 == strcmp(extension, "dat")) {
               index = FLOAT;
               fprintf(stdout, "Assuming the data is floating point.\n\n");
            } else if (0 == strcmp(extension, "sdat")) {
               index = SHORT;
               fprintf(stdout, "Assuming the data is short integers.\n\n");
            } else if (0 == strcmp(extension, "fft")) {
               index = FCPLEX;
               fprintf(stdout, "Assuming the data is single precision complex.\n\n");
            } else if (0 == strcmp(extension, "fits")) {
                if (strstr(short_filenm, "spigot_5") != NULL) {
                    cmd->spigotP = 1;
                    index = SPIGOTHDR;
                    fprintf(stdout,
                            "Assuming the data is from the Caltech/NRAO Spigot.\n\n");
                } else if (is_PSRFITS(cmd->argv[0])) {
                    cmd->psrfitsP = 1;
                    index = PSRFITSHDR;
                    fprintf(stdout,
                            "Assuming the data is in PSRFITS format.\n\n");
                }
            } else if (0 == strcmp(extension, "h5")) {
                cmd->lofarP = 1;
                index = LOFARHDR;
                fprintf(stdout, "Assuming the data is in LOFAR BF format.\n\n");
                
            } else if (0 == strcmp(extension, "bcpm1") ||
                       0 == strcmp(extension, "bcpm2")) {
               cmd->bcpmP = 1;
               index = BCPMHDR;
               fprintf(stdout, "Assuming the data is from a BCPM machine.\n\n");
            } else if (0 == strcmp(extension, "pkmb")) {
               cmd->pkmbP = 1;
               index = PKMBHDR;
               fprintf(stdout,
                       "Assuming the data is from the Parkes Multibeam machine.\n\n");
            } else if (0 == strcmp(extension, "fil") || 0 == strcmp(extension, "fb")) {
               cmd->filterbankP = 1;
               index = SIGPROCHDR;
               fprintf(stdout,
                       "Assuming the data is a SIGPROC filterbank file.\n\n");
            } else if (isdigit(extension[0]) &&
                       isdigit(extension[1]) && isdigit(extension[2])) {
               cmd->wappP = 1;
               index = WAPPHDR;
               fprintf(stdout, "Assuming the data is from a WAPP machine.\n\n");
            } else if (0 == strcmp(extension, "pos")) {
               index = POSITION;
               fprintf(stdout,
                       "Assuming the data contains 'position' structures.\n\n");
            } else if (0 == strcmp(extension, "cand")) {
               /* A binary or RZW search file? */
               if (NULL != (cptr = strstr(cmd->argv[0], "_bin"))) {
                  index = BINCAND;
                  fprintf(stdout,
                          "Assuming the file contains binary candidates.\n\n");
               } else if (NULL != (cptr = strstr(cmd->argv[0], "_rzw"))) {
                  index = RZWCAND;
                  ct = (long) (cptr - cmd->argv[0]);
                  fprintf(stdout, "Assuming the file contains 'RZW' candidates.\n");
                  free(short_filenm);
                  short_filenm = (char *) malloc(ct + 1);
                  short_filenm[ct] = '\0';
                  strncpy(short_filenm, cmd->argv[0], ct);
                  fprintf(stdout, "\nAttempting to read '%s.inf'.  ", short_filenm);
                  readinf(&inf, short_filenm);
                  fprintf(stdout, "Successful.\n");
                  N = (long) (inf.N + DBLCORRECT);
                  dt = inf.dt;
                  if (cmd->nphP)
                     nph = cmd->nph;
                  else
                     nph = 1.0;
                  fprintf(stdout,
                          "\nUsing N = %ld, dt = %g, and DC Power = %f\n\n",
                          N, dt, nph);
               } else if (NULL != (cptr = strstr(cmd->argv[0], "_ACCEL"))) {
                  index = RZWCAND;
                  ct = (long) (cptr - cmd->argv[0]);
                  fprintf(stdout, "Assuming the file contains 'RZW' candidates.\n");
                  free(short_filenm);
                  short_filenm = (char *) malloc(ct + 1);
                  short_filenm[ct] = '\0';
                  strncpy(short_filenm, cmd->argv[0], ct);
                  fprintf(stdout, "\nAttempting to read '%s.inf'.  ", short_filenm);
                  readinf(&inf, short_filenm);
                  fprintf(stdout, "Successful.\n");
                  N = (long) (inf.N + DBLCORRECT);
                  dt = inf.dt;
                  if (cmd->nphP)
                     nph = cmd->nph;
                  else
                     nph = 1.0;
                  fprintf(stdout,
                          "\nUsing N = %ld, dt = %g, and DC Power = %f\n\n",
                          N, dt, nph);
		       	} 	else {
										 printf("need type!\n");
	               	   need_type = 1;
							  }
				}				
				else if (0 == strcmp(extension, "h5")) {
					    printf("LOFAR BFwriter HDF5 file.\n");		/* DEBUG */
            } else
              	need_type = 1;
         }
      }

      /* If no file extension or if we don't understand the extension, exit */

      if (need_type) {
         fprintf(stdout, "You must specify a data type for this file.\n\n");
         free(short_filenm);
         exit(-1);
      }
      free(short_filenm);
      if (has_suffix)
         free(extension);
   }

//	 printf("need_type = %d\n", need_type);						/* DEBUG */
//	 printf("cmd->index[1] = %d\n", cmd->index[1]);		/* DEBUG */
//	 printf("cmd->lofarP = %d\n", cmd->lofarP);				/* DEBUG */
//	 printf("cmd->psrfitsP = %d\n", cmd->psrfitsP); 	/* DEBUG */

   if (cmd->index[1] == -1 || cmd->index[1] == 0)
      cmd->index[1] = INT_MAX;
   if (cmd->index[1] < cmd->index[0]) {
      fprintf(stdout, "\nThe high index must be >= the low index.");
      fprintf(stdout, "  Exiting.\n\n");
      exit(-1);
   }

   if (cmd->psrfitsP) {
       struct spectra_info s;
       // -1 causes the data to determine if we use weights, scales, & offsets
/* TODO: Scott has to fix this       
			s.apply_weight = s.apply_scale = s.apply_offset = -1;
       if (ls_PSRFITS_files(cmd->argv, cmd->argc, &s)) {
           print_PSRFITS_info(&s);
           printf("\n");
       } else {
           printf("\n  Error reading PSRFITS file!\n\n");
       }
*/
       exit(0);
   }

  /* LOFAR BFwriter files */
  if (cmd->lofarP) {
		struct spectra_info s;
		s.num_files = cmd->argc;
    s.datatype = LOFARBF;

 		if (s.num_files > MAXPFITSFILES) {
        printf("Error!: There are more than %d input files!\n", MAXPFITSFILES);
        exit(1);
    }

		// copy filenames into spectra_info
		unsigned i=0;		
		s.filenames = (char**) calloc(sizeof(char*), s.num_files);
		for(i=0; i < s.num_files; i++)
		{
			s.filenames[i] = (char*) calloc(sizeof(char*), strlen(argv[i+1]));
			strncpy(s.filenames[i], argv[i+1], 255);
		}

		read_LOFARBF_files(&s);
  	//print_LOFARBF_info(&s);
		print_spectra_info(&s);

    exit(0);
  }

   if (cmd->spigotP) {
      SPIGOT_INFO spigot;
			/*
      if (read_SPIGOT_header(cmd->argv[0], &spigot)) {
         print_SPIGOT_header(&spigot);
         printf("\n");
      } else {
         printf("\n  Error reading spigot file!\n\n");
      }
      exit(0);
   		*/
	 }

   if (cmd->filterbankP) {
      sigprocfb fb;

      infile = chkfopen(cmd->argv[0], "rb");
			/*
      if (read_filterbank_header(&fb, infile)) {
         print_filterbank_header(&fb);
         printf("\n");
      } else {
         printf("\n  Error reading SIGPROC filterbank file!\n\n");
      }
      exit(0);
		 */   
	}

   if (cmd->wappP) {
       struct HEADERP *hdr = NULL;
       /*
       infile = chkfopen(cmd->argv[0], "rb");
       hdr = head_parse(infile);
       set_WAPP_HEADER_version(hdr);
       if (hdr) {
           print_WAPP_hdr(hdr);
           printf("\n");
       } else {
           printf("\n  Error reading WAPP file!\n\n");
       }
       exit(0);
			*/   
	 }

   /* Open the file */

   infile = chkfopen(cmd->argv[0], "rb");

   if (cmd->fortranP) {
      chkfileseek(infile, 1, sizeof(long), SEEK_SET);
   }

   /* Skip to the correct first object */

   if (cmd->index[0] > 0) {
      chkfileseek(infile, (long) (cmd->index[0]), type_sizes[index], SEEK_CUR);
   }

   /* Read the file */

   objs_to_read = objs_at_a_time[index];
   data = (char *) malloc(type_sizes[index] * objs_at_a_time[index]);

   i = cmd->index[0];
   do {
      if (objs_to_read > cmd->index[1] - i)
         objs_to_read = cmd->index[1] - i;
      objs_read = chkfread(data, type_sizes[index], objs_to_read, infile);
      for (j = 0; j < objs_read; j++)
         print_funct_ptrs[index] (i + j, data + j * type_sizes[index]);
      /* Just print 1 header for BCPM and WAPP files */
      if (index == BCPMHDR || index == WAPPHDR || index == SPIGOTHDR)
         break;
      i += objs_read;
      if (cmd->pageP) {
         fflush(NULL);
         fprintf(stdout, "\nPress ENTER for next page, or any other key and ");
         fprintf(stdout, "then ENTER to exit.\n\n");
         key = getchar();
      }
   } while (!feof(infile) && i < cmd->index[1] && key == '\n');

   fflush(NULL);
   if (feof(infile)) {
      fprintf(stdout, "\nEnd of file.\n\n");
   }

   free(data);
   fclose(infile);
   exit(0);
}


int BYTE_print(long count, char *obj_ptr)
{
   unsigned char *object;

   object = (unsigned char *) obj_ptr;
   printf("%9ld:  %d\n", count, *object);
   return 0;
}

int FLOAT_print(long count, char *obj_ptr)
{
   float *object;

   object = (float *) obj_ptr;
   printf("%9ld:  %11.7g\n", count, *object);
   return 0;
}

int DOUBLE_print(long count, char *obj_ptr)
{
   double *object;

   object = (double *) obj_ptr;
   printf("%9ld:  %20.15g\n", count, *object);
   return 0;
}

int FCPLEX_print(long count, char *obj_ptr)
{
   fcplex *object;

   object = (fcplex *) obj_ptr;
   printf("%9ld:  %11.7g + %11.7g i\n", count, object->r, object->i);
   return 0;
}

int DCPLEX_print(long count, char *obj_ptr)
{
   dcplex *object;

   object = (dcplex *) obj_ptr;
   printf("%9ld:  %20.15g + %20.15g i\n", count, object->r, object->i);
   return 0;
}

int SHORT_print(long count, char *obj_ptr)
{
   short *object;

   object = (short *) obj_ptr;
   printf("%9ld:  %6d\n", count, *object);
   return 0;
}

int INT_print(long count, char *obj_ptr)
{
   int *object;

   object = (int *) obj_ptr;
   printf("%9ld:  %10d\n", count, *object);
   return 0;
}

int LONG_print(long count, char *obj_ptr)
{
   long *object;

   object = (long *) obj_ptr;
   printf("%9ld:  %12ld\n", count, *object);
   return 0;
}

int RZWCAND_print(long count, char *obj_ptr)
{
   fourierprops *object;

   object = (fourierprops *) obj_ptr;
   printf("\n%ld:\n", count + 1);
   print_candidate(object, dt, N, nph, 2);
   return 0;
}

int BINCAND_print(long count, char *obj_ptr)
{
   rawbincand *object;

   object = (rawbincand *) obj_ptr;
   printf("\n%ld:\n", count + 1);
   print_rawbincand(*object);
   return 0;
}

int POSITION_print(long count, char *obj_ptr)
{
   position *object;

   object = (position *) obj_ptr;
   printf("%ld:  pow = %-7g    p1 = %-13.2f   p2 = %-9.2f   p3 = %-9.2f\n",
          count + 1, object->pow, object->p1, object->p2, object->p3);
   return 0;
}

int PKMBHDR_print(long count, char *obj_ptr)
{
   PKMB_tapehdr *object;

   object = (PKMB_tapehdr *) obj_ptr;
   printf("\n%ld:", count + 1);
   /*print_PKMB_hdr(object); */
   return 0;
}

int BCPMHDR_print(long count, char *obj_ptr)
{
   BPP_SEARCH_HEADER *object;

   object = (BPP_SEARCH_HEADER *) obj_ptr;
   printf("\n%ld:", count + 1);
   /*print_BPP_hdr(object);*/
   return 0;
}

int WAPPHDR_print(long count, char *obj_ptr)
{
   char *bogus;

   printf("\n%ld:", count + 1);
   bogus = (char *) obj_ptr;
   return 0;
}

int SPIGOTHDR_print(long count, char *obj_ptr)
{
   char *bogus;

   printf("\n%ld:", count + 1);
   bogus = (char *) obj_ptr;
   return 0;
}

int PSRFITSHDR_print(long count, char *obj_ptr)
{
   char *bogus;
   
   printf("\n%ld:", count + 1);
   bogus = (char *) obj_ptr;
   return 0;
}

/*
  LOFARHDR_print uses the LOFAR DAL2.x to open a HDF5 BFwriter file and
  display header information
*/
int LOFARHDR_print()
{
	struct spectra_info *s;
	print_LOFARBF_info(s);

  return 0;
}

int SIGPROCHDR_print(long count, char *obj_ptr)
{
   printf("\n%ld:", count + 1);
   /* print_filterbank_header((sigprocfb *) obj_ptr); */
   return 0;
}

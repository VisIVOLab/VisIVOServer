
/**
 * \file io_tipsy_header.c
 *
 * Provides functions for reading and writing the header of TIPSY
 * files.
 */


/**********************************************************************\
 *    Includes                                                        * 
 \**********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#ifdef HAVE_STDINT_H
#	include <stdint.h>
#else
#	include <replace_stdint.h>
#endif

#include "io_tipsy_header.h"
#include "io_util.h"


/**********************************************************************\
 *    Local defines, structure definitions and typedefs               * 
 \**********************************************************************/
#define SKIP(f) {fseek(f, 4L, SEEK_CUR);}


/**********************************************************************\
 *    Prototypes of local functions                                   * 
 \**********************************************************************/


/**********************************************************************\
 *    Implementation of global functions                              * 
 \**********************************************************************/
extern io_tipsy_header_t
io_tipsy_header_get(io_logging_t log, io_tipsy_t f)
{
  io_tipsy_header_t dummy;
  long skipsize;
  int iloop;
  FILE *ftipsy;
  char tipsyline[2048];
  
  /* Some sanity checks */
  if ((f == NULL) || (f->file == NULL))
    return NULL;
  
  /* Check if there already is a header, do nothing then */
  if (f->header != NULL)
    return f->header;
  
  /* Create the header structure array */
  dummy = (io_tipsy_header_t)malloc((size_t)TIPSY_HEADER_SIZE+TIPSY_HEADER_EXTRA);
  if (dummy == NULL) {
    io_logging_memfatal(log, "TIPSY header structure");
    return NULL;
  }
  
  /* just to be on the safe side: rewind the file to the start... */
  rewind(f->file);
  
  /* Now start reading the header */
  io_util_readdouble(f->file, &(dummy->time), f->swapped);
  io_util_readint32(f->file, &(dummy->nbodies), f->swapped);
  io_util_readint32(f->file, &(dummy->ndim), f->swapped);
  io_util_readint32(f->file, &(dummy->nsph), f->swapped);
  io_util_readint32(f->file, &(dummy->ndark), f->swapped);
  io_util_readint32(f->file, &(dummy->nstar), f->swapped);
  io_util_readint32(f->file, &(dummy->pad), f->swapped);

  
  
  
  /* read TIPSY relevant (unit) information from an extra user provided file */
  if((ftipsy = fopen("tipsy.info","r")) == NULL)
    {
      io_logging_fatal(log,"Could not open tipsy.info containing the following information:");
      io_logging_fatal(log,"0.3                     TIPSY_OMEGA0");
      io_logging_fatal(log,"0.7                     TIPSY_LAMBDA0");
      io_logging_fatal(log,"20.0                    TIPSY_BOXSIZE           [Mpc/h]");
      io_logging_fatal(log,"690.988298942671        TIPSY_VUNIT             [km/sec]");
      io_logging_fatal(log,"2.2197e15               TIPSY_MUNIT             [Msun/h]");
      io_logging_fatal(log,"Please generate this file. Aborting now!");
      free(dummy);
      return NULL;
      //exit(0);
    }
  fgets(tipsyline, 2048, ftipsy);
  sscanf(tipsyline,"%lf",&(dummy->TIPSY_OMEGA0));
  fgets(tipsyline, 2048, ftipsy);
  sscanf(tipsyline,"%lf",&(dummy->TIPSY_LAMBDA0));
  fgets(tipsyline, 2048, ftipsy);
  sscanf(tipsyline,"%lf",&(dummy->TIPSY_BOXSIZE));
  fgets(tipsyline, 2048, ftipsy);
  sscanf(tipsyline,"%lf",&(dummy->TIPSY_VUNIT));
  fgets(tipsyline, 2048, ftipsy);
  sscanf(tipsyline,"%lf",&(dummy->TIPSY_MUNIT));  
  fclose(ftipsy);
  
  f->header = dummy;

  return dummy;
}

extern void
io_tipsy_header_del(io_logging_t log, io_tipsy_header_t *header)
{
  if ( (header == NULL) || (*header == NULL) )
    return;
  
  free(*header);
  
  *header = NULL;
  
  return;
}

extern void
io_tipsy_header_write(io_logging_t log,
                      io_tipsy_header_t header,
                      io_tipsy_t f)
{
  if ( (header == NULL) || (f == NULL))
    return;
  
  if (f->mode != IO_FILE_WRITE)
    return;
  
  if (f->file == NULL)
    return;
  
  if (header != f->header) {
    io_logging_msg(log, INT32_C(1),
                   "Writing a different header than stored in "
                   "the file object to the file.");
  }
  
  /* TODO: Write the header */
  
  return;
}


extern void
io_tipsy_header_log(io_logging_t log, io_tipsy_header_t header)
{
  io_logging_msg(log, INT32_C(5),
                 "   time:                         %e",
                 header->time);
  io_logging_msg(log, INT32_C(5),
                 "  nbodies:                       %" PRIi32,
                 header->nbodies);
  io_logging_msg(log, INT32_C(5),
                 "  ndim:                          %" PRIi32,
                 header->ndim);
  io_logging_msg(log, INT32_C(5),
                 "  nsph:                          %" PRIi32,
                 header->nsph);
  io_logging_msg(log, INT32_C(5),
                 "  ndark:                         %" PRIi32,
                 header->ndark);
  io_logging_msg(log, INT32_C(5),
                 "  nstar:                         %" PRIi32,
                 header->nstar);  
  io_logging_msg(log, INT32_C(5),
                 "  pad:                           %" PRIi32,
                 header->pad);  
  return;
}


/**********************************************************************\
 *    Implementation of local functions                               * 
 \**********************************************************************/

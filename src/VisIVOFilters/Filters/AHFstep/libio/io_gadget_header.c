/* $Id: io_gadget_header.c,v 1.3 2007/02/05 16:18:01 knolli Exp $ */

/**
 * \file io_gadget_header.c
 *
 * Provides functions for reading and writing the header of Gadget
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

#include "io_gadget_header.h"
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
extern io_gadget_header_t
io_gadget_header_get(io_logging_t log, io_gadget_t f)
{
	io_gadget_header_t dummy;
	long skipsize;

	/* Some sanity checks */
	if ((f == NULL) || (f->file == NULL))
		return NULL;

	/* Check if there already is a header, do nothing then */
	if (f->header != NULL)
		return f->header;

	/* Create the header structure array */
	dummy = (io_gadget_header_t)malloc((size_t)GADGET_HEADER_SIZE);
	if (dummy == NULL) {
		io_logging_memfatal(log, "Gadget header structure");
		return NULL;
	}

	/* Go to the header */
	skipsize = sizeof(int);
	if (f->ver == 2)
		skipsize += (sizeof(int) * 3 + 4*sizeof(char));
	fseek(f->file, skipsize, SEEK_SET);
	
	/* Now start reading the header */
	io_util_readint32(f->file, &(dummy->np[0]), f->swapped);
	io_util_readint32(f->file, &(dummy->np[1]), f->swapped);
	io_util_readint32(f->file, &(dummy->np[2]), f->swapped);
	io_util_readint32(f->file, &(dummy->np[3]), f->swapped);
	io_util_readint32(f->file, &(dummy->np[4]), f->swapped);
	io_util_readint32(f->file, &(dummy->np[5]), f->swapped);
	io_util_readdouble(f->file, &(dummy->massarr[0]), f->swapped);
	io_util_readdouble(f->file, &(dummy->massarr[1]), f->swapped);
	io_util_readdouble(f->file, &(dummy->massarr[2]), f->swapped);
	io_util_readdouble(f->file, &(dummy->massarr[3]), f->swapped);
	io_util_readdouble(f->file, &(dummy->massarr[4]), f->swapped);
	io_util_readdouble(f->file, &(dummy->massarr[5]), f->swapped);
	io_util_readdouble(f->file, &(dummy->expansion), f->swapped);
	io_util_readdouble(f->file, &(dummy->redshift), f->swapped);
	io_util_readint32(f->file, &(dummy->flagsfr), f->swapped);
	io_util_readint32(f->file, &(dummy->flagfeedback), f->swapped);
	io_util_readuint32(f->file, &(dummy->nall[0]), f->swapped);
	io_util_readuint32(f->file, &(dummy->nall[1]), f->swapped);
	io_util_readuint32(f->file, &(dummy->nall[2]), f->swapped);
	io_util_readuint32(f->file, &(dummy->nall[3]), f->swapped);
	io_util_readuint32(f->file, &(dummy->nall[4]), f->swapped);
	io_util_readuint32(f->file, &(dummy->nall[5]), f->swapped);
	io_util_readint32(f->file, &(dummy->flagcooling), f->swapped);
	io_util_readint32(f->file, &(dummy->numfiles), f->swapped);
	io_util_readdouble(f->file, &(dummy->boxsize), f->swapped);
	io_util_readdouble(f->file, &(dummy->omega0), f->swapped);
	io_util_readdouble(f->file, &(dummy->omegalambda), f->swapped);
	io_util_readdouble(f->file, &(dummy->hubbleparameter), f->swapped);
	io_util_readint32(f->file, &(dummy->flagstellarage), f->swapped);
	io_util_readint32(f->file, &(dummy->flagmetals), f->swapped);
	io_util_readuint32(f->file, &(dummy->nallhighw[0]), f->swapped);
	io_util_readuint32(f->file, &(dummy->nallhighw[1]), f->swapped);
	io_util_readuint32(f->file, &(dummy->nallhighw[2]), f->swapped);
	io_util_readuint32(f->file, &(dummy->nallhighw[3]), f->swapped);
	io_util_readuint32(f->file, &(dummy->nallhighw[4]), f->swapped);
	io_util_readuint32(f->file, &(dummy->nallhighw[5]), f->swapped);
	io_util_readint32(f->file, &(dummy->flagentropyu), f->swapped);

	/* Read the fillheader */
	io_util_readstring(f->file, dummy->unused, GADGET_HEADER_FILLHEADER);

	f->header = dummy;

	return dummy;
}

extern void
io_gadget_header_del(io_logging_t log, io_gadget_header_t *header)
{
	if ( (header == NULL) || (*header == NULL) )
		return;

	free(*header);

	*header = NULL;

	return;
}

extern void
io_gadget_header_write(io_logging_t log,
                       io_gadget_header_t header,
                       io_gadget_t f)
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
io_gadget_header_log(io_logging_t log, io_gadget_header_t header)
{
	io_logging_msg(log, INT32_C(5),
	               "Headerobject information:");
	io_logging_msg(log, INT32_C(5),
	               "  No of particles in file:       %" PRIi32,
	               header->np[0]);
	io_logging_msg(log, INT32_C(5),
	               "                                 %" PRIi32,
	               header->np[1]);
	io_logging_msg(log, INT32_C(5),
	               "                                 %" PRIi32,
	               header->np[2]);
	io_logging_msg(log, INT32_C(5),
	               "                                 %" PRIi32,
	               header->np[3]);
	io_logging_msg(log, INT32_C(5),
	               "                                 %" PRIi32,
	               header->np[4]);
	io_logging_msg(log, INT32_C(5),
	               "                                 %" PRIi32,
	               header->np[5]);
	io_logging_msg(log, INT32_C(5),
	               "  Mass of particle species:      %e",
	               header->massarr[0]);
	io_logging_msg(log, INT32_C(5),
	               "                                 %e",
	               header->massarr[1]);
	io_logging_msg(log, INT32_C(5),
	               "                                 %e",
	               header->massarr[2]);
	io_logging_msg(log, INT32_C(5),
	               "                                 %e",
	               header->massarr[3]);
	io_logging_msg(log, INT32_C(5),
	               "                                 %e",
	               header->massarr[4]);
	io_logging_msg(log, INT32_C(5),
	               "                                 %e",
	               header->massarr[5]);
	io_logging_msg(log, INT32_C(5),
	               "  Expansion:                     %e",
	               header->expansion);
	io_logging_msg(log, INT32_C(5),
	               "  Redshift:                      %e",
	               header->redshift);
	io_logging_msg(log, INT32_C(5),
	               "  Flagsfr:                       %" PRIi32,
	               header->flagsfr);
	io_logging_msg(log, INT32_C(5),
	               "  Flag Feedback:                 %" PRIi32,
	               header->flagfeedback);
	io_logging_msg(log, INT32_C(5),
	               "  No of particles in total:      %" PRIu32,
	               header->nall[0]);
	io_logging_msg(log, INT32_C(5),
	               "                                 %" PRIu32,
	               header->nall[1]);
	io_logging_msg(log, INT32_C(5),
	               "                                 %" PRIu32,
	               header->nall[2]);
	io_logging_msg(log, INT32_C(5),
	               "                                 %" PRIu32,
	               header->nall[3]);
	io_logging_msg(log, INT32_C(5),
	               "                                 %" PRIu32,
	               header->nall[4]);
	io_logging_msg(log, INT32_C(5),
	               "                                 %" PRIu32,
	               header->nall[5]);
	io_logging_msg(log, INT32_C(5),
	               "  Flag Cooling:                  %" PRIi32,
	               header->flagcooling);
	io_logging_msg(log, INT32_C(5),
	               "  Number of files:               %" PRIi32,
	               header->numfiles);
	io_logging_msg(log, INT32_C(5),
	               "  Boxsize:                       %e",
	               header->boxsize);
	io_logging_msg(log, INT32_C(5),
	               "  Omega0:                        %e",
	               header->omega0);
	io_logging_msg(log, INT32_C(5),
	               "  OmegaLambda:                   %e",
	               header->omegalambda);
	io_logging_msg(log, INT32_C(5),
	               "  Hubble parameter:              %e",
	               header->hubbleparameter);
	io_logging_msg(log, INT32_C(5),
	               "  Flag Stellar Age:              %" PRIi32,
	               header->flagstellarage);
	io_logging_msg(log, INT32_C(5),
	               "  Flag Metals:                   %" PRIi32,
	               header->flagmetals);
	io_logging_msg(log, INT32_C(5),
	               "  No of particles in total (HW): %" PRIu32,
	               header->nallhighw[0]);
	io_logging_msg(log, INT32_C(5),
	               "                                 %" PRIu32,
	               header->nallhighw[1]);
	io_logging_msg(log, INT32_C(5),
	               "                                 %" PRIu32,
	               header->nallhighw[2]);
	io_logging_msg(log, INT32_C(5),
	               "                                 %" PRIu32,
	               header->nallhighw[3]);
	io_logging_msg(log, INT32_C(5),
	               "                                 %" PRIu32,
	               header->nallhighw[4]);
	io_logging_msg(log, INT32_C(5),
	               "                                 %" PRIu32,
	               header->nallhighw[5]);
	io_logging_msg(log, INT32_C(5),
	               "  Flag Entropy insted U:         %" PRIi32,
	               header->flagentropyu);

	return;
}


/**********************************************************************\
 *    Implementation of local functions                               * 
\**********************************************************************/

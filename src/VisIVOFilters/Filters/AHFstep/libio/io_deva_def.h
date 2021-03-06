#ifndef IO_DEVA_DEF_H
#define IO_DEVA_DEF_H

/**
 * \file io_deva_def.h
 *
 * Provides the structure definition for the DEVA file structure.
 */


/**********************************************************************\
 *    Includes                                                        * 
\**********************************************************************/
#include <stdio.h>
#ifdef HAVE_STDINT_H
#	include <stdint.h>
#else
#	include <replace_stdint.h>
#endif

#include "io_deva_header_def.h"
#include "io_file.h"


/**********************************************************************\
 *    Global defines, structure definitions and typedefs              * 
\**********************************************************************/

/**
 * The file structure itself
 */
struct io_deva_struct {
	/** Holds the filetype ID, *must* be first in the structure */
	io_file_type_t ftype;
#ifdef WITH_MPI
	/** The global rank of the process */
	int rank;
	/** The size of the global communicator */
	int size;
	/** Stores the communicator used for intra libio communication */
	MPI_Comm mycomm;
	/** The size of the intra-library communicator */
	int size_mycomm;
	/** The rank of the local process */
	int rank_mycomm;
#endif
	/** Holds the file stream of the currently opened file part */
	FILE *file;
	/** Holds the filenames of the parts, array of length parts */
	char *fname;
	/** In which mode should this be opened */
	io_file_mode_t mode;
	/** Is this file byteswapped? */
	io_file_swap_t swapped;
	/** DEVA version */
	int8_t ver;
	/** DEVA header size */
	int32_t headsize;
	/** Stores the total number of particles in the file */
	uint64_t no_part;
	/** Stores the total number of particles in the file which have
	 *   their mass specified in a mass block */
	uint64_t no_part_with_mass;
	/** Flag if this is a multimass file */
	int8_t multimass;
	/** Stores the minimal weight of halo particles in file units */
	double mmass;
	/** Stores the minimal particle weight in file units */
	double minweight;
	/** Stores the maximal particle weight in file units */
	double maxweight;
	/** Stores the sum over all particle weights in file units */
	double sumweight;
	/** Stores the number of weight types */
	int32_t no_species;
	/** Stores the minimal position values in file units */
	double minpos[3];
	/** Stores the maximal position values in file units */
	double maxpos[3];
	/** Scaling value to convert from internal units to Mpc */
	double posscale;
	/** Scaling value to convert from internal mass to Msun */
	double weightscale;
	/** Holds the header information for each file part */
	io_deva_header_t header;
};

/** Convenient typedef */
typedef struct io_deva_struct io_deva_struct_t;

/** Convenient typedef */
typedef io_deva_struct_t *io_deva_t;


#endif /* IO_DEVA_DEF_H */


/**
 * \file io_tipsy.c
 *
 * Provides functions for reading and writing TIPSY files.
 */

/*
 *
 *  THIS ONLY CAME INTO EXISTENCE THANKS TO THE GREAT HELP OF JUSTIN READ!
 *
 *
 */

/**********************************************************************\
 *    Includes                                                        * 
\**********************************************************************/
#include <inttypes.h>
#include <stdlib.h>
#include <string.h>
#if (NONC99 == 0)
#	include <math.h>
#else
#	include <replace_math.h>
#endif
#include <stddef.h>
#ifdef WITH_MPI
#	include <mpi.h>
#endif

#include "io_tipsy.h"
#include "io_tipsy_header.h"
#include "io_util.h"


/**********************************************************************\
 *    Local defines, structure definitions and typedefs               * 
\**********************************************************************/
#define H0              100.  // not known as param.h is not included

// internally used bytes for various particle types
#define TIPSY_GASVAR    12
#define TIPSY_DARKVAR    9
#define TIPSY_STARVAR   11

// ptype values as used by AHF
// NOTE: this has to be identical to the definitions found in param.h
#define PGAS               ( 0.0)   /* identifier for gas particles; has to be exactly 0.0!!!! */
#define PDM                (-1.0)   /* identifier for dm particles; whatever negative value */
#define PSTAR              (-4.0)   /* identifier for star particles; whatever negative value */
#define PDMbndry           (-5.0)   /* identifier for star particles; whatever negative value */

/**********************************************************************\
 *    Prototypes of local functions                                   * 
\**********************************************************************/

/**
 * \brief Helper function to open the file
 *
 * This function is supposed to get inlined anyway. Makes
 * io_tipsy_open more readable.
 *
 * \param log   The logging object.
 * \param f     The TIPSY file object sofar
 * \param mode  The mode in which to open the file.
 *
 * \returns In case of an error NULL is returned, otherwise f.
 */
inline static io_tipsy_t
local_openopen(io_logging_t log, io_tipsy_t f, io_file_mode_t mode);

 /**
  * \brief Helper funtion to set the swap-state and write some
  *        messages
  *
  * \param log      The logging object.
  * \param f        The TIPSY file object sofar.
  * \param swapped  The swap state.
  *
  * \return Nothing.
  */
inline static void
local_openswapped(io_logging_t log,
                  io_tipsy_t f,
                  io_file_swap_t swapped);

/**
 * \brief Try to find out swapping status.
 *
 * \param log  The logging object.
 * \param f    The file object.
  *
  * \return Returns the file object or NULL in case of an error.
 */
inline static io_tipsy_t
local_opengetswap(io_logging_t log, io_tipsy_t f);

/**
 * \brief Check before writing particles and place file pointer at the
 *        right spot.
 *
 * \param log
 * \param f
 * \param pskip
 * \param bytes
 *
 * \return 
 */
inline static int32_t
local_write_common(io_logging_t log,
                   io_tipsy_t f,
                   uint64_t pskip,
                   int32_t bytes);


/**********************************************************************\
 *    Implementation of global functions                              * 
\**********************************************************************/
extern io_tipsy_t
io_tipsy_open(io_logging_t log,
               char *fname,
               io_file_swap_t swapped,
               io_file_mode_t mode,
               uint32_t reader)
{
	io_tipsy_t f;

	/* Get memory for the structure */
	f = (io_tipsy_t)malloc(sizeof(io_tipsy_struct_t));
	if (f == NULL) {
		io_logging_memfatal(log,  "io_tipsy structure");
		return NULL;
	}

	/* Start filling the structure */

	/* Store the filename */
	f->fname = (char *)malloc(sizeof(char) * (strlen(fname) + 1));
	if (f->fname == NULL) {
		io_logging_memfatal(log, "filename of TIPSYfile");
		free(f);
		return NULL;
	}
	strncpy(f->fname, fname, strlen(fname)+1);

	/* Okay, we are a TIPSY file */
	f->ftype = IO_FILE_TIPSY;

	/* And we can just copy in the parallel information */
#	ifdef WITH_MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &(f->rank));
	MPI_Comm_size(MPI_COMM_WORLD, &(f->size));
	if (f->size >= reader) {
		/* TODO 
		 * THIS IS JUST A QUICK HACK TO PREVENT TIPSY FILES TO
		 * AGAIN TRY TO SPLIT THE COMMUNICATOR, THAT IS ALREADY DONE
		 * IN mtipsy.c 
		 * TODO 
		 */
		MPI_Comm_split(MPI_COMM_WORLD, 1,
		               f->rank, &(f->mycomm));
		MPI_Comm_size(f->mycomm, &(f->size_mycomm));
		MPI_Comm_rank(f->mycomm, &(f->rank_mycomm));
	} else {
		f->mycomm = MPI_COMM_NULL;
		f->size_mycomm = -1;
		f->rank_mycomm = -1;
	}
#	endif

	/* Try to open the file and set mode */
	if (local_openopen(log, f, mode) == NULL) {
		free(f->fname);
		free(f);
		return NULL;
	}

	/* Set swapping */
	local_openswapped(log, f, swapped);
	if (    (f->mode == IO_FILE_READ)
	     && (f->swapped == IO_FILE_UNKOWN_SWAPPING)) {
		if (local_opengetswap(log, f) != f) {
			io_logging_fatal(log, "Cannot open this file.");
			free(f->fname);
			free(f);
			return NULL;
		}
	}

	/* Identify TIPSY format */
   f->ver = 1; // 1 = refers to raw simulation converted to single precision though...

	/* Nothing for the header for now */
	f->header = NULL;

	/* Set some dummy values */
	f->no_part = UINT64_C(0);
	f->no_part_with_mass = UINT64_C(0);
	f->multimass = INT8_C(0);
	f->mmass = 1e40;
	f->minweight = 1e40;
	f->maxweight = 0.0;
	f->sumweight = 0.0;
	f->no_species = INT32_C(0);
	f->posscale = 1.0;
	f->weightscale = 1.0;
  
	return f;
}

extern void
io_tipsy_close(io_logging_t log,
                io_tipsy_t *f)
{
	/* Catch NULLs */
	if (f == NULL || *f == NULL)
		return;

	/* Put header to the file if necessary */
	if (    ((*f)->mode == IO_FILE_WRITE)
	     && ((*f)->header != NULL)) {
		io_tipsy_header_write(log, (*f)->header, *f);
	}

	/* Close */
	if ((*f)->header != NULL)
		io_tipsy_header_del(log, &((*f)->header));
	if ((*f)->fname != NULL)
		free((*f)->fname);
#	ifdef WITH_MPI
	if ((*f)->mycomm != MPI_COMM_NULL)
		MPI_Comm_free(&((*f)->mycomm));
#	endif

	/* Actually close the file */
	if ((*f)->file != NULL)
		fclose((*f)->file);

	/* Cleaning */
	free(*f);
	*f = NULL;

	return;
}

extern void
io_tipsy_init(io_logging_t log,
               io_tipsy_t f)
{
  if (f == NULL)
    return;
  
  if (f->header != NULL) {
    io_logging_warn(log, INT32_C(1),
                    "Already have the header information! Rereading.");
    io_tipsy_header_del(log, &(f->header));
  }
  
  if (f->mode != IO_FILE_READ) {
    io_logging_warn(log, INT32_C(1),
                    "%s is not opened for reading. "
                    "Will do nothing.",
                    f->fname);
    return;
  }
  
  io_logging_msg(log, INT32_C(5),
                 "Starting to initialize file object from %s",
                 f->fname);
  f->header = io_tipsy_header_get(log, f);
  io_logging_msg(log, INT32_C(5),
                 "Done with initializing file object from %s",
                 f->fname);
  
  f->multimass = 1;
  f->no_part = f->header->nbodies;
  f->no_part_with_mass = f->no_part;

  return;
}

extern uint64_t
io_tipsy_readpart(io_logging_t log,
                   io_tipsy_t f,
                   uint64_t pskip,
                   uint64_t pread,
                   io_file_strg_struct_t strg)
{
	uint64_t particles_read, tmp;
	double box[3], shift[3];
	double scale_pos, scale_mom, scale_weight;

  /* 
	 * First read the particles unscaled. This will set important
	 * scaling information
	 */
	particles_read = io_tipsy_readpart_raw(log, f, pskip, pread, strg);
  
  
  if (particles_read != pread) {
		return UINT64_C(0);
	}

	/* And do the scaling */
#ifdef WITH_MPI
	io_tipsy_scale_global(log, f->mycomm,  f->maxpos,
	                       f->minpos, &(f->mmass), &(f->minweight), &(f->maxweight), &(f->sumweight));
#endif
  
  tmp = io_tipsy_scale_particles(log, f->maxpos, f->minpos,
	                                f->header->TIPSY_BOXSIZE,
	                                f->header->time,
	                                f->posscale, f->header->TIPSY_VUNIT, f->mmass,
	                                particles_read, strg);
	if (tmp != particles_read) {
		return tmp;
	}

  /* Wow, we are done! */
	return particles_read;
}

/* No we define a bunch of macros to make life easier */
#define SKIP {io_util_readint32(f->file, &blocksize, f->swapped);}
#define SKIP2 {io_util_readint32(f->file, &blocksize2, f->swapped);}
#define CHECK_BLOCK {\
	if (blocksize != blocksize2) {\
		io_logging_fatal(log,\
		                 "The block boundaries (beginning: %" PRIi32\
		                 " end: %" PRIi32 ") are not identical. "\
		                 "Corrupt file?", blocksize, blocksize2);\
		return UINT64_C(0);\
	} else {\
		io_logging_msg(log, INT32_C(5),\
		               "Block claimed correctly to be %f MB long.", \
		               (float)(blocksize/1024./1024.));\
	}\
}
#define DESCRIBE_BLOCK {\
	if (f->ver == 2) {\
		SKIP;\
		io_util_readstring(f->file, str, (size_t)4);\
		io_util_readint32(f->file, &nextblocksize, f->swapped);\
		io_logging_msg(log, INT32_C(1),\
		               "Arrived at block %s, size of it will be %" \
		               PRIi32, str, nextblocksize);\
		SKIP2;\
		CHECK_BLOCK;\
	}\
}
#define CHECK_FLOATBYTES(bfile, bstore) {\
	if (bfile == sizeof(float)) { \
		io_logging_msg(log, INT32_C(1), \
		               "Obviously the file uses float for " \
		               "floating point values (%" PRIi32 " bytes).", \
		               bfile); \
	} else if (bfile == sizeof(double)) { \
		io_logging_msg(log, INT32_C(1), \
		               "Obviously the file uses double for " \
		               "floating point values (%" PRIi32 " bytes).", \
		               bfile); \
	} else { \
		io_logging_fatal(log, \
		                 "No clue what kind of floating point uses " \
		                 "%" PRIi32 " bytes. Aborting reading.", \
		                 bfile); \
		return UINT64_C(0); \
	}\
	if (bfile < bstore) { \
		io_logging_msg(log, INT32_C(1), \
		               "The floating point values in the file have " \
		               "less precision than the particle storage " \
		               "(%" PRIi32 " bytes vs. %" PRIi32 " bytes). " \
		               "No problem, will upcast.", \
		                bfile, bstore); \
	} else if (bfile > bstore) { \
		io_logging_warn(log, INT32_C(1), \
		                "The floating point values in the file have " \
		                "a higher precision than the particle storage " \
		                "(%" PRIi32 " bytes vs. %" PRIi32 " bytes). " \
		                "Will downcast, but precision might be lost.", \
		                bfile, bstore); \
	} \
}

/* We are going to need that a few times for book-keeping */
#	define BOOK_KEEPING {\
        if (   isgreater(fweight, oldfweight) || isless(fweight, oldfweight)) { \
            f->no_species++; \
            oldfweight = fweight; \
            if (isless(fweight,f->minweight)) f->minweight = fweight; \
            if (isgreater(fweight, f->maxweight)) f->maxweight = fweight; \
            if (isless(fweight, f->mmass) ) f->mmass = fweight; } \
}

extern uint64_t
io_tipsy_readpart_raw(io_logging_t log,
                       io_tipsy_t f,
                       uint64_t pskip,
                       uint64_t pread,
                       io_file_strg_struct_t strg)
{
	uint64_t i, j, psum;
	int ptype;
	uint32_t bytes_file, bytes_int_file;
	int32_t blocksize, blocksize2, partsize;
	int32_t nextblocksize;
	long skipsize;
	double fposx, fposy, fposz;
	double fmomx, fmomy, fmomz;
	double fweight=0., oldfweight=0.;
	double fu;
   double ftemp;
	uint64_t fid;
	float dummy;
   double e_fac;
	uint32_t dummy_int;
	char str[5];
	/** Used to figure out at which particle type we are */
	uint32_t curprtt;


	/* Check if we actually have to do something */
	if ( (f == NULL) || (f->header == NULL) )
		return UINT64_C(0);

   //============================================================================
   // NOTE: we are scaling the thermal energies to km/sec right here, right now!
   //============================================================================
   //e_fac = f->header->time*f->header->time/2.89440501823307;
   e_fac = 1.0/2.89440501823307;

  
  /* Make sure that we are at the beginning of the file */
	rewind(f->file);

  /*******************************************************************\
   *  Skip HEADER                                                    *
  \*******************************************************************/
  skipsize = TIPSY_HEADER_SIZE;
  fseek(f->file, skipsize, SEEK_SET);

  
  /*******************************************************************\
   *  Parallel I/O: what fraction should the current CPU read?       *
  \*******************************************************************/
  bytes_file = 4; // TIPSY file stores every variable (int, float) as a 4 byte number
  blocksize  = bytes_file * 1*f->no_part;
  CHECK_FLOATBYTES(bytes_file, strg.bytes_float);
  io_logging_msg(log, INT32_C(1),
                 "A total of %" PRIu64 " particle masses "
                 "with %" PRIi32 " bytes per float (%f MB total) "
                 "are stored.",
                 f->no_part, bytes_file,
                 (float)(blocksize/1024./1024.));
  
  /* Set the number of particles to loop over correctly */
  io_logging_msg(log, INT32_C(3),
                 "Asked to read %" PRIu64 " and to skip %" PRIu64
                 " particles. Checking those numbers.",
                 pread, pskip);
  if (pskip > f->no_part) {
    pskip = f->no_part;
    io_logging_msg(log, INT32_C(3),
                   "Cannot skip more than there is, will now "
                   "only skip %" PRIu64 " particles.",
                   pskip);
  }
  if ( f->no_part - pskip < pread ) {
    pread = f->no_part - pskip;
    io_logging_msg(log, INT32_C(3),
                   "Cannot read more than there is left after "
                   "skipping. Will only read %" PRIu64
                   "particles.", pread);
  }
  
  /*******************************************************************\
   *  SKIP to the first particle to be read                          *
  \*******************************************************************/
  // first particle lies in GAS block
  //==================================
  if(pskip <= f->header->nsph)
    {
      // simply skip over gas particles consisting of TIPSY_GASVAR variables
      partsize = TIPSY_GASVAR * bytes_file;
      fseek(f->file, partsize*pskip, SEEK_CUR);
    }
  // first particle lies in DARK block
  //===================================
  else if(pskip <= f->header->nsph+f->header->ndark)
    {
      // skip over all gas particles consisting of TIPSY_GASVAR variables
      partsize = TIPSY_GASVAR * bytes_file;
      fseek(f->file, partsize*f->header->nsph, SEEK_CUR);

      // skip over remaining dark particles consisting of TIPSY_DARKVAR variables
      partsize = TIPSY_DARKVAR * bytes_file;
      fseek(f->file, partsize*(pskip-f->header->nsph), SEEK_CUR);
    }
  // first particle lies in STAR block
  //===================================
  else if(pskip <= f->header->nsph+f->header->ndark+f->header->nstar)
    {
      // skip over all gas particles consisting of TIPSY_GASVAR variables
      partsize = TIPSY_GASVAR * bytes_file;
      fseek(f->file, partsize*f->header->nsph, SEEK_CUR);
      
      // skip over all dark particles consisting of TIPSY_DARKVAR variables
      partsize = TIPSY_DARKVAR * bytes_file;
      fseek(f->file, partsize*f->header->ndark, SEEK_CUR);
      
      // skip over remaining star particles consisting of TIPSY_STARVAR variables
      partsize = TIPSY_STARVAR * bytes_file;
      fseek(f->file, partsize*(pskip-f->header->nsph+f->header->ndark), SEEK_CUR);
    }

  
  /*******************************************************************\
   *  READ particles                                                 *
  \*******************************************************************/
  
  /* keep track of total mass across all particle species */
  f->sumweight = 0.;
  
  
  for (i=pskip; i<pskip+pread; i++)
    {
      /* Read mass */
      io_util_readfloat(f->file, &dummy, f->swapped);
      fweight = (double)dummy;
      
      /* Read positions */
      io_util_readfloat(f->file, &dummy, f->swapped);
      fposx = (double)dummy;
      io_util_readfloat(f->file, &dummy, f->swapped);
      fposy = (double)dummy;
      io_util_readfloat(f->file, &dummy, f->swapped);
      fposz = (double)dummy;
      
      /* Read velocities */
      io_util_readfloat(f->file, &dummy, f->swapped);
      fmomx = (double)dummy;
      io_util_readfloat(f->file, &dummy, f->swapped);
      fmomy = (double)dummy;
      io_util_readfloat(f->file, &dummy, f->swapped);
      fmomz = (double)dummy;
      
      // We are in the GAS block
      if(i < f->header->nsph)
        {
          fid = i + f->header->ndark + f->header->nstar - 1;

          io_util_readfloat(f->file, &dummy, f->swapped);  // rho
          io_util_readfloat(f->file, &dummy, f->swapped);  // temp
          ftemp = (double)dummy*e_fac;
          io_util_readfloat(f->file, &dummy, f->swapped);  // hsmooth
          io_util_readfloat(f->file, &dummy, f->swapped);  // metals
          io_util_readfloat(f->file, &dummy, f->swapped);  // phi
        }
      // We are in the DARK block
      else if(i < f->header->nsph+f->header->ndark)
        {
          fid = i - f->header->nsph;
          ftemp = PDM;

          io_util_readfloat(f->file, &dummy, f->swapped);  // eps
          io_util_readfloat(f->file, &dummy, f->swapped);  // phi
        }
      // We are in the STAR block
      else if(i < f->header->nsph+f->header->ndark+f->header->nstar)
        {
          fid = i - f->header->nsph;
          ftemp = PSTAR;
          
          io_util_readfloat(f->file, &dummy, f->swapped);  // metals
          io_util_readfloat(f->file, &dummy, f->swapped);  // tform
          io_util_readfloat(f->file, &dummy, f->swapped);  // eps
          io_util_readfloat(f->file, &dummy, f->swapped);  // phi
        }
      
      
      /* Store valuable values */
      *((float *)strg.weight.val) = (float)fweight;
      *((float *)strg.posx.val) = (float)fposx;
      *((float *)strg.posy.val) = (float)fposy;
      *((float *)strg.posz.val) = (float)fposz;
      *((float *)strg.momx.val) = (float)fmomx;
      *((float *)strg.momy.val) = (float)fmomy;
      *((float *)strg.momz.val) = (float)fmomz;
      *((float *)strg.u.val) = (float)ftemp;
      *((uint32_t *)strg.id.val) = (uint32_t)fid;
      
      /* Increment the pointers to the next particle */
      strg.posx.val = (void *)(((char *)strg.posx.val) + strg.posx.stride);
      strg.posy.val = (void *)(((char *)strg.posy.val) + strg.posy.stride);
      strg.posz.val = (void *)(((char *)strg.posz.val) + strg.posz.stride);
      strg.momx.val = (void *)(((char *)strg.momx.val) + strg.momx.stride);
      strg.momy.val = (void *)(((char *)strg.momy.val) + strg.momy.stride);
      strg.momz.val = (void *)(((char *)strg.momz.val) + strg.momz.stride);
      strg.weight.val = (void *)(((char *)strg.weight.val) + strg.weight.stride);
      strg.u.val = (void *)(((char *)strg.u.val) + strg.u.stride);
      strg.id.val = (void *)(((char *)strg.id.val) + strg.id.stride);
      
      /* Detect extreme positions */
      if (isless(fposx, f->minpos[0]))      f->minpos[0] = fposx;
      if (isless(fposy, f->minpos[1]))      f->minpos[1] = fposy;
      if (isless(fposz, f->minpos[2]))      f->minpos[2] = fposz;
      if (isgreater(fposx, f->maxpos[0]))   f->maxpos[0] = fposx;
      if (isgreater(fposy, f->maxpos[1]))   f->maxpos[1] = fposy;
      if (isgreater(fposz, f->maxpos[2]))   f->maxpos[2] = fposz;
      
      /* Keep track of total mass (in internal units) and various other important quantities */
      f->sumweight += fweight;
      BOOK_KEEPING;
    }
  

  /*******************************************************************\
   *  Done with reading the TIPSY file, yeah!                       *
  \*******************************************************************/
  
	/* Return the number of particles read */
	return pread;
}

/* Getting rid of the macros */
#undef SKIP
#undef SKIP2
#undef CHECK_BLOCK
#undef DESCRIBE_BLOCK
#undef CHECK_FLOATBYTES


extern uint64_t
io_tipsy_writepart(io_logging_t log,
                    io_tipsy_t f,
                    uint64_t pskip,
                    uint64_t pwrite,
                    io_file_strg_struct_t strg)
{
	return 0;
}

extern uint64_t
io_tipsy_writepart_ord(io_logging_t log,
                        io_tipsy_t f,
                        uint64_t pskip,
                        uint64_t pwrite,
                        void *nxt_part,
                        io_file_strg_struct_t strg)
{
	return UINT64_C(0);
}

extern bool
io_tipsy_get(io_logging_t log,
              io_tipsy_t f,
              io_file_get_t what,
              void *res)
{
  if ( (f == NULL) || (f->header == NULL) )
    return false;
  
  switch (what) {
    case IO_FILE_GET_NOPART_IN_FILE:
    case IO_FILE_GET_NOPART:
      *((long *)res) = (long)f->header->nbodies;
      break;
    case IO_FILE_GET_NOVPART:
      *((double *)res) = (double)f->sumweight;
      break;
    case IO_FILE_GET_NOSPECIES:
      *((int *)res) = (int)f->no_species;
      break;
    case IO_FILE_GET_BOXSIZE:
      *((double *)res) = f->header->TIPSY_BOXSIZE;
      break;
    case IO_FILE_GET_PMASS:
      *((double *)res) = f->header->TIPSY_MUNIT;
      break;
    case IO_FILE_GET_ZINITIAL:
      io_logging_warn(log, INT32_C(1),
                      "zinitial is not set in a TIPSY file, "
                      "using current redshift");
    case IO_FILE_GET_Z:
      *((double *)res) = 1./f->header->time-1.;
      break;
    case IO_FILE_GET_AINITIAL:
      io_logging_warn(log, INT32_C(1),
                      "ainitial is not set in a TIPSY file, "
                      "using current expansion.");
    case IO_FILE_GET_A:
      *((double *)res) = f->header->time;
      break;
    case IO_FILE_GET_OMEGA0:
      *((double *)res) = f->header->TIPSY_OMEGA0;
      break;
    case IO_FILE_GET_OMEGAL:
      *((double *)res) = f->header->TIPSY_LAMBDA0;
      break;
    case IO_FILE_GET_H:
      io_logging_warn(log, INT32_C(1),
                      "TIPSY files don't store the timestep. "
                      "Setting to 100.0");
      *((double *)res) = H0;
      break;
    case IO_FILE_GET_DOUBLE:
      io_logging_warn(log, INT32_C(1),
                      "TIPSY files don't store the use of "
                      "double precision. Assuming it is not "
                      "double precision.");
      *((int *)res) = 0;
      break;
    case IO_FILE_GET_MMASS:
        *((int *)res) = 1;
      break;
    case IO_FILE_GET_NOTSTEP:
      io_logging_warn(log, INT32_C(1),
                      "TIPSY files don't store the timestep. "
                      "Setting to 500");
      *((int32_t *)res) = 500; //TODO: is nstep important for TIPSY files?
      break;
    case IO_FILE_GET_TSTEP:
      io_logging_warn(log, INT32_C(1),
                      "TIPSY files don't store the timestep. "
                      "Setting to 0.0");
      *((double *)res) = 0.0;
      break;
    case IO_FILE_GET_HEADERSTR:
      io_logging_warn(log, INT32_C(1),
                      "TIPSY files don't have a header string. "
                      "Using a dummy one.");
      *((char **)res) = "No header string.";
      break;
    case IO_FILE_GET_MINWEIGHT:
      *((double *)res) = f->minweight;
      break;
    case IO_FILE_GET_MAXWEIGHT:
      *((double *)res) = f->maxweight; 
      break;
    default:
      io_logging_fatal(log, "Requesting something unkown in %s.",
                       __func__);
      return false;
  }
  
  return true;
}


extern void
io_tipsy_log(io_logging_t log, io_tipsy_t f)
{
	io_logging_msg(log, INT32_C(5),
	               "Fileobject information:");
	io_logging_msg(log, INT32_C(5),
	               "  Filetype:             %s",
	               io_file_typestr(f->ftype));
	io_logging_msg(log, INT32_C(5),
	               "  Filename:             %s",
	               f->fname);
	io_logging_msg(log, INT32_C(5),
	               "  Mode:                 %" PRIi8,
	               f->mode);
	io_logging_msg(log, INT32_C(5),
	               "  Swapping:             %" PRIi8,
	               f->swapped);
	io_logging_msg(log, INT32_C(5),
	               "  File version:         %" PRIi8,
	               f->ver);
	io_logging_msg(log, INT32_C(5),
	               "  Header size:          %" PRIi8,
	               f->ver);//TODO: ask Steffen why f->ver here
	io_logging_msg(log, INT32_C(5),
	               "  No. particles:        %" PRIu64,
	               f->no_part);
	io_logging_msg(log, INT32_C(5),
	               "  No. particles w/mass: %" PRIu64,
	               f->no_part_with_mass);
	io_logging_msg(log, INT32_C(5),
	               "  Multimass:            %" PRIi8,
	               f->multimass);
	io_logging_msg(log, INT32_C(5),
	               "  MMass (Halo parts):   %g",
	               f->mmass);
	io_logging_msg(log, INT32_C(5),
	               "  Minimal Weight:       %g",
	               f->minweight);
	io_logging_msg(log, INT32_C(5),
	               "  Maximal Weight:       %g",
	               f->maxweight);
	io_logging_msg(log, INT32_C(5),
	               "  Sum of all weights:   %g",
	               f->sumweight);
	io_logging_msg(log, INT32_C(5),
	               "  No. of species:       %" PRIi32,
	               f->no_species);
	io_logging_msg(log, INT32_C(5),
	               "  Position scale:       %g",
	               f->posscale);
	io_logging_msg(log, INT32_C(5),
	               "  Weight scale:         %g",
	               f->weightscale);
	io_tipsy_header_log(log, f->header);

	return;
}

extern void
io_tipsy_resetscale(io_logging_t log,
                     io_tipsy_t f,
                     double posscale,
                     double weightscale) {
	if (f == NULL)
		return;

	io_logging_msg(log, INT32_C(8),
	               "Old posscale: %g   New posscale: %g",
	               f->posscale, posscale);
	io_logging_msg(log, INT32_C(8),
	               "Old weightscale: %g   New weightscale: %g",
	               f->weightscale, weightscale);
	f->posscale = posscale;
	f->weightscale = weightscale;

	return;
}

extern uint64_t
io_tipsy_scale_particles(io_logging_t log,
                          double maxpos[],
                          double minpos[],
                          double boxsize,
                          double expansion,
                          double posscale,
                          double vunit,
                          double mmass,
                          uint64_t particles_read,
                          io_file_strg_struct_t strg)
{
	double box[3], shift[3];
	double scale_pos, scale_mom, scale_weight;
	uint64_t i;
  
	/* Now we can do the scaling */
	box[0] = fabs(maxpos[0] - minpos[0]);
	box[1] = fabs(maxpos[1] - minpos[1]);
	box[2] = fabs(maxpos[2] - minpos[2]);
	if (isgreater(box[0], boxsize)) {
		io_logging_warn(log, INT32_C(1),
		                "x-Separation of particles exceeds boxsize "
		                "(%g > %g). Aborting!",
		                box[0], boxsize);
     exit(0);
	}
	if (isgreater(box[1], boxsize)) {
		io_logging_warn(log, INT32_C(1),
		                "y-Separation of particles exceeds boxsize "
		                "(%g > %g). Aborting!",
		                box[1], boxsize);
     exit(0);
	}
	if (isgreater(box[2], boxsize)) {
		io_logging_warn(log, INT32_C(1),
		                "z-Separation of particles exceeds boxsize "
		                "(%g > %g). Aborting!",
		                box[2], boxsize);
     exit(0);
	}
	io_logging_msg(log, INT32_C(4),
	               "Extreme positions: xmin = %g  xmax = %g",
	               minpos[0], maxpos[0]);
	io_logging_msg(log, INT32_C(4),
	               "                   ymin = %g  ymax = %g",
	               minpos[1], maxpos[1]);
	io_logging_msg(log, INT32_C(4),
	               "                   zmin = %g  zmax = %g",
	               minpos[2], maxpos[2]);
	shift[0] = (isless(minpos[0], 0.0) ? -minpos[0] : 0.0);
	shift[1] = (isless(minpos[1], 0.0) ? -minpos[1] : 0.0);
	shift[2] = (isless(minpos[2], 0.0) ? -minpos[2] : 0.0);
	io_logging_msg(log, INT32_C(4),
	               "Applying shift: (%g, %g, %g)",
	               shift[0], shift[1], shift[2]);

  
  
  
  
  
  
  
	/*=================================================================*\ 
    *                    TIPSY UNIT SCALING                           *
    *             we only need to scale VELOCTIES!                    *
   \*=================================================================*/
  scale_pos      = 1.0;
  scale_mom      = expansion*expansion * vunit/(H0*boxsize);
  scale_weight   = 1.0;
  
	io_logging_msg(log, INT32_C(3),
	               "Scaling by:  positions:  %g", scale_pos);
	io_logging_msg(log, INT32_C(3),
	               "             velocities: %g", scale_mom);
	io_logging_msg(log, INT32_C(3),
	               "             weights:    %g", scale_weight);
  
  
  
  
  
  
  
  
  

	/* Define the actual scaling calls type independent */
#	define SCALE_CALL(type) {\
		*((type *)strg.posx.val) += (type)(shift[0]); \
		*((type *)strg.posx.val) *= (type)(scale_pos); \
		*((type *)strg.posy.val) += (type)(shift[1]); \
		*((type *)strg.posy.val) *= (type)(scale_pos); \
		*((type *)strg.posz.val) += (type)(shift[2]); \
		*((type *)strg.posz.val) *= (type)(scale_pos); \
		*((type *)strg.momx.val) *= (type)(scale_mom); \
		*((type *)strg.momy.val) *= (type)(scale_mom); \
		*((type *)strg.momz.val) *= (type)(scale_mom); \
		if (strg.weight.val != NULL) *((type *)strg.weight.val) *= (type)(scale_weight); \
		strg.posx.val = (void *)(((char *)strg.posx.val) + strg.posx.stride); \
		strg.posy.val = (void *)(((char *)strg.posy.val) + strg.posy.stride); \
		strg.posz.val = (void *)(((char *)strg.posz.val) + strg.posz.stride); \
		strg.momx.val = (void *)(((char *)strg.momx.val) + strg.momx.stride); \
		strg.momy.val = (void *)(((char *)strg.momy.val) + strg.momy.stride); \
		strg.momz.val = (void *)(((char *)strg.momz.val) + strg.momz.stride); \
		if (strg.weight.val != NULL) strg.weight.val = (void *)(((char *)strg.weight.val) + strg.weight.stride); \
}

	/* Do the scaling depending on the storage type */
	if (strg.bytes_float == sizeof(float)) {
		for (i=0; i<particles_read; i++) {
			SCALE_CALL(float);
		}
	} else if (strg.bytes_float == sizeof(double)) {
		for (i=0; i<particles_read; i++) {
			SCALE_CALL(double);
		}
	} else if (strg.bytes_float == sizeof(long double)) {
		for (i=0; i<particles_read; i++) {
			SCALE_CALL(long double);
		}
	} else {
		io_logging_fatal(log,
		                 "Don't know which floating point types "
		                 "has %" PRIi32 " bytes. Aborting read.",
		                 strg.bytes_float);
		return UINT64_C(0);
	}

	/* Clean the macro away */
#	undef SCALE_CALL
  
	/* And we are done */
	return particles_read;
}

#ifdef WITH_MPI
extern void
io_tipsy_scale_global(io_logging_t log,
                       MPI_Comm comm,
                       double *maxpos,
                       double *minpos,
                       double *mmass,
                       double *minweight,
                       double *maxweight,
                       double *sumweight)
{
  int size, rank;
  double buffer[3];
  
  io_logging_msg(log, INT32_C(5),
                 "Updating local scale values to global values.");
  MPI_Allreduce((void *)maxpos, (void *)buffer, 3,
                MPI_DOUBLE, MPI_MAX, comm);
  io_logging_msg(log, INT32_C(5),
                 "local : maxpos[0] = %g \t"
                 "maxpos[1] = %g \t"
                 "maxpos[2] = %g",
                 maxpos[0], maxpos[1], maxpos[2]);
  maxpos[0] = buffer[0];
  maxpos[1] = buffer[1];
  maxpos[2] = buffer[2];
  io_logging_msg(log, INT32_C(5),
                 "global: maxpos[0] = %g \t"
                 "maxpos[1] = %g \t"
                 "maxpos[2] = %g",
                 maxpos[0], maxpos[1], maxpos[2]);
  
  MPI_Allreduce((void *)minpos, (void *)buffer, 3,
                MPI_DOUBLE, MPI_MIN, comm);
  io_logging_msg(log, INT32_C(5),
                 "local : minpos[0] = %g \t"
                 "minpos[1] = %g \t"
                 "minpos[2] = %g",
                 minpos[0], minpos[1], minpos[2]);
  minpos[0] = buffer[0];
  minpos[1] = buffer[1];
  minpos[2] = buffer[2];
  io_logging_msg(log, INT32_C(5),
                 "global: minpos[0] = %g \t"
                 "minpos[1] = %g \t"
                 "minpos[2] = %g",
                 minpos[0], minpos[1], minpos[2]);
  
  MPI_Allreduce((void *)sumweight, (void *)buffer, 1,
                MPI_DOUBLE, MPI_SUM, comm);
  io_logging_msg(log, INT32_C(5), "local : sumweight = %g", *sumweight);
  *sumweight = buffer[0];
  io_logging_msg(log, INT32_C(5), "global: sumweight = %g", *sumweight);
  
  MPI_Allreduce((void *)minweight, (void *)buffer, 1,
                MPI_DOUBLE, MPI_MIN, comm);
  io_logging_msg(log, INT32_C(5), "local : minweight = %g", *minweight);
  *minweight = buffer[0];
  io_logging_msg(log, INT32_C(5), "global: minweight = %g", *minweight);
  
  MPI_Allreduce((void *)maxweight, (void *)buffer, 1,
                MPI_DOUBLE, MPI_MAX, comm);
  io_logging_msg(log, INT32_C(5), "local : maxweight = %g", *maxweight);
  *maxweight = buffer[0];
  io_logging_msg(log, INT32_C(5), "global: maxweight = %g", *maxweight);
  
  MPI_Allreduce((void *)mmass, (void *)buffer, 1,
                MPI_DOUBLE, MPI_MIN, comm);
  io_logging_msg(log, INT32_C(5), "local : mmass = %g", *mmass);
  *mmass = buffer[0];
  io_logging_msg(log, INT32_C(5), "global: mmass = %g", *mmass);
  
  return;
}
#endif /* WITH_MPI */


/**********************************************************************\
 *    Implementation of local functions                               * 
\**********************************************************************/

inline static io_tipsy_t
local_openopen(io_logging_t log, io_tipsy_t f, io_file_mode_t mode)
{
	if (mode == IO_FILE_READ) {
		f->file = fopen(f->fname, IO_FILE_MODE_READ);
		if (f->file == NULL) {
			io_logging_fatal(log,
			                 "Could not open '%s' for reading.",
			                 f->fname);
			return NULL;
		}
	} else {
		f->file = fopen(f->fname, IO_FILE_MODE_WRITE);
		if (f->file == NULL) {
			io_logging_fatal(log,
			                 "Could not open '%s' for writing.",
			                 f->fname);
			return NULL;
		}
	}

	f->mode = mode;

	return f;
}

inline static void
local_openswapped(io_logging_t log,
                  io_tipsy_t f,
                  io_file_swap_t swapped)
{
	if (f->mode == IO_FILE_READ) {
		switch(swapped) {
			case IO_FILE_ISNOT_SWAPPED:
				io_logging_msg(log, INT32_C(3),
				               "Assuming unswapped file");
				break;
			case IO_FILE_IS_SWAPPED:
				io_logging_msg(log, INT32_C(3),
				               "Assuming swapped file");
				break;
			case IO_FILE_UNKOWN_SWAPPING:
			default:
				io_logging_msg(log, INT32_C(3),
				               "Will try to find out swap status");
		}
	}

	/* Now set the swapping */
	f->swapped = swapped;

	return;
}

inline static io_tipsy_t
local_opengetswap(io_logging_t log, io_tipsy_t f)
{

  /* We simply let the user decide! */
  
#ifdef BYTESWAP
  f->swapped = IO_FILE_IS_SWAPPED;
#else
  f->swapped = IO_FILE_ISNOT_SWAPPED;
#endif
  
  
	return f;
}


inline static int32_t
local_write_common(io_logging_t log,
                   io_tipsy_t f,
                   uint64_t pskip,
                   int32_t bytes)
{
	return 0;
}

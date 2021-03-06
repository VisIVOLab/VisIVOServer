#ifndef SPLOTCH_CUDA2_H
#define SPLOTCH_CUDA2_H

#include "cuda/splotch_cuda.h"
#include "splotch/splotch_host.h"
#include "cxxsupport/walltimer.h"

THREADFUNC host_thread_func(void *pinfo);
THREADFUNC cu_thread_func(void *pinfo);


// struct containing thread task info
struct thread_info
  {
  int devID;                  //index of the device selected
  int startP, endP;           //start and end particles to handle
  long npart_all;             //total number of particles
  arr2<COLOUR> *pPic;         //output image computed 
  wallTimerSet times;
  };

//some global info shared by all threads
extern paramfile       *g_params;
extern std::vector<particle_sim> particle_data; //raw data from file
extern vec3 campos, lookat, sky;
extern std::vector<COLOURMAP> amap;
extern int ptypes;
extern wallTimerSet cuWallTimers;

int check_device(int rank);
void print_device_info(int rank, int dev);

#ifdef SPLVISIVO
void cuda_rendering(int mydevID, int nDev, arr2<COLOUR> &pic,VisIVOServerOptions &opt);
#else
void cuda_rendering(int mydevID, int nDev, arr2<COLOUR> &pic);
#endif
void DevideThreadsTasks(thread_info *tInfo, int nThread, bool bHostThread);
void cu_draw_chunk(void *pinfo, cu_gpu_vars* gv);
int filter_chunk(int StartP, int chunk_dim, int nParticle, int maxRegion,
                 int nFBufInCell, cu_particle_splotch *cu_ps,
                 cu_particle_splotch *cu_ps_filtered, int *End_cu_ps, 
                 int *nFragments2Render);
void render_chunk(int EndP, int nFBufInCell, cu_particle_splotch *cu_ps_filtered,                           void *fragBuf, cu_gpu_vars *gv, bool a_eq_e, float64 grayabsorb,
                  arr2<COLOUR> &pPic, wallTimerSet &times);
void combine_chunk(int StartP, int EndP, cu_particle_splotch *cu_ps_filtered, 
                  void *fragBuf, bool a_eq_e, float64 grayabsorb, arr2<COLOUR> &pPic);
void setup_colormap(int ptypes, cu_gpu_vars* gv);

void GPUReport(wallTimerSet &cuTimers);
void cuda_timeReport(paramfile &params);

#endif

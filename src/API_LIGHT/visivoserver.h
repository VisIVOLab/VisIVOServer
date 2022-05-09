/***************************************************************************
 *   Copyright (C) 2011 by VisIVO team *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#ifndef VISIVOSERVER_H
#define VISIVOSERVER_H

#ifndef WIN32
	#include <pthread.h>
	#include <sys/types.h>
	#include <sys/wait.h>
	#include <errno.h>
	#include <unistd.h>
	#include <stdlib.h>
	#include <stdio.h>
#endif


struct VisIVOFilter
{
  int setatt[1000];  // MUST BE EQUAL TO NPAR in visivodef.h QUI MASSIMO NUMERO PAR
  char op[64];
  char vbt[256], field[256], outcol[256], start[256], out[256], newcolnames[256];
  char outtxt[256];
  char pointcolumns[256], volume[256], infiles1[256], infiles2[256];
  char intfrom[20], intto[256], center[256], radius[256], vo[512], lfnout[512],se[512];
  char newres[256];
  char perc[20], volperc[20], threshold[64], startingcell[64], resolution[64];
  char gridorigin[64], gridspacing[64], interval[64];
  char seed[20], skip[20], box[20], numbin[20], outvalue[20], invalue[20], pad[20];
  char numrows[20], rangerowsfrom[64],rangerowsto[64],dataformat[20],width[20],precision[20];
  char numcells[20], nsigma[20], volumesplit[20], numoftables[20], maxsizetable[20];
  char visualsize[20], histobin[20];
  char **list;
  int nList;
  char appendList[256];
  char limitsList[256];
  char geometryList[256], MRgeometry[256], MRpos[256];
  char mergeList[256];
  char visualList[256];
  char math[512];
  char MRbackground[64];
  char dimvox[256];
  char trackplanedist[64];
  char innerdist[64];
  char outpoints[256];
  char outvol[256];

};

struct VisIVOImporter
{
  int setatt[1000];  // MUST BE EQUAL TO NPAR in visivodef.h QUI MASSIMO NUMERO PAR
  char fformat[64];
  char infile[256], outfile[256], userpwd[256], binaryheader[256];
  char datasetList[512], hyperslab[512], VO[512], lfnout[512], se[512];
  int comp[3];
  float size[3], missing, text;
  unsigned long long int npoints;
};

struct VisIVOViewer
{
  float *internalPointers[20];
  char internalNames[20][256];
  int setatt[1000];  // MUST BE EQUAL TO NPAR in visivodef.h QUI MASSIMO NUMERO PAR
  char path[256];
  char xField[256],yField[256],zField[256];
  char xVectorField[256],yVectorField[256],zVectorField[256];
  char colorScalar[256];
  char vRenderingField[256];
  char isosurfaceField[256];
  char vo[512], lfnout[512],se[512];
  char slicefield[256];
  char radiusscalar[256];
  char heightscalar[256];
  char slicePlane[10];
  char imageName[256];
  char imageSize[20];
  char backColor[20];
  char oneColor[20];
  char colorTable[256];
  char stereoMode[20];
  char anaglyphmask[20];
  char isoSmooth[20];
  char slicePlanePoint[64];
  char slicePlaneNormal[64];
  float colorRangeTo,colorRangeFrom;
  float anaglyphsat;
  int x,y,z;//!colums number associated to each axes 
  int numImage;//! id number of image
  int  numImageToLoad; //!number of image to generete
  int vx,vy,vz;//!colums number associated to each vectore axes
  int hsml,spColor,spIntensity;//!colums number associated used in splotch
  int nRows, nCols;//!numeber of rows and colums 
  int nRadius, nHeight,nColorScalar,nVRenderingField ,nSlicePlane,
  nIsosurfaceField;//!colums number associated to radius scaling,hiegth scaling, lut, rendering, slice and isosurface 
  int nGlyphs,nColorTable;//!set the number of glyphs and luttable 
  int nSliceField;//!set the number of plane that the user want for slice
  
  double elevation,zoom,azimuth;//! value for camera position default il 0, 0,1
  double campos[3], camfp[3];
  char cycleFile[256]; //! set if the system if little or big endian
  int cycleOffset; //! set if the system if little or big endian
//  int cycleSkipFrom; //! skip the first # lines from the begin in cycle file
//  int cycleSkipTo; //! read up to # line number in cycle file
 char glyphs[20];
 double opacity;//!value of opacity dafault is 0.666
  float opacityTF[3];
  double radius,height ;//! value of radius and hieght for glyphs. the usere can use this with scaling or not. default is 1 for both 
  double size[3];//! cell size of volume
  double comp[3];//! resolution of volume 
  char slicePosition[10];//! value for the position of slice 
  double isosurfaceValue;//! value for isosurface
  double vectorScalingFactor;
  int vectorScale;
//Splotch section
  char splotchpar[256];
  char labelIntensity[256];
  char labelColor[256];
  char labelhsml[256];
  double cliprange[2];
};

struct VBT
{
  char datatype[20];
  int nOfFields;
  unsigned long long int nOfRows;
  int nCells[3];
  float cellSize[3];
  char endianity[20];
  char **field;
  float **statistic; //max, min, avg, sum
  char locator[256];
};

struct VisIVOAsynchId
{
	#ifndef WIN32 
		pthread_t threadId;
		pid_t pid;
	#endif
	int state;
	int errorCode;
	void *env;
	void *envSecond;	// Usato solo per directImage
	float random;	// Usato solo per directImage
	int withThread;	// Usata solo per VV_View
};

typedef struct VisIVOAsynchId VisIVOAsynchId;
typedef struct VBT VBT;
typedef struct VisIVOViewer VisIVOViewer;
typedef struct VisIVOImporter VisIVOImporter;
typedef struct VisIVOFilter VisIVOFilter;



#endif

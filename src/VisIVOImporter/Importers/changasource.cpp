/***************************************************************************
 *   Copyright (C) 2008 by Gabriella Caniglia,Roberto Munzone *
 *  gabriella.caniglia@oact.inaf.it *
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
//#include "VisIVOImporterConfigure.h"
#include "changasource.h"

#include "visivoutils.h"
#include "mpi.h"
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <unordered_map>
#include <filesystem>
#include <fcntl.h>
#include <omp.h>
#include <unistd.h>
#include <rpc/xdr.h>
#include <rpc/rpc.h>

//---------------------------------------------------------------------
int ChangaSource::readHeader()
//---------------------------------------------------------------------

{
  char dummy[4]; 
//   char checkType[4];
  //int provided;
  //MPI_Init_thread(0, 0, MPI_THREAD_MULTIPLE, &provided);
  int type=0; unsigned int Npart=0;
  std::string systemEndianism;
  //umBlock=0;
  bool needSwap=false;
#ifdef VSBIGENDIAN
  systemEndianism="big";
#else
  systemEndianism="little";
#endif
  if((m_endian=="b" || m_endian=="big") && systemEndianism=="little")
    needSwap=true;
  if((m_endian=="l" || m_endian=="little") && systemEndianism=="big")
    needSwap=true;
	
  std::string fileName = m_pointsFileName.c_str();
  std::ifstream inFile;

  inFile.open(fileName, std::ios::binary);
  if (!inFile)
  {
    std::cerr<<"Error while opening File"<<std::endl;
    return -1;
  }
  //inFile.read((char *)(dummy), 4*sizeof(char)); //!*** IMPORTANT NOT REMOVE ***//
  //inFile.read((char *)(tmpType), 4*sizeof(char));   //!*** IMPORTANT NOT REMOVE ***//

  //tagType=tmpType;
  //checkType.push_back(tagType);
  
  header curHead;
  inFile.seekg(0, std::ios::beg);
  int pad;



  fpread = fopen( m_pointsFileName.c_str(), "r" );
  fread(&curHead,28,1,fpread);
  rewind(fpread);

  xdrstdio_create(&xdrread, fpread, XDR_DECODE);
  if(xdr_double(&xdrread, &((struct header *)&curHead)->time) != TRUE){
    return 1;
  }
  
  if(xdr_int(&xdrread, &((struct header *)&curHead)->nbodies) != TRUE)
    return 1;
  if(xdr_int(&xdrread, &((struct header *)&curHead)->ndim) != TRUE)
    return 1;
  if(xdr_int(&xdrread, &((struct header *)&curHead)->nsph) != TRUE)
    return 1;
  if(xdr_int(&xdrread, &((struct header *)&curHead)->ndark) != TRUE)
    return 1;
  if(xdr_int(&xdrread, &((struct header *)&curHead)->nstar) != TRUE)
    return 1;
  if(xdr_int(&xdrread, &pad) != TRUE)
    return 1;

  nsph = curHead.nsph;
  ndark = curHead.ndark;
  nstar = curHead.nstar;
  
  return 0;
}

int ChangaSource::readData()
{
  int idx = m_pointsBinaryName.rfind('.');
	std::string pathFileIn = m_pointsBinaryName.erase(idx, idx+4);
	std::string pathFileOut = pathFileIn;
  //gas_particle* gasParticles = (gas_particle*) malloc(sizeof(gas_particle)*nsph);
  std::ofstream outfile((pathFileOut + "GAS" + ".bin").c_str(),std::ofstream::binary );
  float* gasParticles = (float*) malloc(sizeof(float)*12*nsph);
  std::clog << sizeof(gas_particle)*nsph << std::endl;
  if (gasParticles == NULL) {
    std::clog << "Malloc Error" << std::endl;
    return 1;
  }
  float* tmp = gasParticles;
  std::clog << sizeof(gas_particle) << std::endl;
  for (int i=0; i<nsph; ++i) {
 // std::clog << i << std::endl;
    xdr_vector(&xdrread,(char *) tmp, (12),
                sizeof(float),(xdrproc_t) xdr_float);
    tmp+=12;
  }

  std::vector<std::string> gasBlocks; //!block fields names
  gasBlocks.push_back("MASS"); //0
  gasBlocks.push_back("POS_X"); //1
  gasBlocks.push_back("POS_Y");  //2
  gasBlocks.push_back("POS_Z");//3
  gasBlocks.push_back("VEL_X");   //4
  gasBlocks.push_back("VEL_y");//5
  gasBlocks.push_back("VEL_Z"); //6
  gasBlocks.push_back("RHO");  //7
  gasBlocks.push_back("TEMP");  //8
  gasBlocks.push_back("EPS");//9
  gasBlocks.push_back("METALS"); //10
  gasBlocks.push_back("PHI"); //11

  float *bufferBlock = (float*) malloc(sizeof(float)*nsph);
  for(int elem = 0; elem < 12; elem ++){
    for (int part=0; part<nsph; part++) {
      bufferBlock[part] = gasParticles[part*12 + elem];
    }
		outfile.write((char *)(bufferBlock), sizeof(float)*nsph);
  }

	outfile.close();

	std::string pathHeader = pathFileOut+"GAS"+ ".bin";
	makeHeader(nsph, pathHeader, gasBlocks, m_cellSize,m_cellComp,m_volumeOrTable);

  outfile.open((pathFileOut + "DARK" + ".bin").c_str(),std::ofstream::binary );
  
  float* darkParticles = (float*) malloc(sizeof(float)*9*ndark);
  if (darkParticles == NULL) {
    std::clog << "Malloc Error" << std::endl;
    return 1;
  }
  
  float *tmp2 = darkParticles;
  for (int i=0; i<ndark; ++i) {
 // std::clog << i << std::endl;
    xdr_vector(&xdrread,(char *) tmp2, (9),
                sizeof(float),(xdrproc_t) xdr_float);
    tmp2+=9;
  }


  std::vector<std::string> darkBlocks; //!block fields names
  darkBlocks.push_back("MASS"); //0
  darkBlocks.push_back("POS_X"); //1
  darkBlocks.push_back("POS_Y");  //2
  darkBlocks.push_back("POS_Z");//3
  darkBlocks.push_back("VEL_X");   //4
  darkBlocks.push_back("VEL_y");//5
  darkBlocks.push_back("VEL_Z"); //6
  darkBlocks.push_back("EPS");//7
  darkBlocks.push_back("PHI"); //8

  float* bufferBlock2 = (float*) malloc(sizeof(float)*ndark);
  for(int elem = 0; elem < 9; elem ++){
    for (int part=0; part<ndark; part++) {
      bufferBlock2[part] = darkParticles[part*9 + elem];
    }
		outfile.write((char *)(bufferBlock2), sizeof(float)*ndark);
  }
  pathHeader = pathFileOut+"DARK"+ ".bin";
	makeHeader(ndark, pathHeader, darkBlocks, m_cellSize,m_cellComp,m_volumeOrTable);

	outfile.close();

  outfile.open((pathFileOut + "STAR" + ".bin").c_str(),std::ofstream::binary );
  
  float* starParticles = (float*) malloc(sizeof(float)*11*nstar);
  if (starParticles == NULL) {
    std::clog << "Malloc Error" << std::endl;
    return 1;
  }
  
  float *tmp3 = starParticles;
  for (int i=0; i<nstar; ++i) {
 // std::clog << i << std::endl;
    xdr_vector(&xdrread,(char *) tmp3, (11),
                sizeof(float),(xdrproc_t) xdr_float);
    tmp3+=11;
  }


  std::vector<std::string> starBlocks; //!block fields names
  starBlocks.push_back("MASS"); //0
  starBlocks.push_back("POS_X"); //1
  starBlocks.push_back("POS_Y");  //2
  starBlocks.push_back("POS_Z");//3
  starBlocks.push_back("VEL_X");   //4
  starBlocks.push_back("VEL_y");//5
  starBlocks.push_back("VEL_Z"); //6
  starBlocks.push_back("METALS");//7
  starBlocks.push_back("TFORM"); //8
  starBlocks.push_back("EPS");//9
  starBlocks.push_back("PHI"); //10

  float* bufferBlock3 = (float*) malloc(sizeof(float)*nstar);
  for(int elem = 0; elem < 11; elem ++){
    for (int part=0; part<nstar; part++) {
      bufferBlock3[part] = starParticles[part*11 + elem];
    }
		outfile.write((char *)(bufferBlock3), sizeof(float)*nstar);
  }
  pathHeader = pathFileOut+"STAR"+ ".bin";
	makeHeader(nstar, pathHeader, starBlocks, m_cellSize,m_cellComp,m_volumeOrTable);

	outfile.close();
  
  return 1;
}

ChangaSource::~ChangaSource() {
}

ChangaSource::ChangaSource() {
}
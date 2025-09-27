/***************************************************************************
 *   Copyright (C) 2008 by Ugo Becciani   *
 *   ugo.becciani@oact.inaf.it   *
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
#include "visivodef.h"
#include <cstdlib>
#include <cstring>

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include "parametersparser.h"
#include "startFilter.h"
#ifdef VSMPI
#pragma message "MPI-PARALLEL compilation"
#include "mpi.h"
#else
#pragma message "SERIAL compilation"
#endif

extern "C" {

// Dichiarazione delle funzioni principali.
int VF_Filter(VisIVOFilter *env);
void* VF_Filter_Thread(void *id);

/*	Questa funzione crea il thread è ritorna o un numero >= di 0 che rappresenta il thread creato,
 *  oppure un valore minore di zero che rappresenta l'errore che si è verificato.
 * 
 * -1	Massimo numero di thread creati raggiunto. E' necessario liberare qualche risorsa.
 * -2	Creazione del thread fallita.
 * 
 */
//----------------------------
int VA_Filter(VisIVOFilter *e, VisIVOAsynchId *id)
//----------------------------
{
	int iret;
	id->withThread=1;
	#ifdef WIN32
 		iret=VF_Filter(e);
 		return iret;
	#endif
	#ifndef WIN32
	id->env=e;
	id->state=runningThread;
	iret=pthread_create(&(id->threadId), NULL, &VF_Filter_Thread, (void *)id);
	if(iret!=0)	
	{
		id->state=undefined;
		iret=errorThread;
	}
	else
	{
		iret=noError;
	}
	return iret;
	#endif
}

//---------------------------
void *VF_Filter_Thread(void *id)
//---------------------------
{
	int iret;
	VisIVOAsynchId *idTh;
	idTh=(VisIVOAsynchId *)id;
	iret=VF_Filter((VisIVOFilter *)(idTh->env));
	idTh->errorCode=iret;
	idTh->state=successfulEndThread;
	return 0;	
}


//---------------------------
int VF_Filter(VisIVOFilter *env)
//---------------------------
{
std::string filename;
std::map<std::string, std::string> appParameters;
std::map<std::string, std::string>::iterator iter;
int rank=0, size=1;
// FILLING appParameters
std::string key;
std::string value;
#ifdef VSMPI
    MPI_Init (&argc, &argv);
    MPI_Comm_size (MPI_COMM_WORLD, &size);
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    // NO MPI ALLOWED at the moment
    if(rank>0) rank=-1;
    size=1;
	int *ranks;
	ranks= new int[size];
    for(int i=0;i<size;i++) ranks[i]=i;

    // create a new communicator
    MPI_Group origGroup, newGroup;
    MPI_Comm NEW_COMM;    
    MPI_Comm_group(MPI_COMM_WORLD, &origGroup);
    MPI_Group_incl(origGroup,size,ranks,&newGroup);
    MPI_Comm_create(MPI_COMM_WORLD, newGroup, &NEW_COMM); 
    // it is only the newGroup in NEW_COMM that will continue. 
#endif
    
for(int idPar=0; idPar<NPAR; idPar++)
{
  if(env->setatt[idPar]==1)
  {
//    std::clog<<" i="<<idPar<<" setatt "<<env->setatt[idPar]<<std::endl;
    key="",
    value="";
    switch(idPar+VF_PARAM)
    {
      case VF_SET_OPERATION:
      {
	key="op"; value=env->op;
	break;
      }
      case VF_SET_FILEVBT:
      {
	key="file"; value=env->vbt;
	break;
      }
      case VF_SET_RANDOMPERC:
      {
	key="perc"; value=env->perc;
	break;
      }
      case VF_SET_RANDOMSEED:
      {
	key="iseed"; value=env->seed;
	break;
      }
      case VF_SET_FIELD:
      {
	key="field"; value=env->field;
	break;
      }
      case VF_SET_OUTVBT:
      {
	key="out"; value=env->out;
	break;
      }
      case VF_SET_OUTCOL:
      {
	key="outcol"; value=env->outcol;
	break;
      }
      case VF_SET_ADDIDSTART:
      {
	key="start"; value=env->start;
	break;
      }
      case VF_SET_APPENDLIST:
      {
	key="filelist"; value=env->appendList;
	break;
      }
      case VF_SET_APPEND:
      {
	key="append"; value="";
	break;
      }
      case VF_SET_NEWCOLNAMES:
      {
	key="newnames"; value=env->newcolnames;
	break;
      }
      case VF_SET_VOLUMEPERC:
      {
	key="perc"; value=env->volperc;
	break;
      }
      case VF_SET_NEWRES:
      {
	key="newres"; value=env->newres;
	break;
      }
      case VF_SET_LIMITS:
      {
	key="limits"; value=env->limitsList;
	break;
      }
      case VF_SET_OPERATOROR:
      {
	key="operator"; value="OR";
	break;
      }
      case VF_SET_OPERATORAND:
      {
	key="operator"; value="AND";
	break;
      }
      case VF_SET_CUTTHRESHOLD:
      {
	key="threshold"; value=env->threshold;
	break;
      }
      case VF_SET_DECIMATORSKIP:
      {
	key="skip"; value=env->skip;
	break;
      }
      case VF_SET_EXTRACTIONGEOMETRY:
      {
	key="geometry"; value=env->geometryList;
	break;
      }
      case VF_SET_STARTINGCELL:
      {
	key="startingcell"; value=env->startingcell;
	break;
      }
      case VF_SET_RESOLUTION:
      {
	key="resolution"; value=env->resolution;
	break;
      }
      case VF_SET_POINTCOLUMNS:
      {
	key="points"; value=env->pointcolumns;
	break;
      }
      case VF_SET_DENSITY:
      {
	key="density"; value="";
	break;
      }
      case VF_SET_ALGORITHMTSC:
      {
	key="tsc"; value="";
	break;
      }
      case VF_SET_ALGORITHMNGP:
      {
	key="ngp"; value="";
	break;
      }
      case VF_SET_ALGORITHMCIC:
      {
	key="cic"; value="";
	break;
      }
      case VF_SET_VOLUME:
      {
	key="volume"; value=env->volume;
	break;
      }
      case VF_SET_GRIDORIGIN:
      {
	key="gridOrigin"; value=env->gridorigin;
	break;
      }
      case VF_SET_GRIDSPACING:
      {
	key="gridSpacing"; value=env->gridspacing;
	break;
      }
      case VF_SET_BOX:
      {
	key="box"; value=env->box;
	break;
      }
      case VF_SET_PERIODIC:
      {
	key="periodic"; value="";
	break;
      }
      case VF_SET_NUMBIN:
      {
	key="numbin"; value=env->numbin;
	break;
      }
       case VF_SET_INTERVAL:
      {
	key="numbin"; 
	value=env->intfrom;
	value+=" ";
	value+=env->intto;
	break;
      }
       case VF_SET_INFILES:
      {
	key="infiles"; 
	value=env->infiles1;
	value+=" ";
	value+=env->infiles2;
	break;
      }
       case VF_SET_CENTER:
      {
	key="center"; value=env->center;
	break;
      }
       case VF_SET_RADIUS:
      {
	key="radius"; value=env->radius;
	break;
      }
      case VF_SET_INVALUE:
      {
	key="invalue";  value=env->invalue;
	break;
      }
      case VF_SET_OUTVALUE:
      {
	key="outvalue";   value=env->outvalue;
	break;
      }
      case VF_SET_MATEXPRESSION:
      {
	key="compute";   value=env->math;
	break;
      }
      case VF_SET_MERGEHUGE:
      {
	key="size";   value="HUGE";
	break;
      }
      case VF_SET_MERGESMALL:
      {
	key="size";   value="SMALLEST";
	break;
      }
      case VF_SET_MERGEPAD:
      {
	key="size";   value=env->pad;
	break;
      }
      case VF_SET_MERGELIST:
      {
	key="filelist";   value=env->mergeList;
	break;
      }
      case VF_SET_NODENSITY:
      {
	key="nodensity";   value="";
	break;
      }
      case VF_SET_AVG:
      {
	key="avg";   value="";
	break;
      }
      case VF_SET_DELETECOLUMNS:
      {
	key="delete";   value="";
	break;
      }
       case VF_SET_NUMROWS:
      {
	key="numrows";   value=env->numrows;
	break;
      }
       case VF_SET_RANGEROWS:
      {
	key="rangerows";
	value=env->rangerowsfrom;
	value+=" ";
	value+=env->rangerowsto;
	break;
      }
       case VF_SET_WIDTH:
      {
	key="rangerows"; value=env->width;
	break;
      }
       case VF_SET_PRECISION:
      {
	key="rangerows"; value=env->precision;
	break;
      }
       case VF_SET_OUT:
      {
	key="out"; value=env->outtxt;
	break;
      }
       case VF_SET_NUMCELLS:
      {
	key="numcells"; value=env->numcells;
	break;
      }
       case VF_SET_NSIGMA:
      {
	key="nsigma"; value=env->nsigma;
	break;
      }
       case VF_SET_ALLCOLUMNS:
      {
	key="allcolumns"; value="";
	break;
      }
       case VF_SET_VOLUMESPLIT:
      {
	key="volumesplit"; value=env->volumesplit;
	break;
      }
       case VF_SET_NUMOFTABLES:
      {
	key="numoftables"; value=env->numoftables;
	break;
      }
       case VF_SET_MAXSIZETABLE:
      {
	key="numoftables"; value=env->maxsizetable;
	break;
      }
       case VF_SET_HUGESPLIT:
      {
	key="hugesplit"; value="";
	break;
      }    
       case VF_SET_SWAPOVERRIDE:
      {
	key="override"; value="";
	break;
      }    
       case VF_SET_VISUALSIZE:
      {
	key="size"; value=env->visualsize;
	break;
      }    
      case VF_SET_VISUALLIST:
      {
	key="filelist";   value=env->visualList;
	break;
      }
       case VF_SET_WRVOTABLEFORCE:
      {
	key="force"; value="";
	break;
      }    
       case VF_SET_HISTOGRAMBIN:
      {
	key="histogram";
	key="test";
	value=env->histobin;
	if(value=="-1")
	  value="";
	break;
      }    
      case VF_SET_VO:
      {
	key="VO";   value=env->vo;
	break;
      }
      case VF_SET_LFNOUT:
      {
	key="lfnout";   value=env->lfnout;
	break;
      }
      case VF_SET_SE:
      {
	key="se";   value=env->se;
	break;
      }
      case VF_SET_MRPOS:
      {
	key="pos";   value=env->MRpos;
	break;
      }
      case VF_SET_MRGEOMETRY:
      {
	key="geometry"; value=env->MRgeometry;
	break;
      }
      case VF_SET_MRBACKGROUND:
      {
	key="geometry"; value=env->MRbackground;
	break;
      }
      case VF_SET_DIMVOX:
      {
	key="dimvox"; value=env->dimvox;
	break;
      }
      case VF_SET_TRACKPLANEDIST:
      {
	key="trackplanedist"; value=env->trackplanedist;
	break;
	  }
	 case VF_SET_INNERDIST:
      {
	key="innerdist"; value=env->innerdist;
	break;
      }
      case VF_SET_OUTPOINTS:
      {
	key="outpoints"; value=env->outpoints;
	break;
      }
      case VF_SET_OUTVOL:
      {
	key="outvol"; value=env->outvol;
	break;
      }

    }//switch

    appParameters.insert(make_pair(key, value)); 
    
  }// if
}// for

// this part implements the multi processes with MPI. It is the map that contains the nfo for multiprocess    
std::stringstream MpiRank, MpiSize;
MpiRank<<rank;
MpiSize<<size;
appParameters.insert(make_pair("MpiRank",MpiRank.str()));
appParameters.insert(make_pair("MpiSize",MpiSize.str()));
///////    // CheckG appParameters

startFilter startFilter(appParameters);
return EXIT_SUCCESS;
}
} // extern "C"
//#ifdef HAVE_CONFIG_H
//#include <config.h>
//#endif
/***************************************************************************
 *   Copyright (C) 2008 by Gabriella Caniglia *
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
#include "visivodef.h"
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <vector>
#ifdef WIN32
	#include <direct.h>
#endif
#include "commandline.h"
#include "visivoutils.h"
#include "time.h"
#include <sstream>
#include <fstream>
#ifdef WIN32
	#include <io.h>
#else
	#include <unistd.h>
#endif

extern "C"{

// Dichiarazione delle funzioni principali.
int VI_Import(VisIVOImporter *env);
void* VI_Import_Thread(void *t_param);


//----------------------------
int VA_Import(VisIVOImporter *e, VisIVOAsynchId *id)
//----------------------------
{
	int iret;
	id->withThread=1;
	#ifdef WIN32
 		iret=VI_Import(e);
 		return iret;
	#endif
	#ifndef WIN32
	id->env=e;
	id->state=runningThread;
	iret=pthread_create(&(id->threadId), NULL, &VI_Import_Thread, (void *)id);
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
void *VI_Import_Thread(void *id)
//---------------------------
{
	int iret, i;
	VisIVOAsynchId *idTh;
	idTh=(VisIVOAsynchId *)id;
	iret=VI_Import((VisIVOImporter *)(idTh->env));
	idTh->errorCode=iret;
	idTh->state=successfulEndThread;
	return 0;
}

void VI_FreeImporter(VisIVOImporter* env) {
    if (!env) return;
    delete env->pComLine;  
}

//---------------------------
int VI_Import(VisIVOImporter *env)
//---------------------------
{
bool fitsRestore=false;
std::string fitsOriginalName;
std::string tmpOutName;  
env->pComLine=new CommandLine;
std::vector<std::string> args;

if(env->setatt[VI_SET_FFORMAT]==0)
{
  std::cerr<<"VI_Import: Invalid file format. Please set VI_SET_FFORMAT ";
  std::cerr<<"with VI_SettAtt function"<<std::endl;
  return invalidFFormatCode;
}
if(env->setatt[VI_SET_FILEPATH]==0)
{
  std::cerr<<"VI_Import: Invalid input file. Please set VI_SET_FILEPATH ";
  std::cerr<<"with VI_SettAtt function"<<std::endl;
  return invaldInputFile;
}
for(int idPar=0; idPar<NPAR; idPar++)
{
  if(env->setatt[idPar]==1)
  {
//   std::clog<<"VI_Import i="<<idPar<<" setatt "<<env->setatt[idPar]<<std::endl;

   switch(idPar)
   {
    case VI_SET_FFORMAT:
    {
      args.push_back("--fformat");
      args.push_back(env->fformat);
      break;
    }
    case VI_SET_FILEPATH:
    {
      args.push_back(env->infile);
      break;
    }
     case VI_SET_OUTFILEVBT:
    {
      args.push_back("--out");
      args.push_back(env->outfile);
      break;
    }
     case VI_SET_VOLUME:
    {
      args.push_back("--volume");

      std::stringstream sstmp;
      std::string stmp;
      sstmp << env->comp[0];
      sstmp>>stmp;
      args.push_back("--compx");
      args.push_back(stmp.c_str());

      sstmp.clear();
      sstmp << env->comp[1];
      sstmp>>stmp;
      args.push_back("--compy");
      args.push_back(stmp.c_str());

      sstmp.clear();
      sstmp << env->comp[2];
      sstmp>>stmp;
      args.push_back("--compz");
      args.push_back(stmp.c_str());
      
      sstmp.clear();
      sstmp << env->size[0];
      sstmp>>stmp;
      args.push_back("--sizex");
      args.push_back(stmp.c_str());

      sstmp.clear();
      sstmp << env->size[1];
      sstmp>>stmp;
      args.push_back("--sizey");
      args.push_back(stmp.c_str());

      sstmp.clear();
      sstmp << env->size[2];
      sstmp>>stmp;
      args.push_back("--sizez");
      args.push_back(stmp.c_str());

      break;
    }
     case VI_SET_USERPWD:
    {
      args.push_back("--userpwd");
      args.push_back(env->userpwd);
      break;
    }
    case VI_SET_BINARYHEADER:
    {
      args.push_back("--binaryheader");
      args.push_back(env->binaryheader);
      break;
    }
    case VI_SET_MISSINGVALUE:
    {
      std::stringstream sstmp;
      std::string stmp;
      sstmp << env->missing;
      sstmp>>stmp;
      args.push_back("--missingvalue");
      args.push_back(stmp.c_str());
      break;
    }
    case VI_SET_TEXTVALUE:
    {
      std::stringstream sstmp;
      std::string stmp;
      sstmp << env->text;
      sstmp>>stmp;
      args.push_back("--textvalue");
      args.push_back(stmp.c_str());
      break;
    }
    case VI_SET_BIGENDIAN:
    {
      args.push_back("--bigendian");
      break;
    }
    case VI_SET_DOUBLE:
    {
      args.push_back("--double");
      break;
    }
    case VI_SET_NPOINTS:
    {
      std::stringstream sstmp;
      std::string stmp;
      sstmp << env->npoints;
      sstmp>>stmp;
      args.push_back("--npoints");
      args.push_back(stmp.c_str());
      break;
    }
    case VI_SET_DATASETLIST:
    {
      args.push_back("--datasetlist");
      args.push_back(env->datasetList);
      break;
    }
    case VI_SET_HYPERSLAB:
    {
      args.push_back("--hyperslab");
      args.push_back(env->hyperslab);
      break;
    }
    case VI_SET_VO:
    {
      args.push_back("--VO");
      args.push_back(env->VO);
      break;
    }
    case VI_SET_LFNOUT:
    {
      args.push_back("--lfnout");
      args.push_back(env->lfnout);
      break;
    }
     case VI_SET_SE:
    {
      args.push_back("--se");
      args.push_back(env->se);
      break;
    }
     case VI_SET_INMEMORY:
    {
      args.push_back("--inmemory");
      break;
    }
  
   } //switch
  }//if
}
 

    for (int i=0;i<args.size()-1;i++) // known bug: bug correction on FitsTable with very long out name
    {
      std::string arg=args[i];
      std::string arg1=args[i+1];
      if((arg=="--fformat") && arg1=="fitstable" )
      {
	  for (i=0;i<args.size();i++) 
    	{
		std::string arg2=args[i];
		if(arg2=="--out")
		{
  			std::string tmp=args[i+1];
			std::string tmpOutDir=getDir(tmp);
			tmpOutName=getName(tmp);
			chdir(tmpOutDir.c_str());
			if(tmpOutName.size()>28)
			{
				fitsRestore=true;
				fitsOriginalName=tmpOutName;
				std::stringstream fileNameOutputSStream;
				time_t rawtime;
				struct tm *timeinfo;
				char buffer [80];
				time ( &rawtime );
				timeinfo = localtime ( &rawtime );
				strftime (buffer,80,"%Y%m%d%H%M%S",timeinfo);
  				fileNameOutputSStream<<buffer<<".bin"; 
				tmpOutName=fileNameOutputSStream.str();
				args[i+1]=tmpOutName;
				tmp=args[i+1];
				tmp=args[i+1];
			}
		}
	}
      }
    }


// for(int i=0;i<args.size();i++)std::clog<<args[i]<<" "; 
//  std::clog<<std::endl<<"test VisIVO Importer"<<std::endl;
  int ret=0;
  ret=env->pComLine->parseOption (args );
  if(ret<0) return invalidImporterOptions;
  ret=env->pComLine->loadFile();
  if(ret<0) return invalidImporterOperation;
 
  if((env->pComLine -> getRemoteFile() !="noremote") && !((env->pComLine -> getType() =="rawpoints") || (env->pComLine -> getType() =="xml") ) )
  { 
	std::string remoteFile;
	remoteFile=env->pComLine -> getRemoteFile();
	remove(remoteFile.c_str());
	if(env->pComLine -> getType() =="binary")
	{
		remoteFile=remoteFile.append(".head");
		remove(remoteFile.c_str());
	}
	
   }
  if(fitsRestore)
  {
	if(fitsOriginalName.find(".bin") == std::string::npos)
		fitsOriginalName=fitsOriginalName+".bin";
  	rename(tmpOutName.c_str(),fitsOriginalName.c_str());
	tmpOutName=tmpOutName+".head";
	fitsOriginalName=fitsOriginalName+".head";
  	rename(tmpOutName.c_str(),fitsOriginalName.c_str());
  }
//  std::cout<<"VisIVOImporter operation done."<<std::endl;		
  //if ( pComLine)
  //  delete pComLine ;
  if(env->enableInMemory){
    env->memTables = &(env->pComLine->getTable());
  }

  return 0;
}

} //extern "C"

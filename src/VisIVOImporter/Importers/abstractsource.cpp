/***************************************************************************
 *   Copyright (C) 2008 by Gabriella Caniglia, Marco Comparato *
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
#include <cstdlib>
#include <cstring>

#include "abstractsource.h"
#include "historyXmlWriter.h"
#include "visivoutils.h"

#include <iostream>
#include <fstream>
#include <sstream>
const unsigned int AbstractSource::MAX_LOAD=1000000;
const unsigned int AbstractSource::MAX_LARGE_LOAD=100000000;
const  int AbstractSource::MAX_INT=200000000;


//---------------------------------------------------------------------
AbstractSource::AbstractSource()
//---------------------------------------------------------------------
{
  m_pointsFileName = "";
  m_pointsBinaryName="";
  m_nRows = 0;
  m_nCols = 0;
  m_fieldNames.clear();
   
}

// //---------------------------------------------------------------------
// void AbstractSource::releaseResources()
// //---------------------------------------------------------------------
// {
//   if(m_visData.data)
//   {
//     for(int i = 0; i < m_nCols; i++)
//     {
//       if(m_visData.data[i])
//       {
//         delete [] m_visData.data[i];
//         m_visData.data[i] = NULL;
//       }
//     }
// 
//     delete [] m_visData.data;
//     m_visData.data = NULL;
//   }
// 
//   return;
// }

//---------------------------------------------------------------------
void AbstractSource::setPointsFileName(const char* fileName, const char* binaryName, 
				       const char* tableOrVolume, double size[], 
				       double comput[], const char* file, 
				       const char* endian, const char* type,
               const char* aliasParticle, const char* aliasHeader, 
				       long unsigned int points, 
				       const char* login, const char* binaryHeader, 
				       float missing, float text, 
				       std::string datasetdList,
				       std::vector<std::string> hyperslab, int fitshdunum, std::vector<std::string> fields)
//---------------------------------------------------------------------
{
  m_pointsFileName = fileName;
  m_pointsBinaryName=binaryName;
  m_file=tableOrVolume;
  m_login=login;
  m_binaryHeader=binaryHeader;

  m_volumeOrTable=file;
  m_cellSize[0]=size[0];
  m_cellSize[1]=size[1];
  m_cellSize[2]=size[2];
   
  m_cellComp[0]=comput[0];
  m_cellComp[1]=comput[1];
  m_cellComp[2]=comput[2];
  
  m_datasetList=datasetdList;
  m_hyperslab=hyperslab;
  m_fields=fields;
  
  m_nRows=points;
  m_type=type;
  m_endian=endian;
  m_aliasHeader=aliasHeader;
  m_aliasParticle=aliasParticle;
  MISSING_VALUE=missing;
  TEXT_VALUE=text;
    m_fitshdunum=fitshdunum;


  return;
}
//---------------------------------------------------------------------
void AbstractSource::setPointsFileName(const char *fileName,const char *binaryName)
//---------------------------------------------------------------------
{
  m_pointsFileName = fileName;
  m_pointsBinaryName=binaryName;
 
//   m_volumeOrTable=file;
//   m_cellSize[0]=size[0];
//   m_cellSize[1]=size[1];
//   m_cellSize[2]=size[2];
//    
//   m_cellComp[0]=comput[0];
//   m_cellComp[1]=comput[1];
//   m_cellComp[2]=comput[2];
//   
//   m_endian=endian;

  return;
}

//---------------------------------------------------------------------
void AbstractSource::setEndianism(const std::string& s) {  
//---------------------------------------------------------------------
  m_endian = s;  
}

//---------------------------------------------------------------------
int AbstractSource::readHeader()
//---------------------------------------------------------------------
{
  return 1;
}

//---------------------------------------------------------------------
int AbstractSource::readData()
//---------------------------------------------------------------------
{
  return 1;
}

//---------------------------------------------------------------------
int AbstractSource::writeHistory (const char* histFile,const char* format,const char* out,const char* tableOrVolume,double comput[],double size[],
                                   const char* login, const char* binaryHeader, float missing,float text, const char* endian,
                                   const char* type, long unsigned int points,
                                   const char* vo, const char* se, const char* lfnout,const char* inputFile)
//---------------------------------------------------------------------
{
    new HistoryXmlWriter(histFile,format,out,tableOrVolume,comput, size, login, binaryHeader, missing, text, endian, type, points, vo, se, lfnout, inputFile);
    return 1;
}


std::vector<VSTable*>& AbstractSource::getMemTables() {
    if (memTables.empty()) {
      std::cerr << "No available memTable" << std::endl;
    }
    return memTables;
}
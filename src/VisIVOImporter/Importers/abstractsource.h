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

#ifndef ABSTRACTSOURCE_H
#define ABSTRACTSOURCE_H

#include <string>
#include <vector>
#include <gtest/gtest.h>
#include <vstable.h>
#include <vstablemem.h>


class AbstractSource
{
   static const  int MAX_INT;

  public:
    AbstractSource();
    virtual ~AbstractSource() = default;

    void setEndianism(const std::string& s);
    void setPointsFileName(const char* fileName, const char* binaryName, 
				       const char* tableOrVolume, double size[], 
				       double comput[], const char* file, 
				       const char* endian, const char* type,
               const char* aliasParticle, const char* aliasHeader, 
				       long unsigned int points, 
				       const char* login, const char* binaryHeader, 
				       float missing, float text, 
				       std::string datasetdList,
				       std::vector<std::string> hyperslab, int fitshdunum, std::vector<std::string> fields);

    void setPointsFileName(const char *fileName,const char *binaryName);
//     void releaseResources();
    virtual int readHeader() = 0;
    virtual int readData() = 0;
    std::vector<VSTable*>& getMemTables();//{return memTables;}
    
    
    void setUseMem(bool use) { useMemory = use; }
    int writeHistory (const char* histFile,const char* format,const char* out,const char* tableOrVolume,double comput[],double size[],const char* login, const char* binaryHeader, float missing,float text, const char* endian,const char* type, long unsigned int points, const char* vo, const char* se, const char* lfnout,const char* inputFile);


  protected:
   static const unsigned int MAX_LOAD;
   static const unsigned int MAX_LARGE_LOAD;

   float MISSING_VALUE; //! a negative value used in case of missing data
   float TEXT_VALUE; //! a negative value used in case of ascii text
   std::string m_pointsFileName;
    std::string m_pointsBinaryName;
    unsigned long long int m_nRows;
    int m_nCols;

    std::vector<std::string> m_fieldNames;  //!column List
    std::string m_volumeOrTable;
    std::string m_type;
    std::string m_file;
    std::string m_endian;
    std::string m_aliasParticle;
    std::string m_aliasHeader;
    std::string m_login;
    std::string m_binaryHeader;
    std::string m_datasetList;
    double m_cellSize[3], m_cellComp[3];
    int maxInt(){return MAX_INT;};
    std::vector<std::string>  m_hyperslab;
    std::vector<std::string>  m_fields;
    std::vector<VSTable*> memTables;
    int m_fitshdunum;
    bool useMemory;
    private:
};

#endif

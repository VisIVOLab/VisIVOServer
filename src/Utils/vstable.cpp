/***************************************************************************
 *   Copyright (C) 2008 by Marco Comparato, Ugo Becciani                   *
 *   marco.comparato@oact.inaf.it, ugo.becciani@oact.inaf.it               *
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
#include "vstable.h"
// #include "VisIVOFiltersConfigure.h"
#include "visivoutils.h"
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>

const unsigned int VSTable::MAX_NUMBER_TO_SKIP = 1000000000;
//const unsigned int VSTable::MAX_NUMBER_ROW_REQUEST = 2147483647;
//const unsigned int VSTable::MAX_NUMBER_TO_SKIP = 250000000;
const unsigned int VSTable::MAX_NUMBER_ROW_REQUEST = 250000000;

VSTable::VSTable()
{
  m_locator = "";
  m_endiannes = "unknown";
  m_type      = "unknown";

  m_nCols = 0;
  m_nRows = 0;
  m_tableExist=false;

  m_isVolume = false;

  for(int i = 0; i < 3; ++i)
  {
    m_cellNumber[i] = 0;
    m_cellSize[i] = 0;
  }
}

VSTable::VSTable(std::string locator, std::string name /*= ""*/, std::string description /*= ""*/) //QUI modificata per letture header etc
  : VSObject(name, description)
{
  bool headerExist=false;

  setLocator(locator);
  m_endiannes = "unknown";
  m_type      = "unknown";
  m_tableExist=false;

  m_nCols = 0;
  m_nRows = 0;

  m_isVolume = false;

  for(int i = 0; i < 3; ++i)
  {
    m_cellNumber[i] = 0;
    m_cellSize[i] = 0;
  }
  
  headerExist=readHeader();
  if(!headerExist)
  {
    std::cerr<<"Invalid table header "<<m_locator<<std::endl;
    return;
  } else
  {
    std::ifstream inTable(m_locator.c_str());
    if(!inTable)
    {
      std::cerr<<"Invalid table "<<m_locator<<std::endl;
      return;
    }
    inTable.close();
    m_tableExist=true;
  }
}

VSTable::~VSTable()
{
}

bool VSTable::setEndiannes(std::string endiannes)
{
  if(endiannes == "big" || endiannes == "b")
  {
    m_endiannes = "big";

    return true;
  }
  else if(endiannes == "little" || endiannes == "l")
  {
    m_endiannes = "little";

    return true;
  }
     
  m_endiannes = "unknown";
  
  return false;
}

bool VSTable::setType(std::string type)
{
  if(type == "float" || type == "f")
  {
    m_type = "float";
    
    return true;
  }
  else if(type == "double" || type == "d")
  {
    m_type = "double";
    
    return true;
  }
    
  m_type = "unknown";

  return false;
}

bool VSTable::addCol(std::string name)
{
  if(name == "")
    return false;

  unsigned int size = m_colVector.size();

  for(int i = 0; i < size; ++i)
  {
    if(m_colVector[i] == name)
    {
      std::cerr<<"Invalid duplicate column name "<<name<<" in table."<<std::endl;
      return false;
    }  
  }

  m_colVector.push_back(name);
  m_nCols++;

  return true;
}

void VSTable::setCellSize(float xCellSize, float yCellSize, float zCellSize)
{
	m_cellSize[0] = xCellSize;
	m_cellSize[1] = yCellSize;
	m_cellSize[2] = zCellSize;

	return;
}

void VSTable::setLocator(std::string locator) 
{ 
  if(locator.find(".bin") != std::string::npos)
		m_locator=locator;
  else
		m_locator = locator + ".bin";


}   

void VSTable::setCellNumber(unsigned int xCellNumber, unsigned int yCellNumber, unsigned int zCellNumber)
{
	m_cellNumber[0] = xCellNumber;
	m_cellNumber[1] = yCellNumber;
	m_cellNumber[2] = zCellNumber;

	return;
}

bool VSTable::readHeader()
{
  int i = 0;

  std::string headerFileName = m_locator + ".head";
  std::ifstream inHeader(headerFileName.c_str());
  
  
   if(!inHeader)
     return false;

  //reads the number of cols and rows

  std::string dummy;
  std::string type;
  std::string endiannes;
 
  unsigned int nCols;
  
  inHeader >> type;
  if(!(type == "float" || type =="f"))
  {
	std::cerr<<"Invalid input not float table in  VisIVO Filters. Please use VisIVO Importer to convert it in float table "<<std::endl;
	return false; 
  } 
  inHeader >> nCols;

  std::string tmp = "";
  
  getline(inHeader, tmp); //to remove the carrige return character  (\r) from the last line
  getline(inHeader, tmp);

  std::stringstream sstmp(tmp);

  sstmp >> m_nRows;

  if(!(sstmp.eof()))
  {
    m_isVolume=true;

    sstmp >> m_cellNumber[0];
    sstmp >> m_cellNumber[1];
    sstmp >> m_cellNumber[2];
        
    sstmp >> m_cellSize[0];
    sstmp >> m_cellSize[1];
    sstmp >> m_cellSize[2];
  }
  
  inHeader >> endiannes;
  
  setType(type);
  setEndiannes(endiannes);

  m_colVector.clear();

  std::string name = "";

  for(i = 0; i < nCols; ++i)
  {
    inHeader >> name;
    addCol(name);
   }

  inHeader.close();
  if((m_colVector.size() != nCols )||(m_cellNumber[0] >0 && m_cellNumber[0] * m_cellNumber[1] * m_cellNumber[2] != m_nRows))
    {
    m_nCols = 0;
    m_colVector.clear();

    return false;
   }
  return true;
}

void VSTable::printSelf()
{
  VSObject::printSelf();

  std::clog << "Locator: "     << m_locator   << std::endl;
  std::clog << "Columns: "     << m_nCols     << std::endl;
  std::clog << "Rows: "        << m_nRows     << std::endl;
  std::clog << "Type: "        << m_type      << std::endl;
  std::clog << "Endiannes: "   << m_endiannes << std::endl;

  std::clog << "Column names:" << std::endl;

  unsigned int size = m_colVector.size();

  for(int i = 0; i < size; ++i)
    std::clog << m_colVector[i] << std::endl;

  return;
}

std::string VSTable::getColName(unsigned int i)
{
  if(i < m_colVector.size())
  {
    return m_colVector[i];
  }
  
  return "";
}

void VSTable::setColName(unsigned int i, std::string newColName)
{
  if(i < m_colVector.size())
  {
   	m_colVector[i]=newColName;
  }
  
  return;
}


int VSTable::getColId(std::string name)
{
  m_colVector.size();
  unsigned int size = m_colVector.size();
  for(unsigned int i = 0; i < size; ++i)
  {
    if(m_colVector[i] == name)
      return i;
  }

  return -1;
}

bool VSTable::createTable()  
{
  if(m_tableExist)
  {
    std::cerr<<"Table "<<m_locator<<" already exist"<<std::endl;
    return false;
  }
  std::string fileName = m_locator;
  std::ifstream checkTable(fileName.c_str());
  if(checkTable)
  {
    std::cerr<<"Cannot create table. File "<<m_locator<<" already exist"<<std::endl;
    return false;
  }
  std::ofstream outTable(fileName.c_str());
  if(!outTable)
    return false;
	
  std::ios::off_type indexOffset;
  unsigned long long int totNumberToSkip;
  int nSkip;
  float zero=0;

  totNumberToSkip = m_nRows*m_nCols-1;
  while(totNumberToSkip!=0)
  {
    nSkip=MAX_NUMBER_TO_SKIP;
    if(nSkip > totNumberToSkip) nSkip=(int) totNumberToSkip;
    indexOffset=(std::ios::off_type) nSkip*sizeof(float);
    outTable.seekp(indexOffset,std::ios::cur);
    totNumberToSkip=totNumberToSkip-nSkip;
  }
  outTable.write((char *) &zero,1*sizeof(float));
  outTable.close();
  m_tableExist=true;
  return true;
}

bool VSTable::createTable(float **fArray)  
{
  if(m_tableExist)
  {
    std::cerr<<"Table "<<m_locator<<" already exist"<<std::endl;
    return false;
  }
  std::string fileName = m_locator;
  std::ifstream checkTable(fileName.c_str());
  if(checkTable)
  {
    std::cerr<<"Cannot create table. File "<<m_locator<<" already exist"<<std::endl;
    return false;
  }
  std::ofstream outTable(fileName.c_str());
  if(!outTable)
    return false;
	
	
  for(int j=0;j<m_nCols;j++)
  {
    std::ios::off_type indexOffset;
    unsigned long long int totNumberToWrite;
    int nWrite;
    float *fPointer=&fArray[j][0];
    totNumberToWrite = m_nRows;
    while(totNumberToWrite!=0)
    {
      nWrite=MAX_NUMBER_TO_SKIP;
      if(nWrite > totNumberToWrite) nWrite=(int) totNumberToWrite;
      outTable.write((char *) fPointer,nWrite*sizeof(float));
      fPointer=fPointer+nWrite;
      totNumberToWrite=totNumberToWrite-nWrite;
    }
			
  }
  outTable.close();
  return true;
  m_tableExist=true;
	
}


bool VSTable::writeTable()
{
  if(m_tableExist)
  {
    std::cerr<<"Table "<<m_locator<<" already exist"<<std::endl;
    return false;
  }

  bool wrHeader=writeHeader();
  if(!wrHeader)
  {
    std::cerr<<"No Table header was written: "<<m_locator<<".head"<<std::endl;
    return false;
  }
    return false;
}

bool VSTable::writeTable(float **fArray)
{
  if(m_tableExist)
  {
    std::cerr<<"Table "<<m_locator<<" already exist"<<std::endl;
    return false;
  }

  bool wrHeader=writeHeader();
  if(!wrHeader)
  {
    std::cerr<<"No Table header was written: "<<m_locator<<".head"<<std::endl;
    return false;
  }
  bool wrTable=createTable(fArray);
  if(!wrTable)
  {
    std::cerr<<"No Table  was written: "<<m_locator<<".head"<<std::endl;
    return false;
  }
  return true;
}


bool VSTable::writeHeader()
{
  int i = 0;

  std::string headerFileName = m_locator + ".head";
 
  std::ofstream outHeader(headerFileName.c_str());

  if(!outHeader)
    return false;

  outHeader << m_type << std::endl;
  outHeader << m_nCols << std::endl;
  outHeader << m_nRows;
  
  if(m_isVolume)
  {
    for(int i = 0; i < 3; ++i)
      outHeader << " " << m_cellNumber[i];
    
    for(int i = 0; i < 3; ++i)
      outHeader << " " << m_cellSize[i];
  }

  outHeader << std::endl;

  outHeader << m_endiannes << std::endl;

  unsigned int size = m_colVector.size();

  for(int i = 0; i < size; ++i)
    outHeader << m_colVector[i] << std::endl;

  outHeader.close();

//  createTable(); //QUI introdotto MA occorre inserire verifca che lheader non esiste gia'

  return true;
}

//-------------------
int VSTable::getColumn(int colNumber,float *Col, int fromRow, int toRow)
//-------------------
{
  if(fromRow<=-1) fromRow=0;
  if(toRow==-1) toRow=(int) m_nRows-1;
  if(toRow<0)
  {
    std::cerr<<"Invalid toRow parameter:"<<toRow<<std::endl;
    return -4;
  } 
  if((toRow-fromRow+1) > MAX_NUMBER_ROW_REQUEST*4)
  {
    std::cerr<<"Cannot be requested more than "<< MAX_NUMBER_ROW_REQUEST*4 <<" row lines"<<std::endl;
    return -3;
  }
  std::ifstream fileInput(m_locator.c_str(), std::ios::in | std::ios::binary);
  if(!fileInput)
  {
    std::cerr<<"Cannot open binary file"<< m_locator <<std::endl;
    return -2;
  }

  std::ios::off_type indexOffset;
  unsigned long long int totNumberToSkip;
  int nSkip;
  int nLoad;

    if(colNumber > m_nCols-1)
    {
      std::cerr << "Invalid Column Id: "<< colNumber<<std::endl;
      return -4;
    }
    if(toRow > (int)( m_nRows-1))
    {
      std::cerr << "Warning: Invalid request of toRow:" << toRow <<" lowered to "<<m_nRows -1<< std::endl;
      toRow=(int) (m_nRows-1);
    }
    nLoad = (int) (toRow-fromRow+1);
    if(nLoad < 0)
    {
      std::cerr << "Error: Invalid range fromRow: " << fromRow <<" toRow: " << toRow << std::endl;
      return -3;
    }

    totNumberToSkip = colNumber*m_nRows+fromRow;

    while(totNumberToSkip!=0)
    {
      nSkip=MAX_NUMBER_TO_SKIP;
      if(nSkip > totNumberToSkip) nSkip=(int) totNumberToSkip;
      indexOffset=(std::ios::off_type) nSkip*sizeof(float);
      fileInput.seekg(indexOffset,std::ios::cur);
      totNumberToSkip=totNumberToSkip-nSkip;
    }
    fileInput.read((char *) Col,nLoad*sizeof(float));
  
    fileInput.close();
#ifdef VSBIGENDIAN
  std::string endianism="big";
#else	
  std::string endianism="little";
#endif

  if((endianism=="big" && m_endiannes=="little") || (endianism=="little" && m_endiannes=="big"))
  {
	std::cerr<<"Warning: endianism swap executed on table "<<m_locator<<std::endl;
    	for(int k=0;k<nLoad;k++)
	 Col[k]=floatSwap((char *)(&Col[k]));
  }	 

  return nLoad;
}


//--------------------
int VSTable::getColumn(unsigned int *colList, unsigned int nOfCol, unsigned long long int fromRow, unsigned long long int toRow, float **fArray)
// Get all the values of the table columns listed in colList, fromRow toRow
// i list array of columns of nOfCol elements. 0 Is the first column 
// from -to. fArray must be allocated by the calling program and must contain all requested elements.Extreme values fromRow  and to Row are included. fromRom=0 toRow=10 are 11 elements 0,1,..10. 
// the method return the number of elements in each column j of fArray[j][i] (i= to-from)
// colList from 0 to m_nCols -1 
{
  if((toRow-fromRow+1) > MAX_NUMBER_ROW_REQUEST)
  {
    std::cerr<<"Cannot be requested more than "<< MAX_NUMBER_ROW_REQUEST <<" row lines"<<std::endl;
    return -3;
  }

  std::ifstream fileInput(m_locator.c_str(), std::ios::in | std::ios::binary);
  if(!fileInput)
  {
    std::cerr<<"Cannot open binary file"<< m_locator <<std::endl;
    return -2;
  }

  std::ios::off_type indexOffset;
  unsigned long long int totNumberToSkip;
  int nSkip;
  unsigned int column;
  int nLoad;
	
  for(unsigned int j=0; j<nOfCol;j++)
  {
    column=colList[j];
    if(column > m_nCols-1)
    {
      std::cerr << "Invalid Column Id"<< std::endl;
      continue;
    }
    if(toRow > m_nRows-1)
    {
      std::cerr << "Warning: Invalid request of toRow:" << toRow <<" lowered to "<<m_nRows -1<< std::endl;
      toRow=m_nRows-1;
    }
    nLoad = (int) (toRow-fromRow+1);
    if(nLoad < 0)
    {
      std::cerr << "Error: Invalid range fromRow: " << fromRow <<" toRow: " << toRow << std::endl;
      return -3;
    }

    if (j==0) totNumberToSkip = colList[j]*m_nRows+fromRow;
    else
    {
      if(colList[j]>=colList[j-1])
        totNumberToSkip = (colList[j]-colList[j-1]-1)*m_nRows+fromRow+(m_nRows-toRow)-1;
      else
      {
        indexOffset=0;
        fileInput.seekg(indexOffset,std::ios::beg);
        totNumberToSkip = colList[j]*m_nRows+fromRow;
			
      }
    }
    while(totNumberToSkip!=0)
    {
      nSkip=MAX_NUMBER_TO_SKIP;
      if(nSkip > totNumberToSkip) nSkip=(int) totNumberToSkip;
      indexOffset=(std::ios::off_type) nSkip*sizeof(float);
      fileInput.seekg(indexOffset,std::ios::cur);
      totNumberToSkip=totNumberToSkip-nSkip;
    }
    fileInput.read((char *) &fArray[j][0],nLoad*sizeof(float));
  }
  fileInput.close();
#ifdef VSBIGENDIAN
  std::string endianism="big";
#else	
  std::string endianism="little";
#endif

  if((endianism=="big" && m_endiannes=="little") || (endianism=="little" && m_endiannes=="big"))
  {
	std::cerr<<"Warning: endianism swap executed on table "<<m_locator<<std::endl;
  	for(unsigned int j=0; j<nOfCol;j++)
    		for(int k=0;k<nLoad;k++)
			fArray[j][k]=floatSwap((char *)(&fArray[j][k]));
  }	 

  return nLoad;
}


int VSTable::getColumnList(unsigned int *colList,unsigned int nOfCol, unsigned long long int *list, int nOfEle, float **fArray)
// Get only the values listed in list of the table columns listed in colList
// colList list array of columns of nOfCol elements. 0 Is the first column 
//  list is an nOfEle array elements containing requested number of rows
// the method return the number of elements in each column j of fArray[j][i] (i= 0-nOfEle)

{

  if(nOfEle > MAX_NUMBER_ROW_REQUEST)
  {
    std::cerr<<"Cannot be requested more than "<< MAX_NUMBER_ROW_REQUEST <<" row lines"<<std::endl;
    return -3;
  }

  std::ifstream fileInput(m_locator.c_str(), std::ios::in | std::ios::binary);
  if(!fileInput)
  {
    std::cerr<<"Cannot open binary file"<<m_locator<<std::endl;
    return -2;
  }
  unsigned long long int indexLast=0,stepToSkip;
  std::ios::off_type indexOffset;
  unsigned long long int totNumberToSkip;
  int nSkip;
  unsigned int column;
  int nLoad=1;
  int totLoad=0;

	
  for(unsigned int k=0; k<nOfCol;k++)
  {
    column=colList[k];
    if(column > m_nCols-1)
    {
      std::cerr << "Invalid Column Id"<< std::endl;
      continue;
    }

// going at the beginning of column
    if(k==0)
    {
      totNumberToSkip=colList[k]*m_nRows;
      while(totNumberToSkip!=0)
      {
        nSkip=MAX_NUMBER_TO_SKIP;
        if(nSkip > totNumberToSkip) nSkip=(int) totNumberToSkip;
        indexOffset=(std::ios::off_type) nSkip*sizeof(float);
        fileInput.seekg(indexOffset,std::ios::cur);
        totNumberToSkip=totNumberToSkip-nSkip;
      }
    }
    if(k>0 && colList[k]<colList[k-1])
    {
      indexOffset=0;
      fileInput.seekg(indexOffset,std::ios::beg);

      totNumberToSkip=colList[k]*m_nRows;
      while(totNumberToSkip!=0)
      {
        nSkip=MAX_NUMBER_TO_SKIP;
        if(nSkip > totNumberToSkip) nSkip=(int) totNumberToSkip;
        indexOffset=(std::ios::off_type) nSkip*sizeof(float);
        fileInput.seekg(indexOffset,std::ios::cur);
        totNumberToSkip=totNumberToSkip-nSkip;
      }
    }
    if(k>0 && colList[k]>=colList[k-1])
    {
      totNumberToSkip=(colList[k]-colList[k-1]-1)*m_nRows + m_nRows- indexLast-1;
      while(totNumberToSkip!=0)
      {
        nSkip=MAX_NUMBER_TO_SKIP;
        if(nSkip > totNumberToSkip) nSkip=(int) totNumberToSkip;
        indexOffset=(std::ios::off_type) nSkip*sizeof(float);
        fileInput.seekg(indexOffset,std::ios::cur);
        totNumberToSkip=totNumberToSkip-nSkip;
      }
    }
    //

    for(int j=0;j<nOfEle;j++)
    {
      if(list[j] > m_nRows-1)
      {
        std::cerr << "Warning: Invalid list request:" << list[j] << std::endl;
        continue;
      }
      indexLast=list[j];
      if(j==0)
      {
        totNumberToSkip=list[j];
        stepToSkip=0;
      }
      if(j>0)
      {
        stepToSkip=(list[j]-list[j-1]-1);
        if(list[j-1]>list[j]-1)
        {
          totNumberToSkip=colList[k]*m_nRows+list[j];
          indexOffset=0;
          fileInput.seekg(indexOffset,std::ios::beg);
        } else
          totNumberToSkip=stepToSkip;
      }
      while(totNumberToSkip!=0)
      {
        nSkip=MAX_NUMBER_TO_SKIP;
        if(nSkip > totNumberToSkip) nSkip=(int) totNumberToSkip;
        indexOffset=(std::ios::off_type) nSkip*sizeof(float);
        fileInput.seekg(indexOffset,std::ios::cur);
        totNumberToSkip=totNumberToSkip-nSkip;
      }
      fileInput.read((char *) &fArray[k][j],nLoad*sizeof(float));
      totLoad++;
    }
	
  }

#ifdef VSBIGENDIAN
  std::string endianism="big";
#else	
  std::string endianism="little";
#endif

  if((endianism=="big" && m_endiannes=="little") || (endianism=="little" && m_endiannes=="big"))
  {
	std::cerr<<"Warning: endianism swap executed on table "<<m_locator<<std::endl;
  	for(unsigned int k=0; k<nOfCol;k++)
    		for(int j=0;j<nOfEle;j++)
			fArray[k][j]=floatSwap((char *)(&fArray[k][j]));
  }	 
  return totLoad;
}

int VSTable::putColumn(unsigned int *colList, unsigned int nOfCol, unsigned long long int fromRow, unsigned long long int toRow, float **fArray)
// i list array of columns of nOfCol elements. 0 Is the first column 
// from -to. fArray must be allocated by the calling program and must contain all requested elements. Extreme values fromRow  and to Row are included. fromRom=0 toRow=10 are 11 elements 0,1,..10. 
// the method return the number of elements in each column j of fArray[j][i] (i= to-from)
// colList from 0 to m_nCols -1 
{
  if((toRow-fromRow+1) > MAX_NUMBER_ROW_REQUEST)
  {
    std::cerr<<"Cannot be put more than "<< MAX_NUMBER_ROW_REQUEST <<" row lines"<<std::endl;
    return -3;
  }

  std::ifstream fileInput(m_locator.c_str(), std::ios::in | std::ios::binary);
  if(!fileInput)
  {
    std::ofstream fileOutput(m_locator.c_str(), std::ios::out | std::ios::binary);
    fileOutput.close();
  } else
    fileInput.close();

    std::fstream fileOutput(m_locator.c_str(), std::ios::in |  std::ios::out | std::ios::binary);
	
    if(!fileOutput)
    {
      std::cerr<<"Cannot open binary file"<< m_locator <<std::endl;
      return -2;
    }
    unsigned long long int indexSeek;
    std::ios::off_type indexOffset;
    unsigned long long int totNumberToSkip;
    int nSkip;
    unsigned int column;
    int nLoad;
	
    for(unsigned int j=0; j<nOfCol;j++)
    {
      column=colList[j];
      if(column > m_nCols+1)
      {
        std::cerr << "Invalid Column Id"<< std::endl;
        continue;
      }
      if(toRow > m_nRows-1)
      {
        std::cerr << "Warning: Invalid request of toRow:" << toRow <<" lowered to "<<m_nRows -1<< std::endl;
        toRow=m_nRows-1;
      }
      nLoad = (int) (toRow-fromRow+1);
      if(nLoad <= 0)
      {
        std::cerr << "Error: Invalid range fromRow: " << fromRow <<" toRow: " << toRow << std::endl;
        return -3;
      }

      if (j==0) totNumberToSkip = colList[j]*m_nRows+fromRow;
      else
      {
        if(colList[j]>=colList[j-1])
          totNumberToSkip = (colList[j]-colList[j-1]-1)*m_nRows+fromRow+(m_nRows-toRow)-1;
        else
        {
          indexOffset=0;
          fileOutput.seekp(indexOffset,std::ios::beg);
          totNumberToSkip = colList[j]*m_nRows+fromRow;
			
        }
      }
      while(totNumberToSkip!=0)
      {
        nSkip=MAX_NUMBER_TO_SKIP;
        if(nSkip > totNumberToSkip) nSkip=(int) totNumberToSkip;
        indexOffset=(std::ios::off_type) nSkip*sizeof(float);
        fileOutput.seekp(indexOffset,std::ios::cur);
        totNumberToSkip=totNumberToSkip-nSkip;
      }
      fileOutput.write((char *) &fArray[j][0],nLoad*sizeof(float));
    }
    fileOutput.close();
    return nLoad;
}


int VSTable::putColumnList(unsigned int *colList, unsigned int nOfCol, unsigned long long int *list, int nOfEle, float **fArray)
// colList list array of columns of nOfCol elements. 0 Is the first column 
//  list is an nOfEle array elements containing requested number of rows
// the method return the number of elements in each column j of fArray[j][i] (i= 0-nOfEle)

{
  if(nOfEle > MAX_NUMBER_ROW_REQUEST)
  {
    std::cerr<<"Cannot be put more than "<< MAX_NUMBER_ROW_REQUEST <<" row lines"<<std::endl;
    return -3;
  }

  std::fstream fileOutput(m_locator.c_str(), std::ios::in |  std::ios::out | std::ios::binary);
  if(!fileOutput)
  {
    std::cerr<<"Cannot open binary file"<<m_locator<<std::endl;
    return -2;
  }
  unsigned long long int indexSeek, indexLast=0,stepToSkip;
  std::ios::off_type indexOffset;
  unsigned long long int totNumberToSkip;
  int nSkip;
  unsigned int column;
  int nLoad=1;
  int totLoad=0;

	
  for(unsigned int k=0; k<nOfCol;k++)
  {
    column=colList[k];
    if(column > m_nCols-1)
    {
      std::cerr << "Invalid Column Id"<< std::endl;
      continue;
    }

// going at the beginning of column
    if(k==0)
    {
      totNumberToSkip=colList[k]*m_nRows;
      while(totNumberToSkip!=0)
      {
        nSkip=MAX_NUMBER_TO_SKIP;
        if(nSkip > totNumberToSkip) nSkip=(int) totNumberToSkip;
        indexOffset=(std::ios::off_type) nSkip*sizeof(float);
        fileOutput.seekp(indexOffset,std::ios::cur);
        totNumberToSkip=totNumberToSkip-nSkip;
      }
    }
    if(k>0 && colList[k]<colList[k-1])
    {
      indexOffset=0;
      fileOutput.seekp(indexOffset,std::ios::beg);

      totNumberToSkip=colList[k]*m_nRows;
      while(totNumberToSkip!=0)
      {
        nSkip=MAX_NUMBER_TO_SKIP;
        if(nSkip > totNumberToSkip) nSkip=(int) totNumberToSkip;
        indexOffset=(std::ios::off_type) nSkip*sizeof(float);
        fileOutput.seekp(indexOffset,std::ios::cur);
        totNumberToSkip=totNumberToSkip-nSkip;
      }
    }
    if(k>0 && colList[k]>=colList[k-1])
    {
      totNumberToSkip=(colList[k]-colList[k-1]-1)*m_nRows + m_nRows- indexLast-1;
      while(totNumberToSkip!=0)
      {
        nSkip=MAX_NUMBER_TO_SKIP;
        if(nSkip > totNumberToSkip) nSkip=(int) totNumberToSkip;
        indexOffset=(std::ios::off_type) nSkip*sizeof(float);
        fileOutput.seekp(indexOffset,std::ios::cur);
        totNumberToSkip=totNumberToSkip-nSkip;
      }
    }
    //

    for(int j=0;j<nOfEle;j++)
    {
      if(list[j] > m_nRows-1)
      {
        std::cerr << "Warning: Invalid list request:" << list[j] << std::endl;
        continue;
      }
      indexLast=list[j];
      if(j==0)
      {
        totNumberToSkip=list[j];
        stepToSkip=0;
      }
      if(j>0)
      {
        stepToSkip=(list[j]-list[j-1]-1);
        if(list[j-1]>list[j]-1)
        {
          totNumberToSkip=colList[k]*m_nRows+list[j];
          indexOffset=0;
          fileOutput.seekp(indexOffset,std::ios::beg);
        } else
          totNumberToSkip=stepToSkip;
      }
      while(totNumberToSkip!=0)
      {
        nSkip=MAX_NUMBER_TO_SKIP;
        if(nSkip > totNumberToSkip) nSkip=(int) totNumberToSkip;
        indexOffset=(std::ios::off_type) nSkip*sizeof(float);
        fileOutput.seekp(indexOffset,std::ios::cur);
        totNumberToSkip=totNumberToSkip-nSkip;
      }
      fileOutput.write((char *) &fArray[k][j],nLoad*sizeof(float));
      totLoad++;
    }
	
  }
  return totLoad;
}

bool VSTable::addCol(std::vector<std::string> listOfCol)
{
  unsigned int newSize = listOfCol.size();
  if(newSize == 0)
    return false;

  for (const auto& str : listOfCol) {
    if (std::find(m_colVector.begin(), m_colVector.end(), str) == m_colVector.end()) {
      m_colVector.push_back(str);
      m_nCols++;
    }
  }

  return true;
}
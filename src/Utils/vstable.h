/***************************************************************************
 *   Copyright (C) 2008 by Marco Comparato try it                          *
 *   marco.comparato@oact.inaf.it                                          *
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
#ifndef VSTABLE_H
#define VSTABLE_H

/**
	@author Marco Comparato <marco.comparato@oact.inaf.it>
        @author Ugo Becciani <ugo.becciani@oact.inaf.it>
*/

#include "vsobject.h"
#include <string>
#include <vector>
#include <algorithm> 

class VSTable : public VSObject
{
  protected:
  std::string m_locator; //! table path
  std::string m_endiannes;  //! endianism
  std::string m_type; //! data format (float, double ecc..)

  std::vector<std::string> m_colVector; //! vector containing columns names 

  static const unsigned int MAX_NUMBER_TO_SKIP; //! maximum number to skip in a single skip
  static const unsigned int MAX_NUMBER_ROW_REQUEST; //! maximum number of rows for each request
  static const unsigned int MAX_NUMBER_OF_BYTES; //!maximum number of bytes for each request

  unsigned int m_nCols;  //! number of columns
  unsigned long long int m_nRows; //! number of Rows
  bool m_tableExist; //! flag that is true if tabnle already exist (written on disk)

  bool m_isVolume;  //! true if table represent a volumetrivìc data
  unsigned int m_cellNumber[3];  //! number of cell of the volume (mesh)
  float m_cellSize[3]; //!cell geometry

  std::vector<std::vector<float>> fieldArray; //! vector containing the columns content


public:
  VSTable();
  VSTable(std::string locator, std::string name = "", std::string description = "");

  ~VSTable();

  void printSelf();  //!  write in stdlog the header of the table

  std::string getLocator() {return m_locator;};  //! retrn the path
  
  unsigned int getNumberOfColumns() { return m_nCols;}; 
  unsigned long long int getNumberOfRows() { return m_nRows;};
  
  void setIsVolume(bool isVolume) { m_isVolume = isVolume; };
  bool getIsVolume() { return m_isVolume; };
  const unsigned int* getCellNumber() {return m_cellNumber;};
  const float* getCellSize() {return m_cellSize;};

  std::string getEndiannes() { return m_endiannes;};
  std::string getType() { return m_type;};

  bool setEndiannes(std::string endiannes);
  bool setType(std::string type);//! set the table (float or double)
  void setLocator(std::string locator);  //! set the table path (used in writeheader and readheader)
  void setNumberOfRows(unsigned long long int nRows) {m_nRows = nRows;};
  void setCellNumber(unsigned int xCellNumber, unsigned int yCellNumber, unsigned int zCellNumber);
  void setCellSize(float xCellSize, float yCellSize, float zCellSize);
  //void SetNumberOfCols(unsigned int nCols) {m_nCols = nCols;};
  virtual bool addCol(std::string name); //! add an element (a column name) to the vector m_colVector
  virtual bool addCol(std::vector<std::string> listOfCol);
  void setColName(unsigned int i, std::string newColName); //! change the column name in the i position. The first column is 0.

  std::string getColName(unsigned int i); //! return the column name in the i position
  int getColId(std::string name); //! return the position of the column with specified name. The first column is 0.

  bool writeHeader(); //! write a new table header with the given member variables

  virtual int getColumn(unsigned int *colList, unsigned int nOfCol, unsigned long long int fromRow, unsigned long long int toRow, float **fArray); //! colList is a vector of nOfCol Elements, listing the columns ID I want to read. fArray is a 2D allocated array fArray[nOfCol,(toRow-fromRow+1)]: the method fill fArray with the column listed in colList (first column is 0: colList contains data in the range  0 - m_nCols -1) and with the postion specified in fromRow toRow. fromRow starts from 0.
  virtual int getColumn(int colNumber,float *Col, int fromRow=-1, int toRow=-1); //! get only one column. colNumerb is the columns ID. Col is an allocated pointer Col[toRow-fromRow+1]: the method fill Col with the column ID (first column is 0: in the range  0 - m_nCols -1) and with the postion specified in fromRow toRow-1. fromRow starts from 0. Default value -1 will read all column
  virtual int getColumnList(unsigned int *colList, unsigned  int nOfCol, unsigned long long int *list, int nOfEle, float **fArray);//! colList is a vector of nOfCol Elements, listing the columns ID I want to read. First column is 0. list is a vectorn of nOfEle elements I want to read from the table. fArray is a 2D allocated array fArray[nOfCol,nOfEle]: the method fill fArray with the column listed in colList (first column is 0: colList contains data in the range  0 - m_nCols -1) and with elements listed in list
  virtual int putColumn(unsigned int *colList, unsigned int nOfCol, unsigned long long int fromRow, unsigned long long int toRow, float **fArray);  //!  put data (see getColumn)
  virtual int putColumnList(unsigned int *colList,unsigned int nOfCol, unsigned long long int *list, int nOfEle, float **fArray);//! put data (see getColumnList)

  virtual bool writeTable(float **fArray);  //!write header and table (fArray must contain all data table)
  virtual bool createTable(float **fArray); //!write data (header already exist) in the table
  virtual bool writeTable();  //!write only header
  virtual bool createTable(); //!write a zero filled table
  virtual bool tableExist() { return m_tableExist;}; //! table already exist 
  bool readHeader(); //! fuction that read  the header table (filling the above values)

};

#endif
 
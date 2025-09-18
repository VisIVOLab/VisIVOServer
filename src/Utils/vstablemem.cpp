#include "vstablemem.h"
#include <iostream>
#include <cstring>

VSTableMem::VSTableMem()
{
}

VSTableMem::~VSTableMem()
{
    fieldArray.clear();
}

bool VSTableMem::createTable()
{
    try {
        fieldArray.assign(m_nCols, std::vector<float>(m_nRows, 0.0f));
    }
    catch (const std::bad_alloc& e) {
        std::cerr << "Memory allocation failed: " << e.what() << std::endl;
        return false;
    }
    m_tableExist = true;
    return true;
}

bool VSTableMem::createTable(float **fArray)
{
    try {
        fieldArray.assign(m_nCols, std::vector<float>(m_nRows));
        for (int j = 0; j < m_nCols; ++j) {
            std::copy(fArray[j], fArray[j] + m_nRows, fieldArray[j].begin());
        }
    }
    catch (const std::bad_alloc& e) {
        std::cerr << "Memory allocation failed: " << e.what() << std::endl;
        return false;
    }
    m_tableExist = true;
    return true;
}

int VSTableMem::getColumn(int colNumber, float *Col, int fromRow, int toRow)
{
    if (fromRow < 0) fromRow = 0;
    if (toRow < 0 || toRow >= (int)m_nRows) toRow = (int)m_nRows - 1;

    int nLoad = toRow - fromRow + 1;
    if (nLoad <= 0) {
        std::cerr << "Invalid range" << std::endl;
        return -1;
    }

    std::copy( fieldArray[colNumber].begin() + fromRow, fieldArray[colNumber].begin() + fromRow + nLoad, Col);
    return nLoad;
}

int VSTableMem::getColumn(unsigned int *colList, unsigned int nOfCol, unsigned long long int fromRow, unsigned long long int toRow, float **fArray)
{
    int nLoad = (int)(toRow - fromRow + 1);
    if (nLoad <= 0) {
        std::cerr << "Invalid range" << std::endl;
        return -1;
    }

    for (unsigned int j = 0; j < nOfCol; j++) {
        unsigned int column = colList[j];
        if (column >= m_nCols) {
            std::cerr << "Invalid Column Id" << std::endl;
            continue;
        }
        std::memcpy(fArray[j], &fieldArray[column][fromRow], nLoad * sizeof(float));
    }

    return nLoad;
}

int VSTableMem::getColumnList(unsigned int *colList, unsigned int nOfCol, unsigned long long int *list, int nOfEle, float **fArray)
{
    int totLoad = 0;
    for (unsigned int k = 0; k < nOfCol; ++k) {
        unsigned int col = colList[k];
        if (col >= m_nCols) {
            std::cerr << "Invalid Column Id" << std::endl;
            continue;
        }

        for (int j = 0; j < nOfEle; ++j) {
            if (list[j] >= m_nRows) {
                std::cerr << "Warning: Invalid list request: " << list[j] << std::endl;
                continue;
            }
            fArray[k][j] = fieldArray[col][list[j]];
            ++totLoad;
        }
    }
    return totLoad;
}

int VSTableMem::putColumn(unsigned int *colList, unsigned int nOfCol, unsigned long long int fromRow, unsigned long long int toRow, float **fArray)
{
    unsigned long long int nLoad = (toRow - fromRow + 1);
    if (nLoad <= 0) {
        std::cerr << "Invalid range" << std::endl;
        return -1;
    }

    for (unsigned int j = 0; j < nOfCol; j++) {
        unsigned int column = colList[j];
        if (column >= m_nCols) {
            std::cerr << "Invalid Column ID" << std::endl;
            continue;
        }

        std::copy(
            fArray[j], 
            fArray[j] + nLoad, 
            &fieldArray[column][fromRow]
        );

    }
    
    return nLoad;
}

int VSTableMem::putColumnList(unsigned int *colList, unsigned int nOfCol, unsigned long long int *list, int nOfEle, float **fArray)
{
    int totLoad = 0;
    for (unsigned int k = 0; k < nOfCol; ++k) {
        unsigned int col = colList[k];
        if (col >= m_nCols) {
            std::cerr << "Invalid Column Id" << std::endl;
            continue;
        }

        for (int j = 0; j < nOfEle; ++j) {
            if (list[j] >= m_nRows) {
                std::cerr << "Warning: Invalid list request: " << list[j] << std::endl;
                continue;
            }
            fieldArray[col][list[j]] = fArray[k][j];
            ++totLoad;
        }
    }
    return totLoad;
}

bool VSTableMem::addCol(std::vector<std::string> listOfCol)
{
  unsigned int newSize = listOfCol.size();
  if(newSize == 0)
    return false;

  for (const auto& str : listOfCol) {
    if (std::find(m_colVector.begin(), m_colVector.end(), str) == m_colVector.end()) {
      m_colVector.push_back(str);
      m_nCols++;
      fieldArray.push_back(std::vector<float>(m_nRows, 0.0f));
    }
  }

  return true;
}

bool VSTableMem::addCol(std::string name)
{
  if(!VSTable::addCol(name)){
    return false;
  }
  fieldArray.push_back(std::vector<float>(m_nRows, 0.0f));

  return true;
}
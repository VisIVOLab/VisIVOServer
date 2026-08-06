#ifndef VSTABLEMEM_H
#define VSTABLEMEM_H

#include "vstable.h"
#include <vector>

class VSTableMem : public VSTable
{
public:
    VSTableMem();
    virtual ~VSTableMem();
    
    virtual bool addCol(std::string name) override;
    virtual bool addCol(std::vector<std::string> listOfCol) override;

    virtual bool createTable() override;
    virtual bool createTable(float **fArray) override;

    virtual int getColumn(int colNumber, float *Col, int fromRow, int toRow) override;
    virtual int getColumn(unsigned int *colList, unsigned int nOfCol, unsigned long long int fromRow, unsigned long long int toRow, float **fArray) override;
    virtual int getColumnList(unsigned int *colList, unsigned int nOfCol, unsigned long long int *list, int nOfEle, float **fArray) override;

    virtual int putColumn(unsigned int *colList, unsigned int nOfCol, unsigned long long int fromRow, unsigned long long int toRow, float **fArray) override;
    virtual int putColumnList(unsigned int *colList, unsigned int nOfCol, unsigned long long int *list, int nOfEle, float **fArray) override;
};

#endif
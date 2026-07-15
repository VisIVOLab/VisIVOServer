/***************************************************************************
 *   Copyright (C) 2026                                                    *
 *                                                                         *
 *   Add rows from an ASCII file to an existing VBT                        *
 *                                                                         *
 ***************************************************************************/

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <ctime>

#include "vstable.h"
#include "vsaddrowsop.h"
#include "VisIVOFiltersConfigure.h"

//---------------------------------------------------------------------
VSAddRowsOp::VSAddRowsOp()
//---------------------------------------------------------------------
{
}

//---------------------------------------------------------------------
VSAddRowsOp::~VSAddRowsOp()
//---------------------------------------------------------------------
{
}

//---------------------------------------------------------------------
void VSAddRowsOp::printHelp()
//---------------------------------------------------------------------
{
    std::cout<<"Append rows from an ASCII file to an existing VBT"<<std::endl;

    std::cout<<"Usage:"<<std::endl;

    std::cout<<"VisIVOFilters --op addrows "
         <<"--rows rows.txt "
         <<"[--file] inputFile.bin"
         <<std::endl;

    std::cout<<"Options:"<<std::endl;

    std::cout<<"--rows  ASCII file containing rows to append"<<std::endl;
    std::cout<<"--out   Output VBT filename"<<std::endl;

    std::cout<<std::endl;

    std::cout<<"ASCII file format:"<<std::endl;
    std::cout<<"Each line must contain the same number of "
             <<"elements as the VBT columns"<<std::endl;

    std::cout<<std::endl;

    std::cout<<std::endl;
}

//---------------------------------------------------------------------
bool VSAddRowsOp::execute()
//---------------------------------------------------------------------
{
    std::string rowsFile;
    std::string outputVBT;

    // Read parameters
    
    if(m_tables.size()==0 || m_tables[0]==NULL)
    {
        std::cerr<<"VSAddRowsOp: no input table loaded"<<std::endl;
        return false;
    }

    rowsFile = getParameterAsString("rows");

    if(rowsFile.empty() || rowsFile=="unknown")
    {
        std::cerr<<"VSAddRowsOp: missing ASCII input file"<<std::endl;
        return false;
    }

    outputVBT = getParameterAsString("out");

    if(outputVBT.empty())
    {
        time_t rawtime;
        struct tm * timeinfo;
        char buffer[80];

        time(&rawtime);
        timeinfo = localtime(&rawtime);

        strftime(buffer,80,"%Y%m%d%H%M",timeinfo);
        std::string locator = m_tables[0]->getLocator();

        outputVBT =
            locator.substr(0, locator.length()-4)
            + "_addrows_"
            + buffer
            + ".bin";
    }
    std::clog << "Output Name" << outputVBT << std::endl;
    if(outputVBT.find(".bin") == std::string::npos)
        outputVBT.append(".bin");


    unsigned int nCols =
        m_tables[0]->getNumberOfColumns();

    unsigned long long oldRows =
        m_tables[0]->getNumberOfRows();

    // Read ASCII rows

    std::ifstream infile(rowsFile.c_str());

    if(!infile)
    {
        std::cerr<<"VSAddRowsOp: cannot open ASCII file "
                 <<rowsFile<<std::endl;
        return false;
    }

    std::vector< std::vector<float> > newRows;

    std::string line;

    while(std::getline(infile, line))
    {
        if(line=="")
            continue;

        std::stringstream ss(line);

        std::vector<float> row;

        float value;

        while(ss >> value)
            row.push_back(value);

        if(row.size() != nCols)
        {
            std::cerr<<"VSAddRowsOp: wrong number of columns "
                     <<"in line:"<<std::endl;

            std::cerr<<line<<std::endl;

            std::cerr << "Expected " <<nCols << " values, found " <<row.size() <<std::endl;

            return false;
        }

        newRows.push_back(row);
    }

    infile.close();

    unsigned long long addedRows =
        newRows.size();

    if(addedRows == 0)
    {
        std::cerr<<"VSAddRowsOp: no rows found in ASCII file"<<std::endl;
        return false;
    }
    unsigned long long totalRows =
        oldRows + addedRows;

    std::cout<<"Original rows : "<<oldRows<<std::endl;
    std::cout<<"Rows to append: "<<addedRows<<std::endl;
    std::cout<<"Total rows    : "<<totalRows<<std::endl;

    // Create output table
    
    remove(outputVBT.c_str());

    VSTable outTable;

    for(unsigned int i=0;i<nCols;i++)
        outTable.addCol(m_tables[0]->getColName(i));

    #ifdef VSBIGENDIAN
    std::string endianism="big";
    #else
    std::string endianism="little";
    #endif

    outTable.setLocator(outputVBT);
    outTable.setEndiannes(endianism);
    outTable.setType("float");
    outTable.setNumberOfRows(totalRows);

    outTable.writeHeader();

    // Copy original data

    int maxInt = getMaxNumberInt();

    unsigned long long maxEle;

    if(oldRows > (unsigned long long)maxInt)
        maxEle = maxInt;
    else
        maxEle = oldRows;

    if(maxEle == 0)
        maxEle = 1;

    float **fArray = new float*[1];

    try
    {
        fArray[0] = new float[maxEle];
    }
    catch(std::bad_alloc &e)
    {
        std::cerr<<"Memory allocation failed"<<std::endl;
        return false;
    }

    unsigned int colList[1];
    int nOfCol = 1;

    // Copy old table

    for(unsigned int col=0; col<nCols; col++)
    {
        colList[0] = col;

        unsigned long long startCounter = 0;
        unsigned long long remaining = oldRows;

        while(remaining > 0)
        {
            unsigned long long fromRow =
                startCounter;

            unsigned long long toRow =
                oldRows - 1;

            if((toRow - fromRow + 1) > maxEle)
                toRow = fromRow + maxEle - 1;

            m_tables[0]->getColumn(colList, nOfCol, fromRow, toRow, fArray);

            outTable.putColumn(colList, nOfCol, fromRow, toRow, fArray);

            unsigned long long copied =
                toRow - fromRow + 1;

            remaining -= copied;

            startCounter = toRow + 1;
        }
    }

    // Append new rows

    float **appendArray = new float*[1];

    try
    {
        appendArray[0] = new float[addedRows];
    }
    catch(std::bad_alloc &e)
    {
        std::cerr<<"Append allocation failed"<<std::endl;
        return false;
    }

    for(unsigned int col=0; col<nCols; col++)
    {
        for(unsigned long long row=0; row<addedRows; row++) appendArray[0][row] = newRows[row][col];


        colList[0] = col;

        outTable.putColumn(colList, nOfCol, oldRows, totalRows - 1, appendArray);
    }

    // Cleanup

    if(fArray)
    {
        if(fArray[0])
            delete [] fArray[0];

        delete [] fArray;
    }

    if(appendArray)
    {
        if(appendArray[0])
            delete [] appendArray[0];

        delete [] appendArray;
    }

    m_realOutFilename.push_back(outputVBT);

    return true;
}

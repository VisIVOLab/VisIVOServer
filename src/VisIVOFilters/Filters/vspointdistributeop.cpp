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
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#ifdef WIN32
#include <time.h>
#endif
#include "vspointdistributeop.h"
#include "vstable.h"
#include "VisIVOFiltersConfigure.h"

const unsigned int VSPointDistributeOp::MAX_NUMBER_TO_REDUCE_ROW = 100000;
const unsigned int VSPointDistributeOp::MIN_NUMBER_OF_ROW = 100;

//---------------------------------------------------------------------
VSPointDistributeOp::VSPointDistributeOp()
//---------------------------------------------------------------------
{
    m_modelBounds[0] = 0.0;
    m_modelBounds[1] = 0.0;
    m_modelBounds[2] = 0.0;
    m_modelBounds[3] = 0.0;
    m_modelBounds[4] = 0.0;
    m_modelBounds[5] = 0.0;
    
    m_sampleDimensions[0] = 50;
    m_sampleDimensions[1] = 50;
    m_sampleDimensions[2] = 50;
    
    m_nullValue = 0.0;
    
    m_splattedScalar = -1;
    
    m_useConstant = false;
    m_constValue=1.0;
    m_fArray= NULL;
    m_grid=NULL;
    m_executeDone= false;
    m_tsc=false;
    m_cic=true;
    m_ngp=false;
    m_OriginSet=false;
    m_SpacingSet=false;
    m_gridSpacing=false;
    m_avg=false;
}

//---------------------------------------------------------------------
VSPointDistributeOp::~VSPointDistributeOp()
//---------------------------------------------------------------------
{
	if(m_fArray !=NULL)
	{
        for(int i=0;i<3;i++)
            delete [] m_fArray[i];
        delete [] m_fArray;
	}
	if(m_grid != NULL)
		delete [] m_grid;
}
//---------------------------------------------------------------------
void VSPointDistributeOp::printHelp()
//---------------------------------------------------------------------
{
    std::cout<<"It produces a table which represent a volume from selected fields of the input table that are distributed using NGP , CIC (default) or TSC algorithm"<<std::endl<<std::endl;
    std::cout<<"Usage: VisIVOFilters --op pointdistribute  --resolution x_res y_res z_res --points x_col y_col z_col [--field column_names] [--nodensity] [--avg] [--out filename_out.bin] [--tsc] [--ngp] [--gridOrigin xg0 xg1 xg2] [--gridSpacing sg0 sg1 sg2]  [--box length] [--periodic] [--history] [--historyfile filename.xml] [--help] [--file] inputFile.bin"<<std::endl<<std::endl;
    
    std::cout<<"Example: VisIVOFilters --op pointdistribute --resolution 16 16 16 --points X Y Z --field Mass Temperature   --out filename_out.bin --file inputFile.bin"<<std::endl;
    
    std::cout<<"Note:  "<<std::endl;
    std::cout<<"--resolution  3D mesh size."<<std::endl;
    std::cout<<"--points Columns to be assumed for points coordinates."<<std::endl;
    std::cout<<"--field Valid columns name list to be distributed in the grid."<<std::endl;
    std::cout<<"--constant Assign a constant to all points to be distributed in the grid Ignored if field option is given. Default value is a 1.0 for all points."<<std::endl;
    std::cout<<"--nodensity Overrides the default behavior. The field distribution is not divided for the cell volume."<<std::endl;
    std::cout<<"--avg Distributes the first field on the volume grid and compute the aritmethic average value on each cell, of the first field. The output volume table will have three field. For each cell the number of total elements in the cell NumberOfElements, the sum of total field value fieldSum, and the aritmetic average value fieldAvg. Only the ngp algorithm will be applied."<<std::endl;
    std::cout<<"--out Name of the new table. Default name is given."<<std::endl;
    std::cout<<"--tsc. The TSC algorithm is adopted."<<std::endl;
    std::cout<<"--ngp. The NGP algorithm is adopted."<<std::endl;
    std::cout<<"--gridOrigin. It specifies the coordinates of the lower left corner of the grid. Default values are assumed from the box of inputFile.bin"<<std::endl;
    std::cout<<"--gridSpacing. It specifies the length of each cell dimension in arbitray unit. This parameter is ignored if the box option is given. Default vaules are assumed from the box of inputFile.bin"<<std::endl;
    std::cout<<"--box. It specifies the length of a box. Default value is assumed from the box of inputFile.bin if the gridSpacing option is not given"<<std::endl;
    std::cout<<"--periodic. It specifies the box is periodic. Particles outside the box limits are considered inside on the other side."<<std::endl;
    std::cout<<"--history (optional) create an XML file which contains the history of operations performed (default create hist.xml file)"<<std::endl;
    std::cout<<"--historyfile [filename.xml]   (optional) Change default history file name  and or directory "<<std::endl;

    std::cout<<"--file Input table filename."<<std::endl;
    
    std::cout<<"--help produce this output."<<std::endl;
    
    return;
    
}
//---------------------------------------------------------------------
bool VSPointDistributeOp::getOrigin(float *origin)
//---------------------------------------------------------------------
{
	if(!m_executeDone)
		return false;
	origin[0]=m_origin[0];
	origin[1]=m_origin[1];
	origin[2]=m_origin[2];
	return true;
}
//---------------------------------------------------------------------
bool VSPointDistributeOp::getSpacing(float *spacing)
//---------------------------------------------------------------------
{
	if(!m_executeDone)
		return false;
	spacing[0]=m_spacing[0];
	spacing[1]=m_spacing[1];
	spacing[2]=m_spacing[2];
	return true;
    
}
//---------------------------------------------------------------------
bool VSPointDistributeOp::allocateArray(int nField)
//---------------------------------------------------------------------
{
    unsigned long long int tempLL=getMaxNumberInt();
    if(((unsigned long long int)m_nOfRow*m_nOfCol)>tempLL)
        m_nOfRow=(int)tempLL/m_nOfCol;
    
    try
    {
        m_fArray=new  float*[m_nOfCol];
    }
    catch(std::bad_alloc &e)
    {
        m_fArray=NULL;
    }
    
	if(m_fArray==NULL)
		return false;
    if(m_avg) nField=3;
    try
    {
        m_grid=new  float*[nField];
    }
    catch(std::bad_alloc &e)
    {
        m_grid=NULL;
    }
    
	if(m_grid==NULL)
		return false;
    
	bool goodAllocation=false;
	while(!goodAllocation)
	{
		goodAllocation=true;
		for(unsigned int i=0;i<m_nOfCol;i++)
		{
            try
            {
                m_fArray[i]=new  float[m_nOfRow];
            }
            catch(std::bad_alloc &e)
            {
                m_fArray[i]=NULL;
            }
            
			if(m_fArray[i]==NULL)
			{
				goodAllocation=false;
				for(unsigned int j=0;j<i;j++)
					delete [] m_fArray[j];
				if(m_nOfRow==MIN_NUMBER_OF_ROW)
				{
					delete [] m_fArray;
					m_fArray=NULL;
					return false;
				}
				m_nOfRow=m_nOfRow-MAX_NUMBER_TO_REDUCE_ROW;
				if(m_nOfRow<=MAX_NUMBER_TO_REDUCE_ROW) m_nOfRow=MIN_NUMBER_OF_ROW;
				break;
			}
            //		std::clog<<i<<" " <<m_nOfRow<<std::endl;
		}
		if(!goodAllocation)
			continue;
		for(unsigned int i=0;i<nField;i++)
		{
            try
            {
                m_grid[i]=new  float[m_numNewPts];
            }
            catch(std::bad_alloc &e)
            {
                m_grid[i]=NULL;
            }
            
			if(m_grid[i]==NULL)
			{
				goodAllocation=false;
				for(unsigned int j=0;j<i;j++)
					delete [] m_grid[j];
				for(unsigned int j=0;j<m_nOfCol;j++)
					delete [] m_fArray[j];
				if(m_numNewPts==MIN_NUMBER_OF_ROW)
				{
					delete [] m_fArray;
					delete [] m_grid;
					m_fArray=NULL;
					m_grid=NULL;
					return false;
				}
				m_nOfRow=m_nOfRow-MAX_NUMBER_TO_REDUCE_ROW;
				if(m_nOfRow<=MAX_NUMBER_TO_REDUCE_ROW)
					m_nOfRow=MIN_NUMBER_OF_ROW;
				m_numNewPts=m_numNewPts-MAX_NUMBER_TO_REDUCE_ROW;
				if(m_numNewPts<=MAX_NUMBER_TO_REDUCE_ROW)
					m_numNewPts=MIN_NUMBER_OF_ROW;
				break;
			}
            
		}
        
	}
	return true;
}

//---------------------------------------------------------------------
bool VSPointDistributeOp::computeModelBounds()
//---------------------------------------------------------------------
// Compute ModelBounds from input geometry.
// bounds are the lower and upper coordinates of the particles
// in code units
// serach for min and max coordinates in the table
{
    
    unsigned int counterCols=3;
    unsigned long long int totRows=m_tables[0]->getNumberOfRows();
    
    unsigned long long int totEle=totRows;
    unsigned long long int fromRow, toRow, startCounter=0;
    unsigned int nOfValidElement=0;
    float maxValue[3],minValue[3];
    while(totEle!=0)
    {
        fromRow=startCounter;
        toRow=fromRow+m_nOfRow-1;
        if(toRow>totRows-1)
            toRow=totRows-1;
        m_tables[0]->getColumn(m_colList, counterCols, fromRow, toRow, m_fArray);
        if(startCounter==0)
        {
            for(int k=0;k<3;k++)
            {
                maxValue[k]=m_fArray[k][0];
                minValue[k]=m_fArray[k][0];
            }
        }
        for(unsigned int j=0;j<(toRow-fromRow+1);j++)
        {
            for(int k=0;k<3;k++)
            {
                if(maxValue[k]<m_fArray[k][j]) maxValue[k]=m_fArray[k][j];
                if(minValue[k]>m_fArray[k][j]) minValue[k]=m_fArray[k][j];
            }
        }
		startCounter=toRow+1;
		totEle=totEle-(toRow-fromRow+1);
		if(totEle<0) totEle=0;
    }
    //  for(int i=0; i<3; i++) std::clog<<"i="<<i<<" min="<<minValue[i]<<std::endl;//AA
    //  for(int i=0; i<3; i++) std::clog<<"i="<<i<<" max="<<maxValue[i]<<std::endl;//AA
    
    //// END search
    
    for(int i=0; i<3; i++) m_modelBounds[2*i] = minValue[i];
    for(int i=0; i<3; i++) m_modelBounds[2*i+1] = maxValue[i];
    
    // Set volume origin and data spacing
    
    for (int i=0; i<3; i++)
    {
        if(!m_OriginSet)
            m_origin[i] = m_modelBounds[2*i];
        if(!m_SpacingSet)
            m_spacing[i] = (m_modelBounds[2*i+1] - m_origin[i]) / m_sampleDimensions[i];
    }
    return true;
}

//---------------------------------------------------------------------
bool VSPointDistributeOp::setOrigin()
{
    
    if(getParameterAsString("gridOrigin").empty()||getParameterAsString("gridOrigin")=="unknown")
    {
        std::cerr<<"No valid input gridOrigin parameter."<<std::endl;
        return false;
    }
    std::stringstream ssOrigin;
    ssOrigin.str(getParameterAsString("gridOrigin"));
    int count=0;
    while (!ssOrigin.eof())
    {
        ssOrigin>>m_origin[count];
        if(count==2)
            break;
        count++;
    }
    
    return true;
}
//---------------------------------------------------------------------
bool VSPointDistributeOp::setSpacing()
{
    if(isParameterPresent("box"))
    {
        float box=getParameterAsFloat("box");
        if(box<=0.0)
        {
            std::cerr<<"No valid input box parameter."<<std::endl;
            return false;
        }
        for(int i=0;i<3;i++)
        {
            if(m_sampleDimensions[i]<=0)
            {
                std::cerr<<"No valid grid cell number."<<std::endl;
                return false;
            }
            m_spacing[i]=box/m_sampleDimensions[i];
        }
        return true;
    }
    if(isParameterPresent("gridSpacing"))
    {
        if(getParameterAsString("gridSpacing").empty()||getParameterAsString("gridSpacing")=="unknown")
        {
            std::cerr<<"No valid input gridSpacing parameter. Operation aborted"<<std::endl;
            return false;
        }
        std::stringstream ssSpacing;
        ssSpacing.str(getParameterAsString("gridSpacing"));
        int count=0;
        while (!ssSpacing.eof())
        {
            ssSpacing>>m_spacing[count];
            if(m_spacing[count] <=0.)
            {
                std::cerr<<"No valid input gridSpacing parameter."<<std::endl;
                return false;
            }
            if(count==2)
                break;
            count++;
        }
    }
    return true;
}
//---------------------------------------------------------------------

int VSPointDistributeOp::parsePointColumns(unsigned int colList[3]) {
    std::stringstream ss(getParameterAsString("points"));
    std::string paramField;
    int counterCols = 0;

    while (ss >> paramField) {  
        int colId = m_tables[0]->getColId(paramField);
        if (colId >= 0) {
            colList[counterCols] = static_cast<unsigned int>(colId);
            counterCols++;
            if (counterCols == 3) break;  // Stop after finding 3 valid columns
        }
    }

    return counterCols;  // Return the number of valid columns found
}

void VSPointDistributeOp::setAlgorithm() {

    if (isParameterPresent("tsc")) {
        m_tsc = true;
        m_cic = false;
        m_ngp = false;
    }
    if (isParameterPresent("ngp")) {
        if (m_tsc) {
            std::cerr << "Ignored --tsc parameter because --ngp was also specified." << std::endl;
        }
        m_cic = false;
        m_tsc = false;
        m_ngp = true;
    }
    if (m_avg) {
        // Force NGP when averaging
        m_tsc = false;
        m_cic = false;
        m_ngp = true;
    }
}
bool VSPointDistributeOp::setGridResolution() {
    std::stringstream ssResolution(getParameterAsString("resolution"));
    int parsedDimensions = 0;

    while (ssResolution >> m_sampleDimensions[parsedDimensions]) {
        if (m_sampleDimensions[parsedDimensions] <= 0) {
            std::cerr << "VSPointDistributeOp: Invalid resolution value given: " 
                      << m_sampleDimensions[parsedDimensions] << std::endl;
            return false;
        }
        parsedDimensions++;
        if (parsedDimensions == 3) break;  // Stop at 3 dimensions
    }

    if (parsedDimensions < 3) {
        std::cerr << "VSPointDistributeOp: Invalid resolution—must provide exactly 3 values." << std::endl;
        return false;
    }
    return true;
}

bool VSPointDistributeOp::parseFieldList(std::vector<int>& fieldList) {
    std::stringstream fieldNameSStream(getParameterAsString("field"));

    if (fieldNameSStream.str().empty() || fieldNameSStream.str() == "unknown") {
        m_useConstant = true;
        fieldList.push_back(-1);
        if (isParameterPresent("constant")) {
            m_constValue = getParameterAsFloat("constant");
        }
    } else {
        std::string paramField;
        while (fieldNameSStream >> paramField) {
            int colId = m_tables[0]->getColId(paramField);
            if (colId >= 0) {
                fieldList.push_back(colId);
            }
            if (m_avg) break;  // Only one field is needed if averaging
        }
    }

    if (fieldList.empty()) {
        std::cerr << "PointDistribute: Invalid field specified." << std::endl;
        return false;
    }

    return true;
}

bool VSPointDistributeOp::allocateColumnList(unsigned int*& colList, const std::vector<int>& fieldList) {
    m_nOfCol = m_useConstant ? 3 : 3 + fieldList.size();  // Adjust column count

    try {
        colList = new unsigned int[m_nOfCol];  // Allocate memory
    } catch (std::bad_alloc&) {
        std::cerr << "Failed array allocation. VSPointDistributeOp operation terminated." << std::endl;
        return false;
    }
    return true;
}

bool VSPointDistributeOp::initializeGrid(const std::vector<int>& fieldList, unsigned int* colList) {
    unsigned long long int totRows = m_tables[0]->getNumberOfRows();
    int maxInt = getMaxNumberInt();

    m_nOfRow = (totRows > maxInt) ? maxInt : totRows;

    m_gridPts = m_sampleDimensions[0] * m_sampleDimensions[1] * m_sampleDimensions[2];
    m_numNewPts = (m_gridPts > maxInt) ? maxInt : m_gridPts;

    bool allocationArray = allocateArray((int) fieldList.size());

    if (m_fArray == nullptr || m_grid == nullptr || !allocationArray) {
        std::cerr << "Failed Array allocation. vspointdistribute Operation terminated" << std::endl;
        delete[] colList;
        return false;
    }
    return true;
}

bool VSPointDistributeOp::processGridSpacing() {
    m_gridSpacing = true;
    m_SpacingSet = true;

    std::stringstream ssgridSpacing(getParameterAsString("gridSpacing"));
    int counterCols = 0;

    while (!ssgridSpacing.eof()) {
        ssgridSpacing >> m_spacing[counterCols];

        if (m_spacing[counterCols] <= 0.) {
            std::cerr << "Invalid gridSpacing values. Pointdistribute Operation terminated" << std::endl;
            return false;
        }

        counterCols++;
        if (counterCols == 3) break;
    }

    return true;
}

std::string VSPointDistributeOp::generateOutputFileName() {
    std::stringstream fileNameOutputSStream;
    fileNameOutputSStream << getParameterAsString("out");

    if (fileNameOutputSStream.str().empty() || fileNameOutputSStream.str() == "unknown") {
        std::string filenameInputTable = m_tables[0]->getLocator();
        int len = filenameInputTable.length();

        // Generate timestamp
        time_t rawtime;
        struct tm* timeinfo;
        char buffer[80];
        time(&rawtime);
        timeinfo = localtime(&rawtime);
        strftime(buffer, 80, "%Y%m%d%H%M", timeinfo);

        // Create a default filename
        fileNameOutputSStream.str("");
        fileNameOutputSStream << filenameInputTable.substr(0, len - 4) << "_pointdistribute_" << buffer << ".bin";
    }

    std::string fileNameOutput = fileNameOutputSStream.str();
    if (fileNameOutput.find(".bin") == std::string::npos) {
        fileNameOutput.append(".bin");
    }

    return fileNameOutput;
}

void VSPointDistributeOp::configureTableGrid(VSTable& tableGrid, const std::vector<int>& fieldList) {
    #ifdef VSBIGENDIAN
    std::string endianism="big";
#else
    std::string endianism="little";
#endif    
    tableGrid.setEndiannes(endianism);
    tableGrid.setType("float");

    if (m_avg) {
        std::stringstream fileNameColSStream;
        if (m_useConstant)
            fileNameColSStream << "Constant";
        else
            fileNameColSStream << m_tables[0]->getColName(fieldList[0]);

        tableGrid.addCol("NumberOfElements");
        tableGrid.addCol(fileNameColSStream.str() + "Sum");
        tableGrid.addCol(fileNameColSStream.str() + "Avg");
    } 
    else if (m_useConstant) {
        tableGrid.addCol("Constant");
    } 
    else {
        for (int colId : fieldList)
            tableGrid.addCol(m_tables[0]->getColName(colId));
    }

    // Set table dimensions
    tableGrid.setNumberOfRows(m_gridPts);
    tableGrid.setIsVolume(true);
    tableGrid.setCellNumber(m_sampleDimensions[0], m_sampleDimensions[1], m_sampleDimensions[2]);

    // Configure spacing
    float spacing[3];
    if (m_gridSpacing) {
        spacing[0] = m_spacing[0];
        spacing[1] = m_spacing[1];
        spacing[2] = m_spacing[2];
    } else {
        spacing[0] = 1.0;
        spacing[1] = m_spacing[1] / m_spacing[0];
        spacing[2] = m_spacing[2] / m_spacing[0];
    }
    tableGrid.setCellSize(spacing[0], spacing[1], spacing[2]);

    // Write header to the table
    tableGrid.writeHeader();
}

bool VSPointDistributeOp::initializeEmptyGrid(VSTable& tableGrid, const std::vector<int>& fieldList, unsigned int*& gridList) {
    // Allocate grid list
    int numOfField = fieldList.size();
    if (m_avg) numOfField = 3;  // Override field count in case of averaging

    try {
        gridList = new unsigned int[numOfField];
    } catch (std::bad_alloc&) {
        std::cerr << "Memory allocation failed for gridList." << std::endl;
        return false;
    }

    // Initialize gridList and zero-fill m_grid
    for (int k = 0; k < numOfField; k++) {
        gridList[k] = k;
        for (unsigned int i = 0; i < m_numNewPts; i++)
            m_grid[k][i] = 0.0;
    }

    // Fill the table with zeros
    gridHandle.totEle = m_gridPts;

    while (gridHandle.totEle > 0) {
        gridHandle.fromRow = gridHandle.startCounter;
        gridHandle.toRow = std::min(gridHandle.fromRow + m_numNewPts - 1, m_gridPts - 1);

        tableGrid.putColumn(gridList, numOfField, gridHandle.fromRow, gridHandle.toRow, m_grid);

        gridHandle.totEle -= (gridHandle.toRow - gridHandle.fromRow + 1);
        gridHandle.startCounter = gridHandle.toRow + 1;
    }

    return true;
}

void VSPointDistributeOp::applyPeriodicBoundary_CIC(int& i1, int& i2, int& i3, int& i11, int& i21, int& i31) {
    if (!m_periodic) return;

    auto applyWrap = [](int& val, int maxVal) {
        if (val < 0) val += maxVal;
        if (val >= maxVal) val -= maxVal;
    };

    applyWrap(i1, m_sampleDimensions[0]);
    applyWrap(i2, m_sampleDimensions[1]);
    applyWrap(i3, m_sampleDimensions[2]);
    applyWrap(i11, m_sampleDimensions[0]);
    applyWrap(i21, m_sampleDimensions[1]);
    applyWrap(i31, m_sampleDimensions[2]);
}

bool VSPointDistributeOp::processCIC(VSTable& tableGrid, unsigned int* gridList, int nOfField, 
                                     std::vector<int>& fieldList, unsigned long long* gridIndex) {
    float wc = 0.0;
    int nCell = m_sampleDimensions[0] * m_sampleDimensions[1] * m_sampleDimensions[2];
    int jkFactor = m_sampleDimensions[0]*m_sampleDimensions[1];
    int jFactor = m_sampleDimensions[0];
    float cellVolume=m_spacing[0]*m_spacing[1]*m_spacing[2];
    float norm;
    for (int ptId = 0; ptId < gridHandle.toRow - gridHandle.fromRow + 1; ptId++) {
        wc = 0.0;
        float px[3] = { m_fArray[0][ptId], m_fArray[1][ptId], m_fArray[2][ptId] };

        float pos1 = (px[0] - m_origin[0]) / m_spacing[0];
        float pos2 = (px[1] - m_origin[1]) / m_spacing[1];
        float pos3 = (px[2] - m_origin[2]) / m_spacing[2];

        int i1 = floor(pos1);
        int i2 = floor(pos2);
        int i3 = floor(pos3);
        int i11 = i1 + 1;
        int i21 = i2 + 1;
        int i31 = i3 + 1;

        applyPeriodicBoundary_CIC(i1, i2, i3, i11, i21, i31);

        // Compute CIC Weights
        float weights[8];
        computeCICWeights(pos1, pos2, pos3, weights);

        // Linearize coordinates
        unsigned long long int ind[8] = {
            i3 * jkFactor + i2 * jFactor + i1,
            i3 * jkFactor + i2 * jFactor + i11,
            i3 * jkFactor + i21 * jFactor + i1,
            i3 * jkFactor + i21 * jFactor + i11,
            i31 * jkFactor + i2 * jFactor + i1,
            i31 * jkFactor + i2 * jFactor + i11,
            i31 * jkFactor + i21 * jFactor + i1,
            i31 * jkFactor + i21 * jFactor + i11
        };

        // Process density assignment
        for (int n = 0; n < 8; n++) {
            if (ind[n] < 0 || ind[n] >= nCell) continue;

            if (ind[n] < gridIndex[0] || ind[n] > gridIndex[1]) { // Not in cache
                tableGrid.putColumn(gridList, nOfField, gridIndex[0], gridIndex[1], m_grid);
                gridIndex[0] = ind[n];
                gridIndex[1] = gridIndex[0] + m_numNewPts - 1;
                if (gridIndex[1] >= m_gridPts) gridIndex[1] = m_gridPts - 1;
                tableGrid.getColumn(gridList, tableGrid.getNumberOfColumns(), gridIndex[0], gridIndex[1], m_grid);
            }

            for(int j=0;j<fieldList.size();j++)
                    {
                        if(m_useConstant)
                            norm=m_constValue;
                        else
                            norm=m_fArray[3+j][ptId];
                        m_grid[j][ind[n]-gridIndex[0]]+=weights[n]*norm/cellVolume;
                        wc+=weights[n]*norm;
                        //		outpippo<<"ptId="<<ptId<<" GRID j="<<j<<" i="<<ind[n]-gridIndex[0] <<" curr val="<<d[n]*norm<<" acc="<<m_grid[j][ind[n]-gridIndex[0]]<<std::endl; //AA
                    }
        }

        // Error check
        if (wc > 1.1 * norm * fieldList.size()) {
            std::cerr << "Error 2 on CIC schema. Operation Aborted" << std::endl;
            return false;
        }
    }

    return true;
}

void VSPointDistributeOp::applyPeriodicBoundary_TSC(float px[3]) {
    if (!m_periodic) return;

    for (int i = 0; i < 3; i++) {
        if (px[i] < 0) px[i] += m_sampleDimensions[i] * m_spacing[i];
        if (px[i] >= (m_sampleDimensions[i] * m_spacing[i])) px[i] -= m_sampleDimensions[i] * m_spacing[i];
    }
}

int VSPointDistributeOp::getLinearizedIndex(int x, int y, int z) {
    if (m_periodic) {
        x = (x + m_sampleDimensions[0]) % m_sampleDimensions[0];
        y = (y + m_sampleDimensions[1]) % m_sampleDimensions[1];
        z = (z + m_sampleDimensions[2]) % m_sampleDimensions[2];
    }
    if (x < 0 || x >= m_sampleDimensions[0] || y < 0 || y >= m_sampleDimensions[1] || z < 0 || z >= m_sampleDimensions[2])
        return -1;

    return x + m_sampleDimensions[0] * y + m_sampleDimensions[0] * m_sampleDimensions[1] * z;
}

void VSPointDistributeOp::computeCICWeights(float pos1, float pos2, float pos3, float weights[8]) {
    float dd1 = pos1 - floor(pos1);
    float dd2 = pos2 - floor(pos2);
    float dd3 = pos3 - floor(pos3);

    float de1 = 1.0 - dd1;
    float de2 = 1.0 - dd2;
    float de3 = 1.0 - dd3;

    weights[0] = de1 * de2 * de3;
    weights[1] = dd1 * de2 * de3;
    weights[2] = de1 * dd2 * de3;
    weights[3] = dd1 * dd2 * de3;
    weights[4] = de1 * de2 * dd3;
    weights[5] = dd1 * de2 * dd3;
    weights[6] = de1 * dd2 * dd3;
    weights[7] = dd1 * dd2 * dd3;
}

float VSPointDistributeOp::computeTSCWeight(float posX, float posY, float posZ, int gridX, int gridY, int gridZ) {
    float dist_x = fabs((posX / m_spacing[0]) - gridX);
    float dist_y = fabs((posY / m_spacing[1]) - gridY);
    float dist_z = fabs((posZ / m_spacing[2]) - gridZ);

    float w_x = (dist_x <= .5) ? (0.75 - dist_x * dist_x) : (dist_x <= 1.5 ? 0.5 * (1.5 - dist_x) * (1.5 - dist_x) : 0);
    float w_y = (dist_y <= .5) ? (0.75 - dist_y * dist_y) : (dist_y <= 1.5 ? 0.5 * (1.5 - dist_y) * (1.5 - dist_y) : 0);
    float w_z = (dist_z <= .5) ? (0.75 - dist_z * dist_z) : (dist_z <= 1.5 ? 0.5 * (1.5 - dist_z) * (1.5 - dist_z) : 0);

    return w_x * w_y * w_z;
}

void VSPointDistributeOp::processTSC(VSTable& tableGrid, unsigned int* gridList, int nOfField,
                                     std::vector<int>& fieldList, unsigned long long* gridIndex) {
    float wc = 0;

    float cellVolume=m_spacing[0]*m_spacing[1]*m_spacing[2];
    int nCell = m_sampleDimensions[0] * m_sampleDimensions[1] * m_sampleDimensions[2];

    for (int ptId = 0; ptId < gridHandle.toRow - gridHandle.fromRow + 1; ptId++) {
        float px[3] = {
            m_fArray[0][ptId] - m_origin[0],
            m_fArray[1][ptId] - m_origin[1],
            m_fArray[2][ptId] - m_origin[2]
        };

        applyPeriodicBoundary_TSC(px);

        int ind_x = (int)floor(px[0] / m_spacing[0]);
        int ind_y = (int)floor(px[1] / m_spacing[1]);
        int ind_z = (int)floor(px[2] / m_spacing[2]);

        wc = 0.;
        int xList[4], yList[4], zList[4];

        for (int j = 0; j < 4; j++) {
            xList[j] = ind_x + j - 1;
            yList[j] = ind_y + j - 1;
            zList[j] = ind_z + j - 1;
        }

        for (int ix = 0; ix < 4; ix++) {
            for (int iy = 0; iy < 4; iy++) {
                for (int iz = 0; iz < 4; iz++) {
                    float w = computeTSCWeight(px[0], px[1], px[2], xList[ix], yList[iy], zList[iz]);
                    if (w == 0) continue;

                    int ind_test = getLinearizedIndex(xList[ix], yList[iy], zList[iz]);
                    if (ind_test < 0 || ind_test >= nCell) continue;

                    if (ind_test < gridIndex[0] || ind_test > gridIndex[1]) {
                        tableGrid.putColumn(gridList, nOfField, gridIndex[0],gridIndex[1],m_grid);
                        gridIndex[0]=ind_test;
                        gridIndex[1]=gridIndex[0]+m_numNewPts-1;
                        if(gridIndex[1]>=m_gridPts)gridIndex[1]=m_gridPts-1;
                        tableGrid.getColumn(gridList,nOfField,gridIndex[0],gridIndex[1],m_grid);
                    }

                    for (int j = 0; j < fieldList.size(); j++) {
                        float norm = m_useConstant ? m_constValue : m_fArray[3 + j][ptId];
                        m_grid[j][ind_test - gridIndex[0]] += w * norm / cellVolume;
                        wc += w * norm;
                    }
                }
            }
        }

        if (wc > 1.1) {
            std::cerr << "Error in TSC schema. Operation Aborted." << std::endl;
            return;
        }
    }
}

void VSPointDistributeOp::computeNGPIndex(float px[3], int gridPos[3]) {
    float pos[3] = {
        (px[0] - m_origin[0]) / m_spacing[0],
        (px[1] - m_origin[1]) / m_spacing[1],
        (px[2] - m_origin[2]) / m_spacing[2]
    };

    for (int i = 0; i < 3; i++) {
        gridPos[i] = floor(pos[i]);
        if (fabs(pos[i] - gridPos[i]) > 0.5) gridPos[i]++;

        if (m_periodic) {
            if (gridPos[i] < 0) gridPos[i] += m_sampleDimensions[i];
            if (gridPos[i] >= m_sampleDimensions[i]) gridPos[i] -= m_sampleDimensions[i];
        }
    }
}

void VSPointDistributeOp::processNGP(VSTable& tableGrid, unsigned int* gridList, int nOfField,
                                     std::vector<int>& fieldList, unsigned long long* gridIndex) {
    int nCell = m_sampleDimensions[0] * m_sampleDimensions[1] * m_sampleDimensions[2];
    float cellVolume=m_spacing[0]*m_spacing[1]*m_spacing[2];
    int jkFactor = m_sampleDimensions[0]*m_sampleDimensions[1];
    int jFactor = m_sampleDimensions[0];
    for (int ptId = 0; ptId < gridHandle.toRow - gridHandle.fromRow + 1; ptId++) {
        float px[3] = {
            m_fArray[0][ptId],
            m_fArray[1][ptId],
            m_fArray[2][ptId]
        };

        int gridPos[3];
        computeNGPIndex(px, gridPos);

        if (gridPos[0] < 0 || gridPos[0] >= m_sampleDimensions[0] ||
            gridPos[1] < 0 || gridPos[1] >= m_sampleDimensions[1] ||
            gridPos[2] < 0 || gridPos[2] >= m_sampleDimensions[2])
            continue;

        int ind = gridPos[2] * jkFactor + gridPos[1] * jFactor + gridPos[0];

        if (ind >= nCell) continue;

        if (ind < gridIndex[0] || ind > gridIndex[1]) { // Grid index out of cache
            tableGrid.putColumn(gridList, nOfField, gridIndex[0],gridIndex[1],m_grid); //QUI controlla:estremi INCLUSI
            gridIndex[0]=ind;
            gridIndex[1]=gridIndex[0]+m_numNewPts-1;
            if(gridIndex[1]>=m_gridPts)gridIndex[1]=m_gridPts-1;
            tableGrid.getColumn(gridList, tableGrid.getNumberOfColumns(),gridIndex[0],gridIndex[1],m_grid);
        }

        for (int j = 0; j < fieldList.size(); j++) {
            float norm = m_useConstant ? m_constValue : m_fArray[3 + j][ptId];

            if (m_avg) {
                m_grid[0][ind - gridIndex[0]] += 1.0;
                m_grid[1][ind - gridIndex[0]] += norm;
                m_grid[2][ind - gridIndex[0]] = m_grid[1][ind - gridIndex[0]] / m_grid[0][ind - gridIndex[0]];
            } else {
                m_grid[j][ind - gridIndex[0]] += norm / cellVolume;
            }
        }
    }
}

//---------------------------------------------------------------------
bool VSPointDistributeOp::execute()
//---------------------------------------------------------------------
{
    m_avg=isParameterPresent("avg");
    m_periodic=isParameterPresent("periodic");
    VSTable tableGrid;
    std::vector<int> fieldList;
    // check for points  coordinate columns
    unsigned int colLs[3];  // Fixed-size array
    unsigned long long int totRows = m_tables[0]->getNumberOfRows();
    int counterCols = parsePointColumns(colLs);

    if (counterCols != 3) {
        std::cerr << "VSPointDistributeOp: Invalid columns in --points argument" << std::endl;
        return false;
    }

    std::memcpy(m_colList, colLs, 3 * sizeof(unsigned int));  // Assign parsed columns

    setAlgorithm();
    
    //check grid resolution
    
    if (!setGridResolution()) {
        return false;
    }
    
    if (!parseFieldList(fieldList)) {
        return false;
    }

    // Prepare colList: list of columns to be read
    unsigned int nOfField = fieldList.size();
    unsigned int* colList = nullptr;
    if (!allocateColumnList(colList, fieldList)) {
        return false;
    }
    
    // allocate m_arrays
    if (!initializeGrid(fieldList, colList)) return false;
    
    //Parameters Setting: gridOrigin
    if(isParameterPresent("gridOrigin"))
        m_OriginSet=setOrigin();
    //Parameters Setting: box
    if(isParameterPresent("box"))
        m_SpacingSet=setSpacing();
    //Parameters Setting: gridSpacing
    if(isParameterPresent("gridSpacing") && !m_SpacingSet)
    {
        if (!processGridSpacing()) {
            delete[] colList;
            return false;
        }
    }
    if(!m_SpacingSet || !m_OriginSet)
        if(!computeModelBounds())
        {
            delete [] colList;
            return false;
        }

    // initialize the grid
    //open file output
    std::string fileNameOutput = generateOutputFileName();
    m_realOutFilename.push_back(fileNameOutput);

    // Clean existing tab
    remove(fileNameOutput.c_str());
    tableGrid.setLocator(fileNameOutput);
    
    configureTableGrid(tableGrid, fieldList);
    
    //  Create Empty Binary File
    unsigned int* gridList = nullptr;
    if (!initializeEmptyGrid(tableGrid, fieldList, gridList)) {
        delete[] colList;
        return false;
    }
    
    //////////////////////
    /////
    // Start Point Distribution
    ////
    
    // set cashed grid on memory
    unsigned long long int gridIndex[2];  //limits of cashed grid
    gridIndex[0]=0;
    gridIndex[1]=m_numNewPts-1;
    gridHandle.totEle=totRows;
    gridHandle.startCounter=0;
    colList[0]=m_colList[0];
    colList[1]=m_colList[1];
    colList[2]=m_colList[2];
    if(!m_useConstant)
        for(int i=0;i<nOfField;i++)
            colList[3+i]=fieldList[i];
    
    int jkFactor = m_sampleDimensions[0]*m_sampleDimensions[1];
    int jFactor = m_sampleDimensions[0];
    float norm;
    
    
    float cellVolume=m_spacing[0]*m_spacing[1]*m_spacing[2];
    if(isParameterPresent("nodensity") || m_avg)
        cellVolume=1.0;
    
    while(gridHandle.totEle!=0)
    {
        // Table downLoad
        gridHandle.fromRow=gridHandle.startCounter;
        gridHandle.toRow=gridHandle.fromRow+m_nOfRow-1;
        if(gridHandle.toRow>totRows-1)
            gridHandle.toRow=totRows-1;
        m_tables[0]->getColumn(colList,m_nOfCol, gridHandle.fromRow,gridHandle. toRow, m_fArray);
        
        if (m_ngp) {
            processNGP(tableGrid, gridList, nOfField, fieldList, gridIndex);
        }
        
        
        if (m_cic) {
            if (!processCIC(tableGrid, gridList, nOfField, fieldList, gridIndex)) {
                delete[] colList;
                delete[] gridList;
                return false;
            }
        }
        if (m_tsc) {
            processTSC(tableGrid, gridList, nOfField, fieldList, gridIndex);
        }
        gridHandle.startCounter=gridHandle.toRow+1;
        gridHandle.totEle=gridHandle.totEle-(gridHandle.toRow-gridHandle.fromRow+1);
        if(gridHandle.totEle<0) gridHandle.totEle=0;
    } 
    if(m_avg) nOfField=3;
    tableGrid.putColumn(gridList,nOfField,gridIndex[0],gridIndex[1],m_grid);  
    
    m_executeDone=true;
    if(colList!=NULL) delete [] colList;
    if(gridList!=NULL)delete [] gridList;
    return true;
}
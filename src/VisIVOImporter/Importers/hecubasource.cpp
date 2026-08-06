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
#include <cstdlib>
#include <cstring>
#include <StorageDict.h>

#include <StorageNumpy.h>
#include <StorageObject.h>
#include <KeyClass.h>
#include <ValueClass.h>
#include "hecubasource.h"
#include "visivoutils.h"
#include <iostream>
#include <fstream>
#include <sstream>

#define SIZE 3
//---------------------------------------------------------------------
int HecubaSource::readData()
{
    Dict dictTest;

    dictTest.getByAlias(m_aliasParticle);

    const int expectedItems = 1;
    int nit = 0;

    for (auto it = dictTest.begin(); it != dictTest.end(); ++it) {
        Key pk = it->first;
        Value v = it->second;

        StorageNumpy &testNumpy = Value::get<0>(v);

        if (testNumpy.data == nullptr) {
            std::cerr << "[ERROR] received StorageNumpy has null data" << std::endl;
            return -1;
        }

        if (testNumpy.metas.size() < 2) {
            std::cerr << "[ERROR] invalid StorageNumpy metadata size: "
                      << testNumpy.metas.size() << std::endl;
            return -1;
        }

        writeGasParticles(testNumpy);

        // writeDarkParticles(testNumpy);
        // writeStarParticles(testNumpy);

        ++nit;

        if (nit >= expectedItems) {
            break;
        }
    }

    if (nit == 0) {
        std::cerr << "[ERROR] no streamed item received" << std::endl;
        return -1;
    }

    return 0;
}


void HecubaSource::writeGasParticles(const StorageNumpy &s)
{
    std::vector<std::string> types;
    types.insert(types.end(), {
        "MASS", "POS_X", "POS_Y", "POS_Z",
        "VEL_X", "VEL_Y", "VEL_Z",
        "RHO", "TEMP", "EPS", "METALS", "PHI"
    });

    double* p = (double*)s.data;
    const int meta1 = s.metas[0];        // total particles (gas+dark+star combined)
    const int meta2 = s.metas[1];
    const int numCols = (int)types.size();
    const int typeCol = numCols;


    std::vector<std::vector<float>> buffers(numCols);
    for (auto &b : buffers) b.reserve(meta1);

    int gasCount = 0;
    for (int i = 0; i < meta1; i++) {
        const double* row = p + (size_t)i * meta2;
        if (row[typeCol] != 0.0) continue;

        for (int t = 0; t < numCols; t++) {
            buffers[t].push_back(static_cast<float>(row[t]));
        }
        gasCount++;
    }

    std::ofstream outfile;
    std::string pathFileOut;
    VSTable *table = nullptr;

    if (useMemory) {
        table = new VSTableMem();
        table->setType("float");
        table->setNumberOfRows(gasCount);
        for (int i = 0; i < numCols; i++) table->addCol(types[i]);

        for (int t = 0; t < numCols; t++) {
            unsigned int colList[1] = {(unsigned int)t};
            float* dataPtrs[1] = { buffers[t].data() };
            table->putColumn(colList, 1, 0, gasCount - 1, dataPtrs);
        }
    } else {
        std::string fileName = m_pointsBinaryName;
        int idx = fileName.rfind('.');
        std::string pathFileIn = fileName.erase(idx, idx + 4);
        pathFileOut = pathFileIn;
        outfile.open((pathFileOut + "GAS" + ".bin").c_str(), std::ofstream::binary);
        for (int t = 0; t < numCols; t++) {
            outfile.write((char*)buffers[t].data(), sizeof(float) * gasCount);
        }
        outfile.close();
        std::string pathHeader = pathFileOut + "GAS" + ".bin";
        makeHeader(gasCount, pathHeader, types, m_cellSize, m_cellComp, m_volumeOrTable);
    }

    if (useMemory) {
        memTables.push_back(table);
    }
}


void HecubaSource::writeDarkParticles(const StorageNumpy &s){
	std::string fileName = m_pointsBinaryName;
  	int idx = fileName.rfind('.');
	std::string pathFileIn = fileName.erase(idx, idx+4);
	std::string pathFileOut = pathFileIn;
	std::ofstream outfile((pathFileOut + "DARK" + ".bin").c_str(),std::ofstream::binary );
	std::vector<std::string> types;
  	types.push_back("MASS");
  	types.push_back("POS_X");
  	types.push_back("POS_Y");
  	types.push_back("POS_Z");
  	types.push_back("VEL_X");
  	types.push_back("VEL_Y");
  	types.push_back("VEL_Z");
  	types.push_back("EPS");
  	types.push_back("PHI");
	float *buffer=NULL;
    double* p = (double*)s.data;
    const int meta1 = s.metas[0];
    int meta2 = s.metas[1];
	buffer = new float[meta1];
	
	for(int t = 0; t < meta2; t++) {
		for(int i = 0; i < meta1; i++) {
			buffer[i] = static_cast<float>(p[i*meta2 + t]);
		}
		outfile.write((char *)(buffer), sizeof(float)*meta1);
	}
	outfile.close();

	std::string pathHeader = pathFileOut+"DARK"+ ".bin";
	makeHeader(meta1, pathHeader, types, m_cellSize,m_cellComp,m_volumeOrTable);
	return;	
}

void HecubaSource::writeStarParticles(const StorageNumpy &s){

	std::string fileName = m_pointsBinaryName;
  	int idx = fileName.rfind('.');
	std::string pathFileIn = fileName.erase(idx, idx+4);
	std::string pathFileOut = pathFileIn;
	std::ofstream outfile((pathFileOut + "STAR" + ".bin").c_str(),std::ofstream::binary );
	std::vector<std::string> types; //!species block nameset 
  	types.push_back("MASS");
  	types.push_back("POS_X");
  	types.push_back("POS_Y");
  	types.push_back("POS_Z");
  	types.push_back("VEL_X");
  	types.push_back("VEL_Y");
  	types.push_back("VEL_Z");
  	types.push_back("METALS");
  	types.push_back("TFORM");
  	types.push_back("EPS");
  	types.push_back("PHI");
	float *buffer=NULL;
    double* p = (double*)s.data;
    const int meta1 = s.metas[0];
    int meta2 = s.metas[1];
	buffer = new float[meta1];
	for(int t = 0; t < meta2; t++) {
		for(int i = 0; i < meta1; i++) {
			buffer[i] = static_cast<float>(p[i*meta2 + t]);
		}
		outfile.write((char *)(buffer), sizeof(float)*meta1);
	}
	outfile.close();

	std::string pathHeader = pathFileOut+"STAR"+ ".bin";
	makeHeader(meta1, pathHeader, types, m_cellSize,m_cellComp,m_volumeOrTable);
	return;
}

//---------------------------------------------------------------------
  int HecubaSource::readHeader()
//---------------------------------------------------------------------
{
    /*if(m_aliasHeader.empty()){
		std::cerr << "--aliasheader option not defined" << std::endl;
		return -1;
	}
  	headerObject header;
	std::clog << m_aliasHeader.c_str() << std::endl;
    //header.getByAlias(m_aliasHeader.c_str());
	std::clog << "After getByAlias Header" << std::endl;
    /*nsph = header.nsph;
    ndark = header.ndark;
    nstar = header.nstar;*/
	return 0;  
}
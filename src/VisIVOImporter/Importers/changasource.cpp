/***************************************************************************
 *   Copyright (C) 2008 by Gabriella Caniglia,Roberto Munzone *
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
//#include "VisIVOImporterConfigure.h"
#include "changasource.h"

#include "visivoutils.h"
#include "mpi.h"
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <unordered_map>
#include <filesystem>
#include <fcntl.h>
#include <omp.h>
#include <stdexcept>
#include <unistd.h>

namespace {

static inline int32_t readIntBE(std::istream &in)
{
    unsigned char b[4];
    in.read(reinterpret_cast<char *>(b), 4);
    if (!in) throw std::runtime_error("readIntBE: failed to read 4 bytes");
    uint32_t v = (uint32_t(b[0]) << 24) |
                 (uint32_t(b[1]) << 16) |
                 (uint32_t(b[2]) << 8)  |
                  uint32_t(b[3]);
    return int32_t(v);
}

static inline float readFloatBE(std::istream &in)
{
    unsigned char b[4];
    in.read(reinterpret_cast<char *>(b), 4);
    if (!in) throw std::runtime_error("readFloatBE: failed to read 4 bytes");
    uint32_t v = (uint32_t(b[0]) << 24) |
                 (uint32_t(b[1]) << 16) |
                 (uint32_t(b[2]) << 8)  |
                  uint32_t(b[3]);
    float f;
    std::memcpy(&f, &v, 4);
    return f;
}

static inline double readDoubleBE(std::istream &in)
{
    unsigned char b[8];
    in.read(reinterpret_cast<char *>(b), 8);
    if (!in) throw std::runtime_error("readDoubleBE: failed to read 8 bytes");
    uint64_t v = (uint64_t(b[0]) << 56) |
                 (uint64_t(b[1]) << 48) |
                 (uint64_t(b[2]) << 40) |
                 (uint64_t(b[3]) << 32) |
                 (uint64_t(b[4]) << 24) |
                 (uint64_t(b[5]) << 16) |
                 (uint64_t(b[6]) << 8)  |
                  uint64_t(b[7]);
    double d;
    std::memcpy(&d, &v, 8);
    return d;
}

}

//---------------------------------------------------------------------
int ChangaSource::readHeader()
{
    std::string fileName = m_pointsFileName;
    std::ifstream inFile(fileName, std::ios::binary);
    if (!inFile) {
        std::cerr << "Error while opening file in readHeader: " << fileName << std::endl;
        return -1;
    }
    
    double  time   = readDoubleBE(inFile);
    int32_t nbodies = readIntBE(inFile);
    int32_t ndim    = readIntBE(inFile);
    int32_t nsph_   = readIntBE(inFile);
    int32_t ndark_  = readIntBE(inFile);
    int32_t nstar_  = readIntBE(inFile);
    int32_t pad     = readIntBE(inFile);

    this->nsph  = nsph_;
    this->ndark = ndark_;
    this->nstar = nstar_;

    return 0;
}

int ChangaSource::readData()
{
    int idx = m_pointsBinaryName.rfind('.');
    std::string pathFileIn  = m_pointsBinaryName;
    if (idx != std::string::npos) {
        pathFileIn.erase(idx);  // remove extension
    }
    std::string pathFileOut = pathFileIn;
    std::string pathHeader;

    std::ifstream inFile(m_pointsFileName, std::ios::binary);
    if (!inFile) {
        std::cerr << "Error while opening file in readData: " << m_pointsFileName << std::endl;
        return 1;
    }

    // Skip header
    const std::streamoff headerSize = 8 + 6 * 4; // double + 6 ints = 32 bytes
    inFile.seekg(headerSize, std::ios::beg);
    if (!inFile) {
        std::cerr << "Failed to seek to particle data";
        return 1;
    }

    // -------------------------
    // GAS PARTICLES (12 floats)
    // -------------------------
    float *gasParticles = nullptr;
    if (nsph > 0) {
        gasParticles = static_cast<float *>(std::malloc(sizeof(float) * 12 * nsph));
        if (!gasParticles) {
            std::clog << "Malloc Error for gasParticles" << std::endl;
            return 1;
        }
        for (int i = 0; i < nsph; ++i) {
            for (int k = 0; k < 12; ++k) {
                gasParticles[i * 12 + k] = readFloatBE(inFile);
            }
        }
    }

    std::vector<std::string> gasBlocks;
    gasBlocks.push_back("MASS");   //0
    gasBlocks.push_back("POS_X");  //1
    gasBlocks.push_back("POS_Y");  //2
    gasBlocks.push_back("POS_Z");  //3
    gasBlocks.push_back("VEL_X");  //4
    gasBlocks.push_back("VEL_Y");  //5
    gasBlocks.push_back("VEL_Z");  //6
    gasBlocks.push_back("RHO");    //7
    gasBlocks.push_back("TEMP");   //8
    gasBlocks.push_back("EPS");    //9
    gasBlocks.push_back("METALS"); //10
    gasBlocks.push_back("PHI");    //11

    std::ofstream outfile;
    if (!useMemory) {
        outfile.open((pathFileOut + "GAS" + ".bin").c_str(), std::ofstream::binary);
        if (!outfile) {
            std::cerr << "Failed to open GAS.bin for writing" << std::endl;
            return 1;
        }
    } else {
        VSTable *table = new VSTableMem();
        table->setType("float");
        table->setNumberOfRows(nsph);
        for (const auto &blockName : gasBlocks) table->addCol(blockName);
        memTables.push_back(table);
    }

    if (nsph > 0 && gasParticles) {
        float *bufferBlock = static_cast<float *>(std::malloc(sizeof(float) * nsph));
        if (!bufferBlock) {
            std::cerr << "Malloc Error for bufferBlock (gas)" << std::endl;
            return 1;
        }

        for (int elem = 0; elem < 12; ++elem) {
            for (int part = 0; part < nsph; ++part) {
                bufferBlock[part] = gasParticles[part * 12 + elem];
            }

            if (useMemory) {
                unsigned int colId = memTables[0]->getColId(gasBlocks[elem]);
                if (colId == static_cast<unsigned int>(-1)) {
                    std::cerr << "Invalid column: " << gasBlocks[elem] << std::endl;
                    continue;
                }

                unsigned int colList[1] = { colId };
                float *dataPtrs[1] = { bufferBlock };

                unsigned long long globalRowStart = 0;
                unsigned long long globalRowEnd   = static_cast<unsigned long long>(nsph) - 1;
                memTables[0]->putColumn(colList, 1, globalRowStart, globalRowEnd, dataPtrs);
            } else {
                outfile.write(reinterpret_cast<char *>(bufferBlock), sizeof(float) * nsph);
            }
        }

        std::free(bufferBlock);
    }

    if (!useMemory) {
        outfile.close();
        pathHeader = pathFileOut + "GAS" + ".bin";
        makeHeader(nsph, pathHeader, gasBlocks, m_cellSize, m_cellComp, m_volumeOrTable);
    }

    if (gasParticles) std::free(gasParticles);

    // -----------------------------
    // DARK PARTICLES (9 floats)
    // -----------------------------
    float *darkParticles = nullptr;
    if (ndark > 0) {
        darkParticles = static_cast<float *>(std::malloc(sizeof(float) * 9 * ndark));
        if (!darkParticles) {
            std::clog << "Malloc Error for darkParticles" << std::endl;
            return 1;
        }

        for (int i = 0; i < ndark; ++i) {
            for (int k = 0; k < 9; ++k) {
                darkParticles[i * 9 + k] = readFloatBE(inFile);
            }
        }
    }

    std::vector<std::string> darkBlocks;
    darkBlocks.push_back("MASS");   //0
    darkBlocks.push_back("POS_X");  //1
    darkBlocks.push_back("POS_Y");  //2
    darkBlocks.push_back("POS_Z");  //3
    darkBlocks.push_back("VEL_X");  //4
    darkBlocks.push_back("VEL_y");  //5
    darkBlocks.push_back("VEL_Z");  //6
    darkBlocks.push_back("EPS");    //7
    darkBlocks.push_back("PHI");    //8

    if (!useMemory) {
        outfile.open((pathFileOut + "DARK" + ".bin").c_str(), std::ofstream::binary);
        if (!outfile) {
            std::cerr << "Failed to open DARK.bin for writing" << std::endl;
            return 1;
        }
    } else {
        VSTable *table = new VSTableMem();
        table->setType("float");
        table->setNumberOfRows(ndark);
        for (const auto &blockName : darkBlocks) table->addCol(blockName);
        memTables.push_back(table);
    }

    if (ndark > 0 && darkParticles) {
        float *bufferBlock2 = static_cast<float *>(std::malloc(sizeof(float) * ndark));
        if (!bufferBlock2) {
            std::cerr << "Malloc Error for bufferBlock (dark)" << std::endl;
            return 1;
        }

        for (int elem = 0; elem < 9; ++elem) {
            for (int part = 0; part < ndark; ++part) {
            bufferBlock2[part] = darkParticles[part * 9 + elem];
        }

        if (useMemory) {
            unsigned int colId = memTables[1]->getColId(darkBlocks[elem]);
            if (colId == static_cast<unsigned int>(-1)) {
                std::cerr << "Invalid column: " << darkBlocks[elem] << std::endl;
                continue;
            }

                unsigned int colList[1] = { colId };
                float *dataPtrs[1] = { bufferBlock2 };
                unsigned long long globalRowStart = 0;
                unsigned long long globalRowEnd   = static_cast<unsigned long long>(ndark) - 1;
                memTables[1]->putColumn(colList, 1, globalRowStart, globalRowEnd, dataPtrs);
            } else {
                outfile.write(reinterpret_cast<char *>(bufferBlock2), sizeof(float) * ndark);
            }
        }

        std::free(bufferBlock2);
    }

    if (!useMemory) {
        outfile.close();
        pathHeader = pathFileOut + "DARK" + ".bin";
        makeHeader(ndark, pathHeader, darkBlocks, m_cellSize, m_cellComp, m_volumeOrTable);
    }

    if (darkParticles) std::free(darkParticles);

    // -----------------------------
    // STAR PARTICLES (11 floats)
    // -----------------------------
    float *starParticles = nullptr;
    if (nstar > 0) {
        starParticles = static_cast<float *>(std::malloc(sizeof(float) * 11 * nstar));
        if (!starParticles) {
            std::clog << "Malloc Error for starParticles" << std::endl;
            return 1;
        }
        for (int i = 0; i < nstar; ++i) {
            for (int k = 0; k < 11; ++k) {
                starParticles[i * 11 + k] = readFloatBE(inFile);
            }
        }
    }

    std::vector<std::string> starBlocks;
    starBlocks.push_back("MASS");   //0
    starBlocks.push_back("POS_X");  //1
    starBlocks.push_back("POS_Y");  //2
    starBlocks.push_back("POS_Z");  //3
    starBlocks.push_back("VEL_X");  //4
    starBlocks.push_back("VEL_y");  //5
    starBlocks.push_back("VEL_Z");  //6
    starBlocks.push_back("METALS"); //7
    starBlocks.push_back("TFORM");  //8
    starBlocks.push_back("EPS");    //9
    starBlocks.push_back("PHI");    //10

    if (!useMemory) {
        outfile.open((pathFileOut + "STAR" + ".bin").c_str(), std::ofstream::binary);
        if (!outfile) {
            std::cerr << "Failed to open STAR.bin for writing" << std::endl;
            return 1;
        }
    } else {
        VSTable *table = new VSTableMem();
        table->setType("float");
        table->setNumberOfRows(nstar);
        for (const auto &blockName : starBlocks) table->addCol(blockName);
        memTables.push_back(table); // index 2 for STAR
    }

    if (nstar > 0 && starParticles) {
        float *bufferBlock3 = static_cast<float *>(std::malloc(sizeof(float) * nstar));
        if (!bufferBlock3){
            std::cerr << "Malloc Error for bufferBlock3 (star)" << std::endl;
            return 1;
        }

        for (int elem = 0; elem < 11; ++elem) {
            for (int part = 0; part < nstar; ++part) {
                bufferBlock3[part] = starParticles[part * 11 + elem];
            }

            if (useMemory) {
                unsigned int colId = memTables[2]->getColId(starBlocks[elem]);
                if (colId == static_cast<unsigned int>(-1)) {
                    std::cerr << "Invalid column: " << starBlocks[elem] << std::endl;
                    continue;
                }

                unsigned int colList[1] = { colId };
                float *dataPtrs[1] = { bufferBlock3 };

                unsigned long long globalRowStart = 0;
                unsigned long long globalRowEnd   = static_cast<unsigned long long>(nstar) - 1;
                memTables[2]->putColumn(colList, 1, globalRowStart, globalRowEnd, dataPtrs);
            } else {
                outfile.write(reinterpret_cast<char *>(bufferBlock3), sizeof(float) * nstar);
            }
        }

        std::free(bufferBlock3);
    }

    if (!useMemory) {
        outfile.close();
        pathHeader = pathFileOut + "STAR" + ".bin";
        makeHeader(nstar, pathHeader, starBlocks, m_cellSize, m_cellComp, m_volumeOrTable);
    }

    if (starParticles) std::free(starParticles);

    inFile.close();

    return 0;
}

ChangaSource::~ChangaSource() {
}

ChangaSource::ChangaSource() {
}
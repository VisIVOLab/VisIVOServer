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
// #include "VisIVOImporterConfigure.h"
#include "changasource.h"

#include "mpi.h"
#include "visivoutils.h"
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <fstream>
#include <iostream>
#include <mpio.h>
#include <omp.h>
#include <stdexcept>
#include <unistd.h>
#include <vector>
#include <vtkIOStream.h>

namespace {

static inline int32_t readIntBE(std::istream &in) {
  unsigned char b[4];
  in.read(reinterpret_cast<char *>(b), 4);
  if (!in)
    throw std::runtime_error("readIntBE: failed to read 4 bytes");
  uint32_t v = (uint32_t(b[0]) << 24) | (uint32_t(b[1]) << 16) |
               (uint32_t(b[2]) << 8) | uint32_t(b[3]);
  return int32_t(v);
}

static inline float readFloatBE(std::istream &in) {
  unsigned char b[4];
  in.read(reinterpret_cast<char *>(b), 4);
  if (!in)
    throw std::runtime_error("readFloatBE: failed to read 4 bytes");
  uint32_t v = (uint32_t(b[0]) << 24) | (uint32_t(b[1]) << 16) |
               (uint32_t(b[2]) << 8) | uint32_t(b[3]);
  float f;
  std::memcpy(&f, &v, 4);
  return f;
}

static inline double readDoubleBE(std::istream &in) {
  unsigned char b[8];
  in.read(reinterpret_cast<char *>(b), 8);
  if (!in)
    throw std::runtime_error("readDoubleBE: failed to read 8 bytes");
  uint64_t v = (uint64_t(b[0]) << 56) | (uint64_t(b[1]) << 48) |
               (uint64_t(b[2]) << 40) | (uint64_t(b[3]) << 32) |
               (uint64_t(b[4]) << 24) | (uint64_t(b[5]) << 16) |
               (uint64_t(b[6]) << 8) | uint64_t(b[7]);
  double d;
  std::memcpy(&d, &v, 8);
  return d;
}

} // namespace

//---------------------------------------------------------------------
int ChangaSource::readHeader() {
  std::string fileName = m_pointsFileName;
  std::ifstream inFile(fileName, std::ios::binary);
  if (!inFile) {
    std::cerr << "Error while opening file in readHeader: " << fileName
              << std::endl;
    return -1;
  }

  double time = readDoubleBE(inFile);
  int32_t nbodies = readIntBE(inFile);
  int32_t ndim = readIntBE(inFile);
  int32_t nsph_ = readIntBE(inFile);
  int32_t ndark_ = readIntBE(inFile);
  int32_t nstar_ = readIntBE(inFile);
  int32_t pad = readIntBE(inFile);

  this->nsph = nsph_;
  this->ndark = ndark_;
  this->nstar = nstar_;

  inFile.close();

  return 0;
}

float *ChangaSource::readParticles(Particle particleType) {
  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int particlesNumber;

  if (rank == 0) {
    switch (particleType) {
    case GAS:
      particlesNumber = this->nsph;
      break;
    case DARK:
      particlesNumber = this->ndark;
      break;
    case STAR:
      particlesNumber = this->nstar;
      break;
    default:
      break;
    }
  }

  int particleFields;
  std::streamoff headerSize = 8 + 6 * 4;

  switch (particleType) {
  case GAS:
    particleFields = 12;
    break;
  case DARK:
    particleFields = 9;
    headerSize += this->nsph * 12 * sizeof(float);
    break;
  case STAR:
    particleFields = 11;
    headerSize +=
        this->nsph * 12 * sizeof(float) + this->ndark * 9 * sizeof(float);
    break;
  default:
    break;
  }

  int *particlesPerRank = nullptr;
  int *displacements = nullptr;
  int *particlesVarsPerRank = nullptr;
  int *displacementsVars = nullptr;

  if (rank == 0) {
    particlesPerRank = static_cast<int *>(std::malloc(sizeof(int) * size));
    displacements = static_cast<int *>(std::malloc(sizeof(int) * size));
    particlesVarsPerRank = static_cast<int *>(std::malloc(sizeof(int) * size));
    displacementsVars = static_cast<int *>(std::malloc(sizeof(int) * size));
    displacements[0] = 0;
    displacementsVars[0] = 0;

    int base = particlesNumber / size;
    int remainder = particlesNumber % size;

    for (int i = 0; i < size; i++) {
      particlesPerRank[i] = base;
      particlesVarsPerRank[i] = base * particleFields;
      if (i < remainder) {
        particlesPerRank[i] += 1;
        particlesVarsPerRank[i] += particleFields;
      }
      if (i > 0) {
        displacements[i] = displacements[i - 1] + particlesPerRank[i - 1];
        displacementsVars[i] =
            displacementsVars[i - 1] + particlesVarsPerRank[i - 1];
      }
    }
  }

  int *localNumParticlesBuffer = static_cast<int *>(std::malloc(sizeof(int)));
  int *localDisplacementBuffer = static_cast<int *>(std::malloc(sizeof(int)));

  MPI_Scatter(particlesPerRank, 1, MPI_INT, localNumParticlesBuffer, 1, MPI_INT,
              0, MPI_COMM_WORLD);
  MPI_Scatter(displacements, 1, MPI_INT, localDisplacementBuffer, 1, MPI_INT, 0,
              MPI_COMM_WORLD);

  int localNumParticles = localNumParticlesBuffer[0];
  int localDisplacement = localDisplacementBuffer[0];

  float *localGasParticles =
      (float *)std::malloc(sizeof(float) * particleFields * localNumParticles);

  // Each process opens the same file
  std::ifstream inFile(m_pointsFileName, std::ios::binary);
  if (!inFile) {
    std::cerr << "Error while opening file in readData: " << m_pointsFileName
              << std::endl;
    return nullptr;
  }

  // Skip header
  inFile.seekg(headerSize, std::ios::beg);
  if (!inFile) {
    std::cerr << "Failed to seek to particle data";
    return nullptr;
  }

  // Offsets to local displacement
  inFile.seekg(localDisplacement * particleFields * sizeof(float),
               std::ios_base::cur);

  for (int i = 0; i < localNumParticles; ++i) {
    for (int k = 0; k < particleFields; ++k) {
      localGasParticles[i * particleFields + k] = readFloatBE(inFile);
    }
  }

  float *particles = nullptr;

  if (rank == 0) {
    if (particlesNumber > 0) {
      particles = static_cast<float *>(
          std::malloc(sizeof(float) * particleFields * particlesNumber));
      if (!particles) {
        std::clog << "Malloc Error for particles" << std::endl;
        return nullptr;
      }
    }
  }

  MPI_Gatherv(localGasParticles, localNumParticles * particleFields, MPI_FLOAT,
              particles, particlesVarsPerRank, displacementsVars, MPI_FLOAT, 0,
              MPI_COMM_WORLD);

  free(localNumParticlesBuffer);
  localNumParticlesBuffer = nullptr;
  free(localDisplacementBuffer);
  localDisplacementBuffer = nullptr;
  free(localGasParticles);
  localGasParticles = nullptr;
  if (rank == 0) {
    free(particlesPerRank);
    particlesPerRank = nullptr;
    free(particlesVarsPerRank);
    particlesVarsPerRank = nullptr;
    free(displacements);
    displacements = nullptr;
    free(displacementsVars);
    displacementsVars = nullptr;
  }

  inFile.close();

  if (rank == 0)
    return particles;

  return nullptr;
}

int ChangaSource::writeParticles(Particle particleType, float *particles) {
  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Status status;

  int particlesNumber;
  int particleId;
  int particleFields;
  std::string particleStartPath;

  switch (particleType) {
  case GAS:
    particleId = 0;
    particleFields = 12;
    particleStartPath = "GAS";
    particlesNumber = this->nsph;
    break;
  case DARK:
    particleId = 1;
    particleFields = 9;
    particleStartPath = "DARK";
    particlesNumber = this->ndark;
    break;
  case STAR:
    particleId = 2;
    particleFields = 11;
    particleStartPath = "STAR";
    particlesNumber = this->nstar;
    break;
  default:
    break;
  }

  int *particlesPerRank = nullptr;
  int *displacements = nullptr;
  int *particlesVarsPerRank = nullptr;
  int *displacementsVars = nullptr;

  if (rank == 0) {
    particlesPerRank = static_cast<int *>(std::malloc(sizeof(int) * size));
    displacements = static_cast<int *>(std::malloc(sizeof(int) * size));
    particlesVarsPerRank = static_cast<int *>(std::malloc(sizeof(int) * size));
    displacementsVars = static_cast<int *>(std::malloc(sizeof(int) * size));
    displacements[0] = 0;
    displacementsVars[0] = 0;

    int base = particlesNumber / size;
    int remainder = particlesNumber % size;

    for (int i = 0; i < size; i++) {
      particlesPerRank[i] = base;
      particlesVarsPerRank[i] = base * particleFields;
      if (i < remainder) {
        particlesPerRank[i] += 1;
        particlesVarsPerRank[i] += particleFields;
      }
      if (i > 0) {
        displacements[i] = displacements[i - 1] + particlesPerRank[i - 1];
        displacementsVars[i] =
            displacementsVars[i - 1] + particlesVarsPerRank[i - 1];
      }
    }

    std::cout << "particles number: " << particlesNumber << endl;
    std::cout << "base: " << base << endl;
    std::cout << "remainder: " << remainder << endl;
    for (int i = 0; i < size; i++) {
      std::cout << "rank: " << i << " particles: " << particlesPerRank[i]
                << " displacement: " << displacements[i] << endl;
    }
  }

  int *localNumParticlesBuffer = static_cast<int *>(std::malloc(sizeof(int)));
  int *localDisplacementBuffer = static_cast<int *>(std::malloc(sizeof(int)));

  MPI_Scatter(particlesPerRank, 1, MPI_INT, localNumParticlesBuffer, 1, MPI_INT,
              0, MPI_COMM_WORLD);
  MPI_Scatter(displacements, 1, MPI_INT, localDisplacementBuffer, 1, MPI_INT, 0,
              MPI_COMM_WORLD);

  int localNumParticles = localNumParticlesBuffer[0];
  int localDisplacement = localDisplacementBuffer[0];

  float *localParticles =
      (float *)malloc(localNumParticles * particleFields * sizeof(float));

  MPI_Scatterv(particles, particlesVarsPerRank, displacementsVars, MPI_FLOAT,
               localParticles, localNumParticles * particleFields, MPI_FLOAT, 0,
               MPI_COMM_WORLD);

  std::cout << "rank: " << rank << " local num: " << localNumParticles
            << " local displacement: " << localDisplacement
            << " total number: " << particlesNumber << endl;

  int idx = m_pointsBinaryName.rfind('.');
  std::string pathFileIn = m_pointsBinaryName;
  if (idx != std::string::npos) {
    pathFileIn.erase(idx); // remove extension
  }

  std::string pathFileOut = pathFileIn;
  std::string pathHeader;

  std::vector<std::string> blocks;

  switch (particleType) {
  case GAS:
    blocks.push_back("MASS");
    blocks.push_back("POS_X");
    blocks.push_back("POS_Y");
    blocks.push_back("POS_Z");
    blocks.push_back("VEL_X");
    blocks.push_back("VEL_y");
    blocks.push_back("VEL_Z");
    blocks.push_back("RHO");
    blocks.push_back("TEMP");
    blocks.push_back("EPS");
    blocks.push_back("METALS");
    blocks.push_back("PHI");
    break;
  case DARK:
    blocks.push_back("MASS");
    blocks.push_back("POS_X");
    blocks.push_back("POS_Y");
    blocks.push_back("POS_Z");
    blocks.push_back("VEL_X");
    blocks.push_back("VEL_Y");
    blocks.push_back("VEL_Z");
    blocks.push_back("EPS");
    blocks.push_back("PHI");
    break;
  case STAR:
    blocks.push_back("MASS");
    blocks.push_back("POS_X");
    blocks.push_back("POS_Y");
    blocks.push_back("POS_Z");
    blocks.push_back("VEL_X");
    blocks.push_back("VEL_Y");
    blocks.push_back("VEL_Z");
    blocks.push_back("METALS");
    blocks.push_back("TFORM");
    blocks.push_back("EPS");
    blocks.push_back("PHI");
    break;
  default:
    break;
  }

  MPI_Offset fileOffset;
  int tableOffset;
  MPI_File fh;
  int amode;
  amode = MPI_MODE_CREATE | MPI_MODE_RDWR;
  int code;

  if (!useMemory) {
    code = MPI_File_open(MPI_COMM_WORLD,
                         (pathFileOut + particleStartPath + ".bin").c_str(),
                         amode, MPI_INFO_NULL, &fh);
    if (code != MPI_SUCCESS) {
      std::cerr << "Failed to open " << particleStartPath << ".bin for writing"
                << std::endl;
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
      return 1;
    }

  } else {
    VSTable *table = new VSTableMem();
    table->setType("float");
    table->setNumberOfRows(localNumParticles);
    for (const auto &blockName : blocks)
      table->addCol(blockName);
    tableOffset = rank + (particleId * 3) + particleId;
    memTables.insert(memTables.begin() + tableOffset, table);
  }

  if (localNumParticles > 0 && localParticles) {
    float *bufferBlock =
        static_cast<float *>(std::malloc(sizeof(float) * localNumParticles));
    if (!bufferBlock) {
      std::cerr << "Malloc Error for bufferBlock" << std::endl;
      return 1;
    }

    for (int field = 0; field < particleFields; ++field) {
      for (int particle = 0; particle < localNumParticles; ++particle) {
        bufferBlock[particle] =
            localParticles[particle * particleFields + field];
      }

      if (useMemory) {
        unsigned int colId = memTables[tableOffset]->getColId(blocks[field]);
        if (colId == static_cast<unsigned int>(-1)) {
          std::cerr << "Invalid column: " << blocks[field] << std::endl;
          continue;
        }

        unsigned int colList[1] = {colId};
        float *dataPtrs[1] = {bufferBlock};

        unsigned long long globalRowStart = 0;
        unsigned long long globalRowEnd =
            static_cast<unsigned long long>(localNumParticles) - 1;
        memTables[tableOffset]->putColumn(colList, 1, globalRowStart,
                                          globalRowEnd, dataPtrs);
      } else {
        fileOffset =
            (field * particlesNumber + localDisplacement) * sizeof(float);
        std::cout << "rank: " << rank << " file offset: " << fileOffset << endl;

        MPI_File_write_at(fh, fileOffset, bufferBlock, localNumParticles,
                          MPI_FLOAT, &status);
      }
    }

    std::free(bufferBlock);
    bufferBlock = NULL;
  }

  if (!useMemory) {
    MPI_File_close(&fh);
    if (rank == 0) {
      pathHeader = pathFileOut + particleStartPath + ".bin";
      makeHeader(particlesNumber, pathHeader, blocks, m_cellSize, m_cellComp,
                 m_volumeOrTable);
    }
  }

  free(localNumParticlesBuffer);
  localNumParticlesBuffer = nullptr;
  free(localDisplacementBuffer);
  localDisplacementBuffer = nullptr;
  free(localParticles);
  localParticles = nullptr;
  if (rank == 0) {
    free(particlesPerRank);
    particlesPerRank = nullptr;
    free(particlesVarsPerRank);
    particlesVarsPerRank = nullptr;
    free(displacements);
    displacements = nullptr;
    free(displacementsVars);
    displacementsVars = nullptr;
  }

  return 0;
}

int ChangaSource::readData() {
  MPI_Init(nullptr, nullptr);
  int rank, size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0) {
    memTables.reserve(size * 3);
  }

  float *gasParticles = readParticles(GAS);
  writeParticles(GAS, gasParticles);

  // float *darkParticles = readParticles(DARK);
  // writeParticles(DARK, darkParticles);

  // float *starParticles = readParticles(STAR);
  // writeParticles(STAR, starParticles);

  if (rank == 0) {
    free(gasParticles);
    // free(darkParticles);
    // free(starParticles);
    gasParticles = NULL;
    // darkParticles = NULL;
    // starParticles = NULL;
  }

  MPI_Finalize();

  return 0;
}

ChangaSource::~ChangaSource() {}

ChangaSource::ChangaSource() {}
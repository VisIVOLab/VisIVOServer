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
#include <cstddef>
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
  int32_t nsph = readIntBE(inFile);
  int32_t ndark = readIntBE(inFile);
  int32_t nstar = readIntBE(inFile);
  int32_t pad = readIntBE(inFile);

  this->nsph = nsph;
  this->ndark = ndark;
  this->nstar = nstar;

  inFile.close();

  return 0;
}

/**
 * @brief  Distributes the information containing the number of particles and
 * the displacement to each process.
 *
 *
 * @return a vector containing structures which contains said information for
 * each type of particle.
 */
std::vector<mpiProcessInfo> ChangaSource::distributeInfo() {
  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  std::vector<mpiProcessInfo> info;

  int **particlesPerRankPerParticle = NULL;
  int **displacementPerRankPerParticle = NULL;
  int numbersOfParticles[] = {this->nsph, this->ndark, this->nstar};

  if (rank == 0) {
    particlesPerRankPerParticle =
        static_cast<int **>(malloc(sizeof(int *) * this->typesOfParticle));
    displacementPerRankPerParticle =
        static_cast<int **>(malloc(sizeof(int *) * this->typesOfParticle));
    for (int i = 0; i < this->typesOfParticle; i++) {
      particlesPerRankPerParticle[i] =
          static_cast<int *>(malloc(sizeof(int) * size));
      displacementPerRankPerParticle[i] =
          static_cast<int *>(malloc(sizeof(int) * size));

      displacementPerRankPerParticle[i][0] = 0;
      for (int j = 0; j < size; j++) {
        int base = numbersOfParticles[i] / size;
        int remainder = numbersOfParticles[i] % size;

        particlesPerRankPerParticle[i][j] = base;
        if (j < remainder) {
          particlesPerRankPerParticle[i][j] += 1;
        }
        if (j > 0) {
          displacementPerRankPerParticle[i][j] =
              displacementPerRankPerParticle[i][j - 1] +
              particlesPerRankPerParticle[i][j - 1];
        }
      }
    }
  }

  int localNumParticlesBuffer;
  int localDisplacementBuffer;

  for (int i = 0; i < this->typesOfParticle; i++) {
    MPI_Scatter(rank == 0 ? particlesPerRankPerParticle[i] : NULL, 1, MPI_INT,
                &localNumParticlesBuffer, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(rank == 0 ? displacementPerRankPerParticle[i] : NULL, 1,
                MPI_INT, &localDisplacementBuffer, 1, MPI_INT, 0,
                MPI_COMM_WORLD);
    info.push_back({localNumParticlesBuffer, localDisplacementBuffer});
  }

  if (rank == 0) {
    for (int i = 0; i < this->typesOfParticle; i++) {
      free(particlesPerRankPerParticle[i]);
      free(displacementPerRankPerParticle[i]);
      particlesPerRankPerParticle[i] = NULL;
      displacementPerRankPerParticle[i] = NULL;
    }
    free(particlesPerRankPerParticle);
    free(displacementPerRankPerParticle);
    particlesPerRankPerParticle = NULL;
    displacementPerRankPerParticle = NULL;
  }

  return info;
}

/**
 * @brief Elaborates all the particles within the standard files.
 *
 * This function chunk-reads the information regarding the particles within the
 * standard files using MPI and OpenMP for different tasks inside the function.
 * Specifically, MPI is used for reading and writing, OpenMP is used for the
 * endianness swap and the columnization.
 *
 * @param info the vector containing the structs that contains the start
 * particle and the displacement.
 * @param particleType an enum that specifies the type of the particle to work
 * with.
 *
 * @return an int = 0 if successful, 1 otherwise.
 */
int ChangaSource::elaborateParticles(std::vector<mpiProcessInfo> info,
                                     Particle particleType) {
  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int particleFields;
  int particlesNumber;
  std::string particleStartPath;
  std::streamoff headerSize = 8 + 6 * 4;

  switch (particleType) {
  case GAS:
    particleFields = 12;
    particleStartPath = "GAS";
    particlesNumber = this->nsph;
    break;
  case DARK:
    particleFields = 9;
    particleStartPath = "DARK";
    particlesNumber = this->ndark;
    headerSize += this->nsph * 12 * sizeof(float);
    break;
  case STAR:
    particleFields = 11;
    particleStartPath = "STAR";
    particlesNumber = this->nstar;
    headerSize +=
        this->nsph * 12 * sizeof(float) + this->ndark * 9 * sizeof(float);
    break;
  default:
    break;
  }

  MPI_File readFileHandle;
  int readFileAmode = MPI_MODE_RDONLY;
  int code;

  code = MPI_File_open(MPI_COMM_WORLD, m_pointsFileName.c_str(), readFileAmode,
                       MPI_INFO_NULL, &readFileHandle);
  if (code != MPI_SUCCESS) {
    std::cerr << "Failure in opening the file.\n";
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    return 1;
  }

  code = MPI_File_seek(readFileHandle,
                       headerSize + info[particleType].localDisplacement *
                                        particleFields * sizeof(float),
                       MPI_SEEK_SET);
  if (code != MPI_SUCCESS) {
    std::cerr << "Failed to seek." << std::endl;
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    return 1;
  }

  int tableOffset;
  MPI_File writeFileHandle;
  int writeFileAmode = MPI_MODE_CREATE | MPI_MODE_WRONLY;
  int idx = m_pointsBinaryName.rfind('.');
  std::string pathFileIn = m_pointsBinaryName;
  if (idx != std::string::npos)
    pathFileIn.erase(idx);

  std::string pathFileOut = pathFileIn;
  std::string pathHeader;

  std::vector<std::string> blocks;
  if (rank == 0) {
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
      blocks.push_back("VEL_y");
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
      blocks.push_back("VEL_y");
      blocks.push_back("VEL_Z");
      blocks.push_back("METALS");
      blocks.push_back("TFORM");
      blocks.push_back("EPS");
      blocks.push_back("PHI");
      break;
    default:
      break;
    }
    pathHeader = pathFileOut + particleStartPath + ".bin";
    makeHeader(particlesNumber, pathHeader, blocks, m_cellSize, m_cellComp,
               m_volumeOrTable);
  }

  if (useMemory) {
    std::cout << "rank: " << rank << " about to create table" << endl;
    VSTable *table = new VSTableMem();
    table->setType("float");
    table->setNumberOfRows(info[particleType].localNumParticles);
    for (const auto &blockName : blocks)
      table->addCol(blockName);
    tableOffset = rank * this->typesOfParticle + particleType;
    std::cout << "rank: " << rank << " table offset: " << tableOffset << endl;
    memTables.insert(memTables.begin() + tableOffset, table);
  }

  else {
    code = MPI_File_open(MPI_COMM_WORLD,
                         (pathFileOut + particleStartPath + ".bin").c_str(),
                         writeFileAmode, MPI_INFO_NULL, &writeFileHandle);
    if (code != MPI_SUCCESS) {
      std::cerr << "Failed to open " << particleStartPath << ".bin for writing."
                << std::endl;
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
      return 1;
    }
  }

  const size_t RAW_CHUNK_SIZE = 16 * 1024 * 1024; // 16 MB
  const size_t CHUNK_SIZE = RAW_CHUNK_SIZE / (particleFields * sizeof(float)) *
                            (particleFields * sizeof(float)); // floor value

  MPI_Status status;
  MPI_Offset writeFileOffset;

  size_t bytesLeft =
      info[particleType].localNumParticles * particleFields * sizeof(float);
  size_t bytesToRead;
  size_t structsLeft;
  int particlesRead;
  int particlesProcessedSoFar = 0;
  uint8_t *rawBuffer = new uint8_t[CHUNK_SIZE];
  float *buffer = new float[CHUNK_SIZE / sizeof(float)];
  float *columnizedBuffer = new float[CHUNK_SIZE / sizeof(float)];
  float f;
  uint32_t v;

  // chunked read
  while (bytesLeft > 0) {
    if (bytesLeft >= CHUNK_SIZE) {
      bytesToRead = CHUNK_SIZE;
    } else {
      structsLeft = bytesLeft / (particleFields * sizeof(float));
      bytesToRead = structsLeft * (particleFields * sizeof(float));
      if (bytesToRead == 0)
        break;
    }
    MPI_File_read(readFileHandle, rawBuffer, bytesToRead, MPI_UINT8_T, &status);
    bytesLeft -= bytesToRead;
    particlesRead = bytesToRead / (particleFields * sizeof(float));

    // endianness swap
#pragma omp parallel for
    for (int i = 0; i < bytesToRead; i += 4) {
      v = (uint32_t(rawBuffer[i]) << 24) | (uint32_t(rawBuffer[i + 1]) << 16) |
          (uint32_t(rawBuffer[i + 2]) << 8) | (uint32_t(rawBuffer[i + 3]));
      memcpy(&f, &v, 4);
      buffer[i / 4] = f;
    }

    // columnization
    for (int field = 0; field < particleFields; ++field) {
#pragma omp parallel for
      for (int particle = 0; particle < particlesRead; ++particle) {
        columnizedBuffer[particle] = buffer[particle * particleFields + field];
      }

      if (useMemory) {
        unsigned int colId = memTables[tableOffset]->getColId(blocks[field]);
        if (colId == static_cast<unsigned int>(-1)) {
          std::cerr << "Invalid column: " << blocks[field] << std::endl;
          continue;
        }

        unsigned int colList[1] = {colId};
        float *dataPtrs[1] = {columnizedBuffer};

        unsigned long long globalRowStart = 0;
        unsigned long long globalRowEnd =
            static_cast<unsigned long long>(particlesRead) - 1;

        memTables[tableOffset]->putColumn(colList, 1, globalRowStart,
                                          globalRowEnd, dataPtrs);
      } else {
        writeFileOffset =
            sizeof(float) *
            (field * particlesNumber +
             (info[particleType].localDisplacement + particlesProcessedSoFar));

        MPI_File_write_at(writeFileHandle, writeFileOffset, columnizedBuffer,
                          particlesRead, MPI_FLOAT, &status);
      }
    }
    particlesProcessedSoFar += particlesRead;
  }

  delete[] rawBuffer;
  delete[] buffer;
  delete[] columnizedBuffer;

  code = MPI_File_close(&readFileHandle);
  if (code != MPI_SUCCESS) {
    std::cerr << "Failed to close read file." << std::endl;
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    return 1;
  }

  code = MPI_File_close(&writeFileHandle);
  if (code != MPI_SUCCESS) {
    std::cerr << "Failed to close write file." << std::endl;
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    return 1;
  }

  return 0;
}

int ChangaSource::readData() {
  int provided;
  MPI_Init_thread(nullptr, nullptr, MPI_THREAD_FUNNELED, &provided);
  int rank, size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0) {
    memTables.reserve(size * 3);
    std::cout << m_changaDen << endl;
    std::cout << useMemory << endl;
  }

  std::vector<mpiProcessInfo> info = distributeInfo();

  elaborateParticles(info, GAS);
  elaborateParticles(info, DARK);
  elaborateParticles(info, STAR);

  MPI_Finalize();

  return 0;
}

ChangaSource::~ChangaSource() {}

ChangaSource::ChangaSource() {}
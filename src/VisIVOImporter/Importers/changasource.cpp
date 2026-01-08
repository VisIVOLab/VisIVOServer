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
#include "changasource.h"

#include "mpi.h"
#include "visivoutils.h"
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <fstream>
#include <ios>
#include <iostream>
#include <omp.h>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <vector>
#include <vtkIOStream.h>

#if __has_include(<mpio.h>)
#include <mpio.h>
#endif

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

  std::vector<mpiProcessInfo> info(this->typesOfParticle);

  std::vector<mpiProcessInfo> allData;

  MPI_Datatype MPI_mpiProcessInfo;
  MPI_Type_contiguous(2, MPI_INT, &MPI_mpiProcessInfo);
  MPI_Type_commit(&MPI_mpiProcessInfo);

  std::vector<int> numbersOfParticles = {this->nsph, this->ndark, this->nstar};

  if (rank == 0) {
    allData.reserve(size * this->typesOfParticle);
    int base, remainder;
    int count, displacement;
    for (int i = 0; i < size; i++) {
      if (i == 0)
        displacement = 0;
      for (int j = 0; j < this->typesOfParticle; j++) {
        base = numbersOfParticles[j] / size;
        remainder = numbersOfParticles[j] % size;
        count = base + (i < remainder ? 1 : 0);
        if (i > 0) {
          displacement =
              allData[(i - 1) * this->typesOfParticle + j].localDisplacement +
              allData[(i - 1) * this->typesOfParticle + j].localNumParticles;
        }
        allData.push_back({count, displacement});
      }
    }
  }

  MPI_Scatter(rank == 0 ? allData.data() : nullptr, typesOfParticle,
              MPI_mpiProcessInfo, info.data(), typesOfParticle,
              MPI_mpiProcessInfo, 0, MPI_COMM_WORLD);
  MPI_Type_free(&MPI_mpiProcessInfo);

  return info;
}

/**
 * @brief  Elaborates the additional information basing on the variables.
 *
 *
 * @return a vector containing structures which contains said information.
 */
std::vector<additionalMpiInfo> ChangaSource::elaborateAdditionalInfo() {
  std::vector<additionalMpiInfo> additionalInfo;
  if (m_changaDen) {
    std::string denExt = ".den";
    std::string fileName = m_pointsFileName + denExt;
    additionalInfo.push_back({1, {"DENSITY"}, fileName, 4});
  }

  return additionalInfo;
}

/**
 * @brief  Populates a vector of strings basing on the type of particle and on
 * additional fields provided.
 *
 * @param particleType an enum that specifies the type of the particle to
 * work with.
 * @param additionalInfo an array of additionaMpiInfo structs to extract
 * additional field names from.
 * @param numberOfFields the total number of fields.
 *
 * @return a vector of string containing the names of the fields.
 */
std::vector<std::string>
ChangaSource::populateBlocks(Particle particleType,
                             std::vector<additionalMpiInfo> additionalInfo,
                             int numberOfFields) {
  std::vector<std::string> blocks;
  blocks.reserve(numberOfFields);
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

  for (int i = 0; i < additionalInfo.size(); i++) {
    for (int j = 0; j < additionalInfo[i].fieldsName.size(); j++) {
      blocks.push_back(additionalInfo[i].fieldsName[j]);
    }
  }

  return blocks;
}

bool ChangaSource::readNextChunk(particleChunk &chunk,
                                 particleReadContext &ctx) {
  MPI_Status status;
  size_t structsLeft;
  size_t totalBytesToRead;

  if (ctx.bytesLeftPerFile[0] >= ctx.chunkSizes[0]) {
    chunk.bytesToReadPerFile[0] = ctx.chunkSizes[0];
  } else {
    structsLeft =
        ctx.bytesLeftPerFile[0] / (chunk.particleFields * sizeof(float));
    chunk.bytesToReadPerFile[0] =
        structsLeft * chunk.particleFields * sizeof(float);
    if (chunk.bytesToReadPerFile[0] == 0)
      return false;
  }
  totalBytesToRead = chunk.bytesToReadPerFile[0];

  MPI_File_read(ctx.readFileHandle, chunk.rawFileBuffers[0],
                chunk.bytesToReadPerFile[0], MPI_UINT8_T, &status);
  chunk.particlesRead =
      chunk.bytesToReadPerFile[0] / (chunk.particleFields * sizeof(float));

  for (int i = 0; i < chunk.additionalInfo.size(); i++) {
    chunk.bytesToReadPerFile[i + 1] = chunk.particlesRead *
                                      chunk.additionalInfo[i].particleFields *
                                      sizeof(float);
    totalBytesToRead += chunk.bytesToReadPerFile[i + 1];
    MPI_File_read(ctx.additionalReadFilesHandles[i],
                  chunk.rawFileBuffers[i + 1], chunk.bytesToReadPerFile[i + 1],
                  MPI_UINT8_T, &status);
  }

  ctx.totalBytesLeft -=
      chunk.particlesRead *
      (chunk.particleFields + chunk.additionalParticleFields) * sizeof(float);

  ctx.bytesLeftPerFile[0] -= chunk.bytesToReadPerFile[0];
  for (int i = 0; i < chunk.additionalInfo.size(); i++) {
    ctx.bytesLeftPerFile[i + 1] -= chunk.bytesToReadPerFile[i + 1];
  }

  return true;
}

void ChangaSource::elaborateChunk(particleChunk &chunk) {
  swapEndianness(chunk.rawFileBuffers[0], chunk.bytesToReadPerFile[0],
                 chunk.fileBuffers[0]);
  for (int i = 0; i < chunk.additionalInfo.size(); i++) {
    swapEndianness(chunk.rawFileBuffers[i + 1], chunk.bytesToReadPerFile[i + 1],
                   chunk.fileBuffers[i + 1]);
  }
}

void ChangaSource::writeChunk(particleChunk &chunk, particleWriteContext &ctx) {
  int localField = 0;
  MPI_Offset writeFileOffset;
  MPI_Status status;

  for (int field = 0;
       field < chunk.particleFields + chunk.additionalParticleFields; ++field) {
    if (field < chunk.particleFields) {
      chunk.currentFile = 0;
      localField = field;
    } else {
      for (int j = 0; j < chunk.additionalInfo.size() + 1; j++) {
        if (field < (ctx.fieldOffsets[j] + chunk.particleFields)) {
          chunk.currentFile = j;
          break;
        }
      }
      localField = field - (ctx.fieldOffsets[chunk.currentFile - 1] +
                            chunk.particleFields);
    }

    columnizeBuffer(
        chunk.columnizedBuffers[chunk.currentFile],
        chunk.fileBuffers[chunk.currentFile], chunk.particlesRead,
        field < chunk.particleFields
            ? chunk.particleFields
            : chunk.additionalInfo[chunk.currentFile - 1].particleFields,
        localField);

    if (useMemory) {
      unsigned int colId =
          memTables[ctx.tableOffset]->getColId(ctx.blocks[field]);
      if (colId == static_cast<unsigned int>(-1)) {
        std::cerr << "Invalid column: " << ctx.blocks[field] << std::endl;
        continue;
      }

      unsigned int colList[1] = {colId};
      float *dataPtrs[1] = {chunk.columnizedBuffers[chunk.currentFile]};

      unsigned long long globalRowStart = ctx.particlesProcessedSoFar;
      unsigned long long globalRowEnd =
          static_cast<unsigned long long>(ctx.particlesProcessedSoFar +
                                          chunk.particlesRead) -
          1;

      memTables[ctx.tableOffset]->putColumn(colList, 1, globalRowStart,
                                            globalRowEnd, dataPtrs);
    } else {
      writeFileOffset =
          sizeof(float) * (field * ctx.particlesNumber +
                           (ctx.info[ctx.particleType].localDisplacement +
                            ctx.particlesProcessedSoFar));

      MPI_File_write_at(ctx.writeFileHandle, writeFileOffset,
                        chunk.columnizedBuffers[chunk.currentFile],
                        chunk.particlesRead, MPI_FLOAT, &status);
    }
  }
  ctx.particlesProcessedSoFar += chunk.particlesRead;
}

/**
 * @brief Swaps the endianness of an array of bytes representing big-endian
 * float values.
 *
 * This function reads a buffer of bytes, converts each value to the host
 * endianness, and stores the resulting floats into a separate output
 * buffer.
 *
 * The conversion is performed in parallel using OpenMP.
 *
 * @param buffer Pointer to the input byte buffer containing big-endian
 * float representations.
 * @param bytes Total number of bytes in the input buffer.
 * @param newBuffer Pointer to the output buffer where the converted floats
 * will be stored.
 */
void ChangaSource::swapEndianness(uint8_t *buffer, size_t bytes,
                                  float *newBuffer) {
  float f;
  uint32_t v;
#pragma omp parallel for
  for (int i = 0; i < bytes; i += 4) {
    v = (uint32_t(buffer[i]) << 24) | (uint32_t(buffer[i + 1]) << 16) |
        (uint32_t(buffer[i + 2]) << 8) | (uint32_t(buffer[i + 3]));
    memcpy(&f, &v, 4);
    newBuffer[i / 4] = f;
  }
}

/**
 * @brief Extracts a column from a row-major buffer.
 *
 * This function copies the values of a specified column from a buffer
 * stored in row-major order into a linear buffer.
 * The operation is parallelized using OpenMP.
 *
 * @param columnBuffer   Output buffer containing the extracted column.
 * @param originalBuffer Input buffer storing the data in row-major order.
 * @param length         Number of elements to extract.
 * @param rows           Number of rows in the original buffer.
 * @param currentColumn  Index of the column to extract.
 */
void ChangaSource::columnizeBuffer(float *columnBuffer, float *originalBuffer,
                                   int length, int rows, int currentColumn) {
#pragma omp parallel for
  for (int column = 0; column < length; ++column) {
    columnBuffer[column] = originalBuffer[column * (rows) + currentColumn];
  }
}

/**
 * @brief Closes all MPI file handles used by the ChangaSource importer.
 *
 * This function closes the main MPI read and write file handles and MPI read
 * every file handle specified in the additionalInfo vector.
 *
 * @param writeFileHandle Pointer to the MPI file handle used for writing.
 * @param readFileHandle  Pointer to the MPI file handle used for reading.
 * @param additionalInfo  Vector containing metadata for any additional MPI read
 * files. Its size determines how many additional file handles must be closed.
 * @param additionalReadFilesHandles Array of MPI file handles corresponding to
 * the entries in additionalInfo.
 *
 * @return an int = 0 if successful, 1 otherwise.
 */
int ChangaSource::closeFiles(MPI_File *writeFileHandle,
                             MPI_File *readFileHandle,
                             std::vector<MPI_File> additionalReadFilesHandles) {
  int code;
  code = MPI_File_close(readFileHandle);
  if (code != MPI_SUCCESS) {
    std::cerr << "Failed to close read file." << std::endl;
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    return 1;
  }

  code = MPI_File_close(writeFileHandle);
  if (code != MPI_SUCCESS) {
    std::cerr << "Failed to close write file." << std::endl;
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    return 1;
  }

  for (int i = 0; i < additionalReadFilesHandles.size(); i++) {
    code = MPI_File_close(&additionalReadFilesHandles[i]);
    if (code != MPI_SUCCESS) {
      std::cerr << "Failed to close read file." << std::endl;
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
      return 1;
    }
  }
  return 0;
}

/**
 * @brief Processes all the particles within the standard files.
 *
 * This function chunk-reads, -elaborates, and -writes the information regarding
 * the particles within the standard extension files using MPI and OpenMP for
 * different tasks inside the function.
 *
 * @param particleType an enum that specifies the type of the particle to
 * work with.
 * @param info the vector containing the structs that contain the start
 * particle and the displacement.
 * @param additionalInfo the vector containing the structs that contain
 * the additional info to work with.
 *
 * @return an int = 0 if successful, 1 otherwise.
 */
int ChangaSource::processParticles(
    Particle particleType, std::vector<mpiProcessInfo> info,
    std::vector<additionalMpiInfo> additionalInfo) {
  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int particleFields;
  int additionalParticleFields = 0;
  int particlesNumber;
  std::string particleStartPath;
  std::streamoff headerSize = 8 + 6 * 4;
  int startParticle = 0;

  switch (particleType) {
  case GAS:
    particleFields = 12;
    particleStartPath = "GAS";
    particlesNumber = this->nsph;
    break;
  case DARK:
    particleFields = 9;
    particleStartPath = "DARK";
    startParticle = this->nsph;
    particlesNumber = this->ndark;
    headerSize += this->nsph * 12 * sizeof(float);
    break;
  case STAR:
    particleFields = 11;
    particleStartPath = "STAR";
    startParticle = this->nsph + this->ndark;
    particlesNumber = this->nstar;
    headerSize +=
        this->nsph * 12 * sizeof(float) + this->ndark * 9 * sizeof(float);
    break;
  default:
    break;
  }

  MPI_File readFileHandle;
  std::vector<MPI_File> additionalReadFilesHandles(additionalInfo.size());
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

  for (int i = 0; i < additionalInfo.size(); i++) {
    additionalParticleFields += additionalInfo[i].particleFields;
    code = MPI_File_open(MPI_COMM_WORLD, additionalInfo[i].filename.c_str(),
                         readFileAmode, MPI_INFO_NULL,
                         &additionalReadFilesHandles[i]);
    if (code != MPI_SUCCESS) {
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
      return 1;
    }

    headerSize =
        additionalInfo[i].headerSize +
        (startParticle * additionalInfo[i].particleFields * sizeof(float));
    headerSize += info[particleType].localDisplacement *
                  additionalInfo[i].particleFields * sizeof(float);

    code =
        MPI_File_seek(additionalReadFilesHandles[i], headerSize, MPI_SEEK_SET);
    if (code != MPI_SUCCESS) {
      std::cerr << "Failed to seek." << std::endl;
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
      return 1;
    }
  }

  std::vector<std::string> blocks = populateBlocks(
      particleType, additionalInfo, particleFields + additionalParticleFields);

  int tableOffset;
  MPI_File writeFileHandle;
  int writeFileAmode = MPI_MODE_CREATE | MPI_MODE_WRONLY;
  int idx = m_pointsBinaryName.rfind('.');
  std::string pathFileIn = m_pointsBinaryName;
  if (idx != std::string::npos)
    pathFileIn.erase(idx);

  std::string pathFileOut = pathFileIn;
  std::string pathHeader;

  if (rank == 0) {
    pathHeader = pathFileOut + particleStartPath + ".bin";
    makeHeader(particlesNumber, pathHeader, blocks, m_cellSize, m_cellComp,
               m_volumeOrTable);
  }

  if (useMemory) {
    VSTable *table = new VSTableMem();
    table->setType("float");
    table->setNumberOfRows(info[particleType].localNumParticles);
    for (const auto &blockName : blocks)
      table->addCol(blockName);
    tableOffset = rank * this->typesOfParticle + particleType;
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
  size_t totalChunkSize = 0;

  std::vector<std::size_t> chunkSizes(1 + additionalInfo.size());
  chunkSizes[0] = (RAW_CHUNK_SIZE / (particleFields * sizeof(float))) *
                  (particleFields * sizeof(float));
  totalChunkSize += chunkSizes[0];

  std::vector<std::size_t> bytesLeftPerFile(1 + additionalInfo.size());
  bytesLeftPerFile[0] =
      info[particleType].localNumParticles * particleFields * sizeof(float);

  std::vector<int> fieldOffsets(1 + additionalInfo.size());
  fieldOffsets[0] = 0;

  uint8_t **rawFileBuffers = new uint8_t *[1 + additionalInfo.size()];
  rawFileBuffers[0] = new uint8_t[chunkSizes[0]];

  float **fileBuffers = new float *[1 + additionalInfo.size()];
  fileBuffers[0] = new float[chunkSizes[0] / sizeof(float)];

  float **columnizedBuffers = new float *[1 + additionalInfo.size()];
  columnizedBuffers[0] = new float[chunkSizes[0] / sizeof(float)];

  for (int i = 0; i < additionalInfo.size(); i++) {
    chunkSizes[i + 1] = (chunkSizes[0] / (particleFields * sizeof(float))) *
                        additionalInfo[i].particleFields * sizeof(float);
    totalChunkSize += chunkSizes[i + 1];
    bytesLeftPerFile[i + 1] = info[particleType].localNumParticles *
                              additionalInfo[i].particleFields * sizeof(float);
    fieldOffsets[i + 1] = fieldOffsets[i] + additionalInfo[i].particleFields;
    rawFileBuffers[i + 1] = new uint8_t[chunkSizes[i + 1]];
    fileBuffers[i + 1] = new float[chunkSizes[i + 1] / sizeof(float)];
    columnizedBuffers[i + 1] = new float[chunkSizes[i + 1] / sizeof(float)];
  }

  std::vector<size_t> bytesToReadPerFile(1 + additionalInfo.size());
  size_t totalBytesLeft = info[particleType].localNumParticles *
                          (particleFields + additionalParticleFields) *
                          sizeof(float);
  size_t totalBytesToRead;
  size_t structsLeft;
  int particlesRead = 0;
  int particlesProcessedSoFar = 0;
  unsigned int currentFile = 0;

  particleChunk chunk = {
      particleFields,     additionalParticleFields, particlesRead,
      bytesToReadPerFile, rawFileBuffers,           fileBuffers,
      columnizedBuffers,  additionalInfo,           currentFile};
  particleReadContext readCtx = {totalBytesLeft, bytesLeftPerFile, chunkSizes,
                                 readFileHandle, additionalReadFilesHandles};
  particleWriteContext writeCtx = {particlesNumber,
                                   particlesProcessedSoFar,
                                   particleType,
                                   fieldOffsets,
                                   info,
                                   writeFileHandle,
                                   blocks,
                                   tableOffset};

  while (readNextChunk(chunk, readCtx)) {
    elaborateChunk(chunk);
    writeChunk(chunk, writeCtx);
  }

  closeFiles(&writeFileHandle, &readFileHandle, additionalReadFilesHandles);

  delete[] rawFileBuffers[0];
  delete[] fileBuffers[0];
  delete[] columnizedBuffers[0];
  for (int i = 0; i < additionalInfo.size(); i++) {
    delete[] rawFileBuffers[i + 1];
    delete[] columnizedBuffers[i + 1];
    delete[] fileBuffers[i + 1];
  }
  delete[] rawFileBuffers;
  delete[] columnizedBuffers;
  delete[] fileBuffers;

  return 0;
}

int ChangaSource::readData() {
  int provided;
  MPI_Init_thread(nullptr, nullptr, MPI_THREAD_FUNNELED, &provided);
  int rank, size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0)
    memTables.reserve(size * 3);

  std::vector<mpiProcessInfo> info = distributeInfo();
  std::vector<additionalMpiInfo> additionalInfo = elaborateAdditionalInfo();

  processParticles(GAS, info, additionalInfo);
  processParticles(DARK, info, additionalInfo);
  processParticles(STAR, info, additionalInfo);

  MPI_Finalize();

  return 0;
}

ChangaSource::~ChangaSource() {}

ChangaSource::ChangaSource() {}
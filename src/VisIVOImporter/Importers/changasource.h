/***************************************************************************
 *   Copyright (C) 2008 by Gabriella Caniglia, Roberto Munzone *
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

#ifndef CHANGASOURCE_H
#define CHANGASOURCE_H

#include "abstractsource.h"

#include <mpi.h>
#include <string>
#include <vector>

enum Particle { GAS = 0, DARK = 1, STAR = 2 };

typedef struct {
  int particleFields;
  std::vector<std::string> fieldsName;
  std::string filename;
  std::streamoff headerSize;
} additionalMpiInfo;

typedef struct {
  int localNumParticles;
  int localDisplacement;
} mpiProcessInfo;

typedef struct {
  int particleFields;
  int additionalParticleFields;
  int particlesRead;
  std::vector<size_t> bytesToReadPerFile;
  uint8_t **rawFileBuffers;
  float **fileBuffers;
  float **columnizedBuffers;
  std::vector<additionalMpiInfo> additionalInfo;
  unsigned int currentFile;
} particleChunk;

typedef struct {
  size_t totalBytesLeft;
  std::vector<size_t> bytesLeftPerFile;
  std::vector<std::size_t> chunkSizes;
  MPI_File readFileHandle;
  std::vector<MPI_File> additionalReadFilesHandles;
} particleReadContext;

typedef struct {
  int particlesNumber;
  int particlesProcessedSoFar;
  Particle particleType;
  std::vector<int> fieldOffsets;
  std::vector<mpiProcessInfo> info;
  MPI_File writeFileHandle;
  std::vector<std::string> blocks;
  int tableOffset;
} particleWriteContext;

class ChangaSource : public AbstractSource {
public:
  int readHeader();
  int readData();
  ~ChangaSource();
  ChangaSource();
  std::vector<mpiProcessInfo> distributeInfo();
  std::vector<additionalMpiInfo> elaborateAdditionalInfo();
  std::vector<std::string>
  populateBlocks(Particle particleType,
                 std::vector<additionalMpiInfo> additionalInfo,
                 int numberOfFields);
  bool readNextChunk(particleChunk &chunk, particleReadContext &ctx);
  void elaborateChunk(particleChunk &chunk);
  void writeChunk(particleChunk &chunk, particleWriteContext &ctx);
  void swapEndianness(uint8_t *buffer, size_t bytes, float *newBuffer);
  void columnizeBuffer(float *columnBuffer, float *originalBuffer, int length,
                       int rows, int currentRow);
  int closeFiles(MPI_File *writeFileHandle, MPI_File *readFileHandle,
                 std::vector<MPI_File>);
  int processParticles(Particle particleType, std::vector<mpiProcessInfo> info,
                       std::vector<additionalMpiInfo> additionalInfo);

private:
  const int typesOfParticle = 3;
  int nsph;
  int ndark;
  int nstar;
};

#endif
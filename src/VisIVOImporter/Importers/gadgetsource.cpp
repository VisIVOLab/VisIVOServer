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
#include "gadgetsource.h"

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
#include <unistd.h>

#ifndef off64_t
#define off64_t off_t
#endif

bool isNumeric(std::string const &str)
{
  auto it = str.begin();
  while (it != str.end() && std::isdigit(*it)) {
    it++;
  }
  return !str.empty() && it == str.end();
}

//---------------------------------------------------------------------
int GadgetSource::readHeader()
//---------------------------------------------------------------------

{
  char dummy[4]; 
//   char checkType[4];
  int provided;
  MPI_Init_thread(0, 0, MPI_THREAD_MULTIPLE, &provided);
  int type=0; unsigned int Npart=0;
  std::string systemEndianism;
  numBlock=0;
  bool needSwap=false;
#ifdef VSBIGENDIAN
  systemEndianism="big";
#else
  systemEndianism="little";
#endif
  if((m_endian=="b" || m_endian=="big") && systemEndianism=="little")
    needSwap=true;
  if((m_endian=="l" || m_endian=="little") && systemEndianism=="big")
    needSwap=true;
	
  std::string fileName = m_pointsFileName.c_str();
  std::ifstream inFile;

  inFile.open(fileName, std::ios::binary);
  if (!inFile)
  {
    std::clog << fileName << std::endl;
    std::cerr<<"Error while opening File"<<std::endl;
    return -1;
  }
  inFile.read((char *)(dummy), 4*sizeof(char)); //!*** IMPORTANT NOT REMOVE ***//
  inFile.read((char *)(tmpType), 4*sizeof(char));   //!*** IMPORTANT NOT REMOVE ***//

  tagType=tmpType;
  checkType.push_back(tagType);
  if (iCompare (checkType[0].c_str(), "HEAD") == 0)  //!read header Type2 else Type1
  {
    headerType2 curHead;
    inFile.seekg(4, std::ios::beg);
    
    inFile.read((char *)(&curHead), 256);   //!*** COPY DATA in STRUCT m_pHeaderType2 ***//
    m_snapformat = 2;
    if (needSwap)
      swapHeaderType2();
    
    int nFiles = curHead.num_files[0];
    if(nFiles > 1){
      fileName = m_pointsFileName.c_str();
      std::filesystem::path p(fileName);
      if(p.extension().generic_string().size() > 0 && isNumeric(p.extension().generic_string().substr(1))){
        for(int j = 0; j < p.extension().generic_string().substr(1).size(); j++){
          fileName.pop_back();
        }
        
        if(checkMultipleFiles(nFiles, fileName)){
          numFiles = nFiles;
          if (readMultipleHeaders(nFiles, fileName, needSwap) == 0){
            
            inFile.close();
            MPI_Finalize();
            return -1;
          }
          else return 0;
        }
        else {
          inFile.close();         
          MPI_Finalize();
          return -1;
        }
      } 
      else{
        inFile.close();
        std::cerr << "Expecting multiple files" << std::endl;
      }
    }
    m_pHeaderType2.push_back(curHead);
    numFiles = 1;
  }
  else  //Type1 
  {
    inFile.seekg(0, std::ios::beg);	
    inFile.read((char *)(&m_pHeaderType1), 268);
    m_snapformat = 1;
    
    if (needSwap)
      swapHeaderType1();
    
    if (m_pHeaderType1.sizeFirstBlock[0] != (3*Npart*sizeof(float)))
    {
      std::cerr<<"The size File is not than expected"<<std::endl;
      return -1;
      
    }
      
  }
  inFile.close();
  return 0;
}

int GadgetSource::readMultipleHeaders(int nFiles, std::string fileName, bool needSwap)
{
  long long int sum = 0; 
  for(int i = 0; i < nFiles; i++){ 
    std::ifstream inFile; 
    headerType2 curHead;
    inFile.open(fileName + std::to_string(i), std::ios::binary);
    
    inFile.seekg(4, std::ios::beg);
    
    inFile.read((char *)(&curHead), 256);   //!*** COPY DATA in STRUCT m_pHeaderType2 ***//
    if (needSwap)
      swapHeaderType2();
    m_pHeaderType2.push_back(curHead);
    sum += curHead.npart[1];
    inFile.close();
  }
  return 1;
}

int GadgetSource::checkMultipleFiles(int nFiles, std::string fName){  
  std::ifstream inFile;
  
  for(int i = 0; i < nFiles; i++){
    inFile.open(fName + std::to_string(i), std::ios::binary);
    if (!inFile)
    {
      std::cerr<<"Error while opening multiple files, read only input file"<<std::endl;
      return 0;  
    }
    inFile.close();
    
  }
  
  return 1; 
}

bool GadgetSource::determineEndianism() {
    std::string systemEndianism;
#ifdef VSBIGENDIAN
    systemEndianism = "big";
#else
    systemEndianism = "little";
#endif
    return ((m_endian == "b" || m_endian == "big") && systemEndianism == "little") ||
           ((m_endian == "l" || m_endian == "little") && systemEndianism == "big");
}

int GadgetSource::readData()
{
  
  char dummy[4]; 
  //int inFile[numFiles];
  std::string systemEndianism;
  needSwap = determineEndianism();

  int type=0;
  unsigned int i=0;
  unsigned int j=0;

  // Get the number of processes
  MPI_Comm_size(MPI_COMM_WORLD, &processingConfig.num_proc);
  MPI_Comm_rank (MPI_COMM_WORLD, &processingConfig.proc_id);

  char tagTmp[5]="";
  std::string tag; 
//   const char point = '.';
  int idx = m_pointsBinaryName.rfind('.');
  std::string pathFileOut = m_pointsBinaryName.erase(idx, idx+4); 
   
  std::string bin = ".bin";
  std::string X = "_X"; std::string Y = "_Y"; std::string Z = "_Z";

  std::vector<std::string> tagTypeForNameFile = {"GAS", "HALO", "DISK", "BULGE", "STARS", "BNDRY"};

  initializeParticleCounts(npartTotal64);

  generateMap(blockData.mapBlockNamesToFields, blockNamesToCompare, blockNamesToFields);
  generateMap(blockData.mapBlockSize, blockNamesToCompare, blockSize);

  //std::vector<std::vector<std::string> > namesFields;
  std::vector<std::string> tmpNamesFields;
  std::string pathHeader = "";  int KK=0;
  
  std::string fileName = processFileName(m_pointsFileName);
  //create list of blocks to read
  blockData.listOfBlocks = extractBlockList(fileName, needSwap);

  std::vector<std::string> blockNames;

  // Declare typePosition and fileStartPosition using vectors
  //std::vector<std::vector<int>> typePosition;
  //std::vector<std::vector<long long>> fileStartPosition;

  // Compute type positions
  computeTypePositions();

  // Compute file start positions
  computeFileStartPositions(numFiles, m_pHeaderType2);
  
  int nTotalBlocks = blockData.listOfBlocks.size() * numFiles;
  double t1, t2; 
  
  bool alreadyOpen = false;
  int totBlocks = blockData.listOfBlocks.size();

  extractHeaderFields(blockData.listOfBlocks, blockData.mapBlockSize, blockData.mapBlockNamesToFields, fieldTypeNames);
  openInputFiles(fileName, numFiles);

  //int outFileBin[6]; 
  openOutputFiles(pathFileOut, tagTypeForNameFile, bin);

  processBlocksParallel(totBlocks);

  for(int nFile = 0; nFile < numFiles; nFile++)close(fileData.inFile[nFile]);
  for(type = 0; type < 6; type++)close(fileData.outFileBin[type]);
   
  //WRITE HEADER FILES
  if (processingConfig.proc_id == 0) {
    writeHeaderFiles(pathFileOut, fieldTypeNames, tagTypeForNameFile);
  }  
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  if (!memTables.empty()) {
    VSTable* table = memTables[0];
  }
    
  return 0;
}

void GadgetSource::processBlocksParallel(int totBlocks) {
    #pragma omp parallel for collapse(2)
    for (int nBlock = 0; nBlock < totBlocks; nBlock++) {
        for (int nFile = processingConfig.proc_id; nFile < numFiles; nFile += processingConfig.num_proc) {
            long long unsigned int offset = findBlockOffset(fileData.inFile[nFile], blockData.listOfBlocks[nBlock], needSwap);
            unsigned long long minPart[6];

            computeMinPart(nFile, minPart);

            for (int type = 0; type < 6; type++) {
              unsigned long long int startingPoint = 0;
              if(m_pHeaderType2[nFile].npart[type] == 0) continue;
              for(int i = 0; i<type; i++) 
              {
                startingPoint+=m_pHeaderType2[nFile].npart[i];
              }
              processParticle(type, nBlock, nFile, minPart, offset+(startingPoint*sizeof(float)*blockData.mapBlockSize.at(blockData.listOfBlocks[nBlock])/*blockData.mapBlockSize.at(blockData.listOfBlocks[nBlock]*/));
            }
        }
    }

}

void GadgetSource::processParticle(int type, int nBlock, int nFile, 
                     unsigned long long* minPart, long long unsigned int offset) {
    // Skip invalid particle types
    if (!isValidParticleType(type, nFile, nBlock, blockData.listOfBlocks, blockData.mapBlockNamesToFields))return;

    unsigned long long chunk = minPart[type];
    unsigned long long n = m_pHeaderType2[nFile].npart[type]/ chunk;
    unsigned long long Resto = m_pHeaderType2[nFile].npart[type] - (chunk * n);
    //unsigned long long pToStart = blockData.typePosition[nBlock][type] * npartTotal64[type] + blockData.fileStartPosition[type][nFile];
    unsigned long long newOffSet = offset;

    int blockSize = blockData.mapBlockSize.at(blockData.listOfBlocks[nBlock]);
    float* bufferBlock = new float[blockSize * chunk];

    std::vector<float*> buffers;
    allocateBuffers(blockSize, chunk, buffers);

    std::string blockName;

    // Process chunks
    if(useMemory) blockName = blockData.listOfBlocks[nBlock];
    for (int k = 0; k < n; k++) {
      processChunk(fileData.inFile[nFile], bufferBlock, buffers, blockSize, chunk, newOffSet);

      newOffSet += blockSize * chunk *sizeof(float);      
      if(useMemory) writeChunkData(fileData.outFileBin[type], buffers, chunk, chunk, nBlock, nFile, k, blockSize, type, blockName);
      else writeChunkData(fileData.outFileBin[type], buffers, chunk, chunk, nBlock, nFile, k, blockSize, type, blockName);

    }
    // Handle Resto
    if (Resto > 0) {
      processChunk(fileData.inFile[nFile], bufferBlock, buffers, blockSize, Resto, offset);
      
      if(useMemory)writeChunkData(fileData.outFileBin[type], buffers, chunk, Resto, nBlock, nFile, n, blockSize, type, blockName);
      else writeChunkData(fileData.outFileBin[type], buffers, chunk, Resto, nBlock, nFile, n, blockSize, type, blockName);

    }

    delete[] bufferBlock;
    for (int dim = 0; dim < blockSize; ++dim) {
        delete[] buffers[dim];
    }
}

void GadgetSource::allocateBuffers(int blockSize, unsigned long long chunk, std::vector<float*>& buffers) {
    buffers.resize(blockSize, nullptr);
    for (int dim = 0; dim < blockSize; ++dim) {
        buffers[dim] = new float[chunk];
    }
}

void GadgetSource::processChunk(int fileDescriptor, float* bufferBlock, std::vector<float*>& buffers, 
                  int blockSize, unsigned long long chunk, long long unsigned int offset) {
    // Read raw data
    pread(fileDescriptor, (char*)(bufferBlock), blockSize * chunk * sizeof(float), offset);
    
    // Process the data
    for (unsigned long long i = 0; i < chunk; i++) {
        for (int dim = 0; dim < blockSize; dim++) {
            buffers[dim][i] = needSwap ? floatSwap((char*)(&bufferBlock[blockSize * i + dim]))
                                       : bufferBlock[blockSize * i + dim];
        }
    }
    
}

void GadgetSource::writeChunkData(int outputFile, const std::vector<float*>& buffers, 
                    unsigned long long chunk, unsigned long long writeSize, 
                    int nBlock, int nFile,
                    unsigned long long chunkIndex, 
                    int blockSize, int type, std::string blockName) {
  if (useMemory) {
        for (int dim = 0; dim < blockSize; ++dim) {
              unsigned int colId;
              if (blockSize == 1) {

                  colId = memTables[type]->getColId(blockName);
              } else if (blockSize == 3) {
                  if (dim == 0) colId = memTables[type]->getColId(blockName + "_X");
                  else if (dim == 1) colId = memTables[type]->getColId(blockName + "_Y");
                  else colId = memTables[type]->getColId(blockName + "_Z");
              } else {
                  colId = memTables[type]->getColId(blockName + "_" + std::to_string(dim));
              }

              if (colId == (unsigned int)-1) {
                  std::cerr << "Invalid column name for block " << blockName << ", dim " << dim << std::endl;
                  continue;
              }
              unsigned int colList[1] = {colId};

              float* dataPtrs[1] = {buffers[dim]};
              unsigned long long int globalRowStart = blockData.fileStartPosition[type][nFile] + chunk * chunkIndex;
              unsigned long long int  globalRowEnd = globalRowStart + writeSize - 1;
              memTables[type]->putColumn(colList, 1, globalRowStart, globalRowEnd, dataPtrs);
          }
      }
  else{
    for (int dim = 0; dim < blockSize; ++dim) {
        unsigned long long pToStart = blockData.typePosition[nBlock][type] * npartTotal64[type] + blockData.fileStartPosition[type][nFile];
        unsigned long long pWrite = (pToStart * sizeof(float)) + 
                                    (chunkIndex * chunk * sizeof(float)) + 
                                    (npartTotal64[type] * dim * sizeof(float));
        ssize_t written = pwrite(outputFile, (char*)(buffers[dim]), writeSize * sizeof(float), pWrite);
    }
  }
}

bool GadgetSource::isValidParticleType(int type, int nFile, int nBlock,
                         const std::vector<std::string>& listOfBlocks,
                         const std::unordered_map<std::string, int>& mapBlockNamesToFields) {
    // Check if the particle type has no particles
    if (m_pHeaderType2[nFile].npart[type] == 0) return false;

    // Check if the block field is valid for this type
    if (!blocksFields[mapBlockNamesToFields.at(listOfBlocks[nBlock])][type]) return false;

    // "MASS" block should only be processed if mass[type] != 0
    if (iCompare(listOfBlocks[nBlock], "MASS") == 0 && m_pHeaderType2[nFile].mass[type] != 0)
        return false;

    return true;
}

void GadgetSource::computeMinPart(int fileIndex, unsigned long long minPart[6]) {
    unsigned int param = 1, esp = 32;
    unsigned long long maxULI = ldexp((float)param, esp);
    std::fill_n(minPart, 6, maxULI);
    for (int type = 0; type < 6; type++) {
        if (m_pHeaderType2[fileIndex].npart[type] != 0 && m_pHeaderType2[fileIndex].npart[type] <= 2500000)
            minPart[type] = m_pHeaderType2[fileIndex].npart[type];
        else
            minPart[type] = 2500000;
    }
}


long long GadgetSource::findBlockOffset(int fileDescriptor, const std::string& targetBlock, bool needSwap) {
    std::string tag;
    char tagTmp[5] = "";
    long long unsigned int offset = 0;
    int sizeBlock[1] = {0};

    while (iCompare(tag, targetBlock) != 0) {
        offset += 4;
        pread(fileDescriptor, (char *)(tagTmp), 4 * sizeof(char), offset);
        offset += 4;
        tag = strtok(tagTmp, " ");
        if (iCompare(tag, targetBlock) == 0) break;
        
        pread(fileDescriptor, (char *)(sizeBlock), sizeof(int), offset);
        offset += (4 + sizeBlock[0] + 4);
        
        if (needSwap) sizeBlock[0] = intSwap((char *)(&sizeBlock[0]));
    }

    return offset + 12;
}


std::string GadgetSource::processFileName(std::string s) {
    std::filesystem::path p(s);
    std::string fileName = s;

    if (p.has_extension()) {
        std::string ext = p.extension().string();
        
        if (!ext.empty() && ext[0] == '.' && isNumeric(ext.substr(1))) {  
            fileName.erase(fileName.size() - ext.size() + 1);  // Remove only numeric part of the extension
        }
    }
    return fileName;
}

void GadgetSource::extractHeaderFields(const std::vector<std::string>& listOfBlocks, 
                                       const std::unordered_map<std::string, int>& mapBlockSize,
                                       const std::unordered_map<std::string, int>& mapBlockNamesToFields,
                                       std::vector<std::vector<std::string>>& namesFields) {
    std::vector<std::string> tmpNamesFields;

    for (int type = 0; type < 6; type++) {
        if (m_pHeaderType2[0].npart[type] == 0) {  
            tmpNamesFields.clear();
            tmpNamesFields.push_back(" ");
            namesFields.push_back(tmpNamesFields);
            tmpNamesFields.clear();
        } else {
            for (const auto& block : listOfBlocks) {
                if (mapBlockSize.at(block) == 3 && blocksFields[mapBlockNamesToFields.at(block)][type]) {
                    tmpNamesFields.push_back(block + "_X");
                    tmpNamesFields.push_back(block + "_Y");
                    tmpNamesFields.push_back(block + "_Z");
                } else if (mapBlockSize.at(block) > 1 && blocksFields[mapBlockNamesToFields.at(block)][type]) {
                    for (int j = 0; j < mapBlockSize.at(block); j++) {
                        tmpNamesFields.push_back(block + "_" + std::to_string(j));
                    }
                } else if (mapBlockSize.at(block) == 1 && blocksFields[mapBlockNamesToFields.at(block)][type] &&
                           (iCompare(block, "MASS") != 0 || m_pHeaderType2[0].mass[type] == 0)) {
                    tmpNamesFields.push_back(block);
                }
            }
            namesFields.push_back(tmpNamesFields);
            tmpNamesFields.clear();
        }
    }
}

void GadgetSource::writeHeaderFiles(const std::string& pathFileOut, 
                                    const std::vector<std::vector<std::string>>& namesFields, 
                                    const std::vector<std::string>& tagTypeForNameFile) {
    std::string pathHeader;
    for (int type = 0; type < 6; type++) {
        if (npartTotal64[type] != 0) {
            for (const auto& field : namesFields[type]) {
                m_fieldsNames.push_back(field);
            }
            pathHeader = pathFileOut + tagTypeForNameFile[type] + ".bin";
            makeHeader(npartTotal64[type], pathHeader, m_fieldsNames, m_cellSize, m_cellComp, m_volumeOrTable);
            m_fieldsNames.clear();
            pathHeader.clear();
        }
    }
}


void GadgetSource::generateMap(std::unordered_map<std::string, int>& targetMap,
                               const std::vector<std::string>& keys,
                               const std::vector<int>& values) {
    for (size_t i = 0; i < keys.size(); ++i) {
        targetMap[keys[i]] = values[i];
    }
}

void GadgetSource::initializeParticleCounts(unsigned long long npartTotal64[6]) {
    for (int type = 0; type < 6; type++) {
        npartTotal64[type] = (static_cast<unsigned long long>(m_pHeaderType2[0].NallWH[type]) << 32) 
                             + m_pHeaderType2[0].npartTotal[type];
    }
}

void GadgetSource::computeTypePositions() {
    int numTypes = 6;
    blockData.typePosition.assign(blockData.listOfBlocks.size() + 1, std::vector<int>(numTypes, 0));

    for (int j = 0; j < numTypes; j++) {
        for (size_t i = 0; i < blockData.listOfBlocks.size(); i++) {
            if (iCompare(blockData.listOfBlocks[i], "MASS") != 0 || 
                (iCompare(blockData.listOfBlocks[i], "MASS") == 0 && m_pHeaderType2[0].mass[j] == 0)) {
                blockData.typePosition[i + 1][j] = blocksFields[blockData.mapBlockNamesToFields.at(blockData.listOfBlocks[i])][j] * 
                                         blockData.mapBlockSize.at(blockData.listOfBlocks[i]);
            } else {
                blockData.typePosition[i + 1][j] = 0;
            }
        }
    }

    for (int j = 0; j < numTypes; j++) {
        for (size_t i = 1; i < blockData.listOfBlocks.size(); i++) {
            blockData.typePosition[i][j] += blockData.typePosition[i - 1][j];
        }
    }
    
}

void GadgetSource::computeFileStartPositions(int numFiles, const std::vector<headerType2>& m_pHeaderType2) {
    int numTypes = 6;  
    blockData.fileStartPosition.assign(numTypes, std::vector<long long>(numFiles + 1, 0));
    
    for (int file = 1; file < numFiles; file++) {
        for (int type = 0; type < numTypes; type++) {
            blockData.fileStartPosition[type][file] = blockData.fileStartPosition[type][file - 1] + 
                                            static_cast<unsigned long long>(m_pHeaderType2[file - 1].npart[type]);
          }
    }
}

void GadgetSource::setNumFiles(int n) {
    if (n < 1) {
        //std::cerr << "Error: numFiles must be at least 1." << std::endl;
        return;
    }
    this->numFiles = n;
}


std::vector<std::string> GadgetSource::extractBlockList(const std::string& fileName, bool needSwap) {
    std::ifstream fStream(fileName + "0", std::ios::binary);
    if (!fStream.is_open()) {
        throw std::runtime_error("Error opening file: " + fileName + "0");
    }

    fStream.seekg(280, std::ios::beg);
    std::vector<std::string> listOfBlocks;
    char tagTmp[5] = "";
    std::string tag;

    while (fStream.peek() != EOF) {
        fStream.seekg(4, std::ios::cur);
        fStream.read(tagTmp, 4);
        tagTmp[4] = '\0'; 
        tag = strtok(tagTmp, " ");
        if (std::find(blockNamesToCompare.begin(), blockNamesToCompare.end(), tag) != blockNamesToCompare.end() &&
            iCompare(tag, "Zs") != 0) {
            if (m_fields.empty() || std::find(m_fields.begin(), m_fields.end(), tag) != m_fields.end()) {
                listOfBlocks.push_back(tag);
            }

        }

        fStream.read(reinterpret_cast<char*>(m_sizeBlock), sizeof(int));
        if (needSwap) {
            m_sizeBlock[0] = intSwap(reinterpret_cast<char*>(&m_sizeBlock[0]));
        }

        fStream.seekg(m_sizeBlock[0], std::ios::cur);
        fStream.read(reinterpret_cast<char*>(m_sizeBlock), sizeof(int));
    }

    fStream.close();
    return listOfBlocks;
}

void GadgetSource::openOutputFiles(const std::string& pathFileOut, 
                                   const std::vector<std::string>& tagTypeForNameFile, 
                                   const std::string& bin) {
    for (int type = 0; type < 6; type++) {
        fileData.outFileBin[type] = -1; 
        if(useMemory){
          VSTable* table = new VSTableMem();
          std::string nameFileBinOut = pathFileOut + tagTypeForNameFile[type] + bin;
          table->setLocator(nameFileBinOut);
          table->setType("float");
          table->setNumberOfRows(npartTotal64[type]);
          for(int i = 0; i < fieldTypeNames[type].size(); i++) table->addCol(fieldTypeNames[type][i]);
          memTables.push_back(table);
        }
        else{
          if (npartTotal64[type] != 0) {	
              std::string nameFileBinOut = pathFileOut + tagTypeForNameFile[type] + bin;
              fileData.outFileBin[type] = open(nameFileBinOut.c_str(), O_RDWR | O_CREAT | O_TRUNC, S_IRWXU);
              // Check for errors in file creation
              if (fileData.outFileBin[type] == -1) {
                  std::cerr << "Error: Failed to create output file " << nameFileBinOut << std::endl;
                  std::exit(EXIT_FAILURE);
              }
          }
        }
    }
}

void GadgetSource::openInputFiles(const std::string& fileName, int numFiles) {
    std::filesystem::path p(fileName);
    std::string baseFileName = fileName;
    if (numFiles > 1 && p.has_extension() && isNumeric(p.extension().string().substr(1))) {
        baseFileName = fileName.substr(0, fileName.rfind('.')); // Remove numerical extension
    }

    for (int i = 0; i < numFiles; i++) {
        std::string fullPath = baseFileName + std::to_string(i);
        fileData.inFile.push_back(open(fullPath.c_str(), O_RDONLY));
        if (fileData.inFile[i] == -1) {
            std::cerr << "Error opening input file: " << fullPath << std::endl;
            exit(EXIT_FAILURE);
        }
    }
}



//--------------------------------------
void GadgetSource::swapHeaderType2()
//--------------------------------------
{
  char first_fortran_legacy[4]; //m_snapformat=2 : you have to "jump" 4+2*4 =12 bytes and you'll find the the Block size
  m_pHeaderType2[0].boh[0] = intSwap((char*)(&m_pHeaderType2[0].boh[0]));
  m_pHeaderType2[0].nBlock[0] = intSwap((char*)(&m_pHeaderType2[0].nBlock[0]));
  m_pHeaderType2[0].size[0] = intSwap((char*)(&m_pHeaderType2[0].size[0]));
	  
  for (int i=0; i<6; i++)
    m_pHeaderType2[0].npart[i] = intSwap((char*)(&m_pHeaderType2[0].npart[i]));
	  
  for (int i=0; i<6; i++)
    m_pHeaderType2[0].mass[i] = intSwap((char*)(&m_pHeaderType2[0].mass[i]));
	  
  m_pHeaderType2[0].time[0] = doubleSwap((char*)(&m_pHeaderType2[0].time[0]));
  m_pHeaderType2[0].redshift[0] = doubleSwap((char*)(&m_pHeaderType2[0].redshift[0]));
  m_pHeaderType2[0].flag_sfr[0] = intSwap((char*)(&m_pHeaderType2[0].flag_sfr[0]));
  m_pHeaderType2[0].flag_feedback[0] = intSwap((char*)(&m_pHeaderType2[0].flag_feedback[0]));
	  
  for (int i=0; i<6; i++)
    m_pHeaderType2[0].npartTotal[i] = intSwap((char*)(&m_pHeaderType2[0].npartTotal[i]));
   
  m_pHeaderType2[0].foling[0] = intSwap((char*)(&m_pHeaderType2[0].foling[0]));
  m_pHeaderType2[0].num_files[0] = intSwap((char*)(&m_pHeaderType2[0].num_files[0]));
  m_pHeaderType2[0].BoxSize[0] = doubleSwap((char*)&(m_pHeaderType2[0].BoxSize[0]));
  m_pHeaderType2[0].Omega0[0] = doubleSwap((char*)(&m_pHeaderType2[0].Omega0[0]));
  m_pHeaderType2[0].OmegaLambda[0] = doubleSwap((char*)(&m_pHeaderType2[0].OmegaLambda[0]));
  m_pHeaderType2[0].HubbleParam[0] = doubleSwap((char*)(&m_pHeaderType2[0].HubbleParam[0]));
  m_pHeaderType2[0].FlagAge[0] = intSwap((char*)(&m_pHeaderType2[0].FlagAge[0]));
  m_pHeaderType2[0].FlagMetals[0] = intSwap((char*)(&m_pHeaderType2[0].FlagMetals[0]));
	  
  for (int i=0; i<6; i++)
    m_pHeaderType2[0].NallWH[i] = intSwap((char*)(&m_pHeaderType2[0].NallWH[i]));
   
  m_pHeaderType2[0].flag_entr_ics[0] = intSwap((char*)(&m_pHeaderType2[0].flag_entr_ics[0]));
   
  char fill[256- 6*sizeof(int)- 6*sizeof(double)- 2*sizeof(double)- 2*sizeof(int)- 6*sizeof(int)- 2*sizeof(int)- 
      4*sizeof(double)- 9*sizeof(int)]; /* fills to 256 Bytes */
  m_pHeaderType2[0].final_boh[0] = intSwap((char*)(&m_pHeaderType2[0].final_boh[0]));
  m_pHeaderType2[0].final_nBlock[0] = intSwap((char*)(&m_pHeaderType2[0].final_nBlock[0]));
  
  char tagFirstBlock[4];
  m_pHeaderType2[0].first_boh[0] = intSwap((char*)(&m_pHeaderType2[0].first_boh[0]));
  m_pHeaderType2[0].first_nBlock[0] = intSwap((char*)(&m_pHeaderType2[0].first_nBlock[0]));
  m_pHeaderType2[0].sizeFirstBlock[0] = intSwap((char*)(&m_pHeaderType2[0].sizeFirstBlock[0]));
}

  //--------------------------------------
  void GadgetSource::swapHeaderType1()
  //--------------------------------------
{	
  m_pHeaderType1.size[0] = intSwap((char*)(&m_pHeaderType1.size[0]));
	  
  for (int i=0; i<6; i++)
    m_pHeaderType1.npart[i] = intSwap((char*)(&m_pHeaderType1.npart[i]));
	  
  for (int i=0; i<6; i++)
    m_pHeaderType1.mass[i] = intSwap((char*)(&m_pHeaderType1.mass[i]));
	  
  m_pHeaderType1.time[0] = doubleSwap((char*)(&m_pHeaderType1.time[0]));
  m_pHeaderType1.redshift[0] = doubleSwap((char*)(&m_pHeaderType1.redshift[0]));
  m_pHeaderType1.flag_sfr[0] = intSwap((char*)(&m_pHeaderType1.flag_sfr[0]));
  m_pHeaderType1.flag_feedback[0] = intSwap((char*)(&m_pHeaderType1.flag_feedback[0]));
	  
  for (int i=0; i<6; i++)
    m_pHeaderType1.npartTotal[i] = intSwap((char*)(&m_pHeaderType1.npartTotal[i]));
   
  m_pHeaderType1.foling[0] = intSwap((char*)(&m_pHeaderType1.foling[0]));
  m_pHeaderType1.num_files[0] = intSwap((char*)(&m_pHeaderType1.num_files[0]));
  m_pHeaderType1.BoxSize[0] = doubleSwap((char*)(&m_pHeaderType1.BoxSize[0]));
  m_pHeaderType1.Omega0[0] = doubleSwap((char*)(&m_pHeaderType1.Omega0[0]));
  m_pHeaderType1.OmegaLambda[0] = doubleSwap((char*)(&m_pHeaderType1.OmegaLambda[0]));
  m_pHeaderType1.HubbleParam[0] = doubleSwap((char*)(&m_pHeaderType1.HubbleParam[0]));
  m_pHeaderType1.FlagAge[0] = intSwap((char*)(&m_pHeaderType1.FlagAge[0]));
  m_pHeaderType1.FlagMetals[0] = intSwap((char*)(&m_pHeaderType1.FlagMetals[0]));
	  
  for (int i=0; i<6; i++)
    m_pHeaderType1.NallWH[i] = intSwap((char*)(&m_pHeaderType1.NallWH[i]));
   
  m_pHeaderType1.flag_entr_ics[0] = intSwap((char*)(&m_pHeaderType1.flag_entr_ics[0]));
   
  char fill[256- 6*sizeof(int)- 6*sizeof(double)- 2*sizeof(double)- 2*sizeof(int)- 6*sizeof(int)- 2*sizeof(int)- 
      4*sizeof(double)- 9*sizeof(int)]; /* fills to 256 Bytes */
  m_pHeaderType1.sizeFirstBlock[0] = intSwap((char*)(&m_pHeaderType1.sizeFirstBlock[0]));
   
}

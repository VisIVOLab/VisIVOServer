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
#ifndef VSPOINTDISTRIBUTEOP_H
#define VSPOINTDISTRIBUTEOP_H
#include "vstableop.h"

/**
	@author Ugo Becciani <ugo.becciani@oact.inaf.it>
*/

struct GridHandle {
    unsigned long long int totEle = 0;
    unsigned long long startCounter = 0;
    unsigned long long fromRow = 0;
    unsigned long long toRow = 0;
};
class VSPointDistributeOp: public VSTableOp

{  
  static const unsigned int MAX_NUMBER_TO_REDUCE_ROW;
  static const unsigned int MIN_NUMBER_OF_ROW;
  unsigned int m_sampleDimensions[3];
  double m_modelBounds[6];
  double m_nullValue;
  double m_massUnity;
  unsigned int m_colList[3];	

  int m_splattedScalar;
  unsigned long long int m_numNewPts;
  bool m_useConstant;
  unsigned int m_nOfCol;
  unsigned int m_nOfRow;
  float **m_fArray;
  float **m_grid;
  float m_origin[3];
  float m_spacing[3];
  float m_constValue;
  bool m_executeDone;
  bool m_tsc;
  bool m_cic;
  bool m_ngp;
  bool setOrigin();
  bool setSpacing();
  unsigned long long int m_gridPts;
  bool m_OriginSet;
  bool m_SpacingSet;
  bool m_gridSpacing;
  bool m_avg;
  bool m_periodic;
  GridHandle gridHandle;

  bool allocateArray(int nField);
  bool allocateArray(int nField, bool isMP);
  bool computeModelBounds();
  int parsePointColumns(unsigned int colList[3]) ;
  void setAlgorithm();
  bool setGridResolution();
  bool parseFieldList(std::vector<int>& fieldList);
  bool allocateColumnList(unsigned int*& colList, const std::vector<int>& fieldList);
  bool initializeGrid(const std::vector<int>& fieldList, unsigned int* colList);
  bool processGridSpacing();
  std::string generateOutputFileName();
  void configureTableGrid(VSTable& tableGrid, const std::vector<int>& fieldList);
  bool initializeEmptyGrid(VSTable& tableGrid, const std::vector<int>& fieldList, unsigned int*& gridList);
  void computeCICWeights(float pos1, float pos2, float pos3, float weights[8]);
  bool processCIC(VSTable& tableGrid, unsigned int* gridList, int nOfField, 
                    std::vector<int>& fieldList, unsigned long long* gridIndex);
  void applyPeriodicBoundary_CIC(int& i1, int& i2, int& i3, int& i11, int& i21, int& i31);
  void applyPeriodicBoundary_TSC(float px[3]);
  int getLinearizedIndex(int x, int y, int z);
  float computeTSCWeight(float posX, float posY, float posZ, int gridX, int gridY, int gridZ);
  void processTSC(VSTable& tableGrid, unsigned int* gridList, int nOfField,
                  std::vector<int>& fieldList, unsigned long long* gridIndex);
  void computeNGPIndex(float px[3], int gridPos[3]);  
  void processNGP(VSTable& tableGrid, unsigned int* gridList, int nOfField,
                  std::vector<int>& fieldList, unsigned long long* gridIndex);

    public:
    VSPointDistributeOp();
    ~VSPointDistributeOp();
    void printHelp();
    bool execute();
    bool getOrigin(float *origin); 
    bool getSpacing(float *spacing); 
};

#endif
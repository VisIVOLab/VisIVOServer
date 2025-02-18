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


#ifndef SUBFINDSOURCE_H
#define SUBFINDSOURCE_H

#include "abstractsource.h"

#include <vector>
#include <string>


struct headerSubFind
{
   long long Ngroups[1];       // number of groups in this file
   long long Nsubhalos[1];     // number of subhalos in the present file
   long long Nids[1];          // number of particles in groups in this file
   long long NgroupsTot[1];    // total number of groups
   long long NsubhalosTot[1];  // total number of subhalos
   long long NidsTot[1];       // total number of particles in groups
   int       NumFiles[1];      // number of subfiles in this catalogue
   double    Time[1];          // output time/scale factor
   double    Redshift[1];      // output redshift
};

class SubFindSource : public AbstractSource
   
{
  public: //! Read the headerType2 file and set the basic table parameters
    int readHeader();
    int readData();
        
  private:
    std::vector <std::string> m_fieldsNames;   
    unsigned int      npart_total[6];
    int numFiles = 1;
    int m_nRows;
    char m_dataType, m_Endian;
    
     int m_snapformat;
     char tmpType[4]; int numBlock; int m_sizeBlock[1];
    
    std::vector<std::string> checkType;
    std::string tagType;
    	
    void swapHeaderType2();
    void swapHeaderType1();

    bool isNumeric(std::string const &str);
    
    std::vector<headerSubFind > m_pHeaderType2;
    int readMultipleHeaders(int, std::string, bool);
    int checkMultipleFiles(int, std::string);
    void updateNpart2(int);
        const std::vector<std::string> blockNamesToCompare = {
        "POS",  // 0
        "VEL",  // 1
        "ID",   // 2
        "MASS", // 3
        "U",    // 4
        "TEMP", // 5
        "RHO",  // 6
        "NE",   // 7
        "NH",   // 8
        "HSML", // 9
        "SFR",  // 10
        "AGE",  // 11
        "Z",    // 12
        "Zs",   // 13
        "iM",   // 14
        "ZAGE", // 15
        "ZALV", // 16
        "CLDX", // 17
        "TSTP", // 18
        "POT",  // 19
        "ACCE", // 20
        "ENDT", // 21
  //blockNamesToCompare.push_back("TSTP");//10
        "IDU",  // 22
        "HOTT", // 23
        "MHOT", // 24
        "MCLD", // 25
        "EHOT", // 26
        "MSF",  // 27
        "MFST", // 28
        "NMF",  // 29
        "EOUT", // 30
        "EREC", // 31
        "EOLD", // 32
        "TDYN", // 33
        "SFRo", // 34
        "CLCK", // 35
        "Egy0", // 36
        "GRAD", // 37
        "BHMA", // 38
        "BHMD", // 39
        "BHPC", // 40
        "ACRB",  // 41
        //New blocks
        "GLEN",  // 42
        "GOFF", // 43
        "MTOT", // 44
        "GPOS", // 45
        "MVIR", // 46
        "RVIR",  // 47
        "M25K", // 48
        "R25K",  // 49
        "M500", // 50
        "R500", // 51
        "MGAS", // 52
        "MSTR", // 53
        "TGAS", // 54
        "LGAS", // 55
        "NCON", // 56
        "MCON", // 57
        "BGPO", // 58
        "BGMA", // 59
        "BGRA", // 60
        "NSUB", // 61
        "FSUB", // 62
        "SLEN", // 63
        "SOFF", // 64
        "SSUB",  // 65
        "MSUB", // 66
        "SPOS",  // 67
        "SVEL", // 68
        "SCM", // 69
        "SPIN", // 70
        "DSUB", // 71
        "VMAX", // 72
        "RMAX", // 73
        "RHMS", // 74
        "MBID", // 75
        "GRNR", // 76
        "SMST", // 77
        "SLUM", // 78
        "SLAT", // 79
        "SLOB", // 80
        "DUST", // 81
        "SAGE", // 82
        "SZ",  // 83
        "SSFR" // 84
        //"PID"   // 85 long long?
    };
      
  const std::vector<std::vector<bool>> blocksFields = 
  { 
    {1,1,1,1,1,1}, // 0: POS, VEL, ID, MASS, IDU, TSTP, POT, ACCE
    {1,0,0,0,0,0}, // 1: U, TEMP, RHO, NE, NH, HSML, SFR, CLDX, ENDT, HOTT, MHOT, MCLD, EHOT, MSF, MFST, NMF, EOUT, EREC, EOLD, TDYN, SFRo, CLCK, Egy0, GRAD, GLEN, GOFF, MTOT, GPOS, MVIR, RVIR, M25K, R25K, M500, R500, MGAS, MSTR, TGAS, LGAS, NCON, MCON, NSUB, FSUB
    {0,0,0,0,1,1}, // 2: AGE
    {1,0,0,0,1,0}, // 3: Z, Zs, ZAGE, ZALV
    {0,0,0,0,1,0}, // 4: iM
    {0,0,0,0,0,1},  // 5: BHMA, BHMD, BHPC, ACRB
    {0,0,0,1,0,0},  // 6: BGPO, BGMA, BGRA
    {0,1,0,0,0,0}  // 7: SLEN, SOFF, SSUB, MSUB, SPOS, SVEL, SCM, SPIN, DSUB, VMAX, RMAX, RHMS, MBID, GRNR, SMST, SLUM, SLAT, SLOB, DUST, SAGE, SZ, SSFR
  };

  const std::vector<int> blockNamesToFields 
  {
    0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 2, 3, 3, 4, 3, 3, 1, 0, 0, 0,
    1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 5, 5, 5, 5,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 6, 6, 6, 1, 1,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7
  };

  const std::vector<int> blockSize 
  {
    3, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1,
    1, 1, 1, 3, 1, 1, 1, 1, 1, 1, 3, 3, 3, 3, 1, 1, 3, 1, 1, 1, 1,
    1, 1, 1, 1, 3, 3, 3, 3, 1, 1, 1, 1, 1, 1, 6, 12, 12, 12, 11, 1, 1, 1, 1
  };
  
};
  

#endif

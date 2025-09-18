
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

#ifndef COMMANDLINE_H
#define COMMANDLINE_H


#include <string>
#include <vector>
#include <abstractsource.h>
  

class CommandLine
{
  AbstractSource* pSource;
  public:
    CommandLine ( );
    ~CommandLine ( );
    int loadFile ();
    int parseOption (const std::vector<std::string> arguments );
    void showHelp();
    std::string getRemoteFile() {return m_remoteFile;} ; 
    std::string getType() {return m_type;} ; 

    std::vector<VSTable*>& getTable() ; 

    std::string   m_currentPath, m_type,m_ext, m_binaryDir, m_binaryPath,  m_binaryName,m_binaryHeader, m_file, m_endian, m_dataType, m_out, m_login, m_remoteFile,m_datasetList, m_aliasHeader, m_aliasParticle;
    bool m_historyEnabled;
    float m_missing, m_text;
    double m_size[3],m_comput[3]; //!volume data
    unsigned long long int m_npoints;
    std::vector<std::string> m_hyperslab;
    std::vector<std::string> m_fields;
    bool m_gLiteOut, m_inMemory = false;
    std::string m_lfn, m_VO, m_outlfn, m_se, m_outPath,m_historyFile;
    int m_fitshdunum;
 
/*!
m_currentPath = name including path of input file
m_binaryDir = output directory fof binary internal data format of VisIVO
m_binaryPath= filename including path of binary internal data format of VisIVO 
m_binaryName= only name of file (exluding path) of binary internal data format of VisIVO 
*/  
};


#endif


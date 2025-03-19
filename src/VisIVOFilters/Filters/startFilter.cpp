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
#include <cstdlib>
#include <cstring>

#include <iostream>
#include <sstream>
#include <string>

#include "vstable.h"

#include "vsselectcolumnsop.h"
#include "vsrandomizertableop.h"
#include "parametersparser.h"
#include "vsmergeop.h"
#include "vsappend.h"
#include "vsselectfieldop.h"
#include "vsmathop.h"
#include "vsdecimatorop.h"
#include "vsvisualop.h"
#include "vssphereop.h"
#include "vsstatisticop.h"
#include "vsshowtableop.h"
#include "vspointdistributeop.h"
#include "vscoarsevolumeop.h"
#include "vsextractsubvolumeop.h"
#include "vsexampleop.h"
#include "vspointpropertyop.h"
#include "vsinterpolateop.h"
#include "vsmoduleop.h"
#include "vssigmacontoursop.h"
#include "vspolarop.h"
#include "vsgrid2pointdistr.h"
#include "vsswapop.h"
#include "vsextractlistrowsop.h"
#include "vsaddidop.h"
#include "vschangecolnameop.h"
#include "vssplittableop.h"
#include "vswrvotableop.h"
#include "vscutop.h"
#include "vsvollimit.h"
#include "vsincludeop.h"
#include "visivoutils.h"
//∫#include "glitegw.h"
#include "startFilter.h"
#include "vsmrcampos.h"
//#include "vsclustercellop.h"
#include "vsmuportalop.h"

#ifdef AHF
#include "vsahfhalolistop.h"
#include "AHFstep/main_AHFstep.h"
#include "vsvbt2ahf.h"
#include "vsahfhalogalaxyextop.h"
#endif

#ifdef VSMPI
#include "mpi.h"
startFilter::startFilter(std::map<std::string,std::string> appParameters,MPI_Comm newcomm)
#else
startFilter::startFilter(std::map<std::string,std::string> appParameters)
#endif

{
    std::map<std::string,std::string>::iterator iter;
    int rank=0, size=1, ierr;
#ifdef VSMPI
    m_VS_COMM=newcomm;
    MPI_Errhandler_set(MPI_COMM_WORLD,MPI_ERRORS_RETURN);
    ierr=MPI_Comm_size (m_VS_COMM, &size);
    if(ierr!=MPI_SUCCESS) return;
    ierr=MPI_Comm_rank (m_VS_COMM, &rank);
    if(ierr!=MPI_SUCCESS) return;
    
    std::clog<<rank<<" TEST "<<size<<std::endl;
#endif
    
    iter =appParameters.find("history");
    {
        m_historyEnabled=true;
    }
    
    iter =appParameters.find("op");
    if( iter == appParameters.end())
    {
        if (rank==0) {
            iter =appParameters.find("help");
            if( iter == appParameters.end())
                std::cerr <<"No operation is requested"<<std::endl;
            std::cerr<<"VisIVOFilters version  3.5.0 March 19th 2025"<<std::endl<<std::endl;
            std::cerr <<"Syntax1: VisIVOFilters --op operation  [PARAMETERS] [--help]"<<std::endl;
            std::cerr <<"Syntax2: VisIVOFilters parameterFile"<<std::endl;
            std::cerr <<"valid operations: randomizer selcolumns merge append selfield mathop decimator extraction visualop showtable statistic pointdistribute pointproperty coarsevolume extractsubvolume interpolate module sigmacontours cartesian2polar AHFstep vbt2ahf grid2point extractlist addId ahfhalolist ahfhalogalaxyext changecolname splittable wrvotable include mres poca"<<std::endl;
        }
        return;
    }
    
    m_opName=iter->second.c_str();

    std::stringstream sstreamOp(iter->second);
    appParameters.erase(iter);
    int idOp=-1;
    

    
    if(sstreamOp.str()=="randomizer") idOp=1;
    if(sstreamOp.str()=="selcolumns") idOp=2;
    if(sstreamOp.str()=="merge") idOp=3;
    if(sstreamOp.str()=="append") idOp=4;
    if(sstreamOp.str()=="selfield") idOp=5;
    if(sstreamOp.str()=="mathop") idOp=6;
    if(sstreamOp.str()=="decimator") idOp=7;
    if(sstreamOp.str()=="visualop") idOp=8;
    if(sstreamOp.str()=="extraction") idOp=9;
    if(sstreamOp.str()=="showtable") idOp=10;
    if(sstreamOp.str()=="statistic") idOp=11;
    if(sstreamOp.str()=="pointdistribute") idOp=12;
    if(sstreamOp.str()=="coarsevolume") idOp=13;
    if(sstreamOp.str()=="extractsubvolume") idOp=14;
    if(sstreamOp.str()=="example") idOp=15;
    if(sstreamOp.str()=="pointproperty") idOp=16;
    if(sstreamOp.str()=="interpolate") idOp=17;
    if(sstreamOp.str()=="module") idOp=18;
    if(sstreamOp.str()=="sigmacontours") idOp=19;
    if(sstreamOp.str()=="cartesian2polar") idOp=20;//QUI
    if(sstreamOp.str()=="AHFstep") idOp=21;
    if(sstreamOp.str()=="vbt2ahf") idOp=22;
    if(sstreamOp.str()=="grid2point") idOp=23;
    if(sstreamOp.str()=="swap") idOp=24;
    if(sstreamOp.str()=="extractlist") idOp=25;
    if(sstreamOp.str()=="addId") idOp=26;
    if(sstreamOp.str()=="ahfhalolist") idOp=27;
    if(sstreamOp.str()=="changecolname") idOp=28;
    if(sstreamOp.str()=="ahfhalogalaxyext") idOp=29;
    if(sstreamOp.str()=="splittable") idOp=30;
    if(sstreamOp.str()=="wrvotable") idOp=31;
    if(sstreamOp.str()=="cut") idOp=32;
    if(sstreamOp.str()=="showvol") idOp=33;
    if(sstreamOp.str()=="include") idOp=34;
    if(sstreamOp.str()=="mres") idOp=35;
    //if(sstreamOp.str()=="clustercell") idOp=36;
    if(sstreamOp.str()=="poca") idOp=37;
    
#ifdef GLITE
    if(idOp==21 || idOp==22 ||idOp==27 ||idOp==29)
    {
        std::cerr <<"Filter not available for gLite version"<<std::endl;
        return;
    }
    bool fileIsLocal=true;
    gLiteGw gLInterface;
    std::string file1="";
    std::string file2="";
    std::string volFile=""; //grid2point --volume
    std::string localVolFile=""; //grid2point --volume
    std::string inputlfn;
    iter=appParameters.find("file");
    if( iter != appParameters.end() && iter->second.substr(0,6) =="lfn://")
    {
        fileIsLocal=false;
        inputlfn=iter->second;  //used for --append
    }
    iter=appParameters.find("infiles");
    if( iter != appParameters.end())
    {
        std::stringstream ssInfileparameters;
        ssInfileparameters.str(iter->second);
        ssInfileparameters>>file1;
        ssInfileparameters>>file2;
        if(file1.substr(0,6) !="lfn://") file1="";
        if(file2.substr(0,6) !="lfn://") file2="";
    }
    std::string localFilename=gLInterface.readVF(appParameters,idOp);
    if(localFilename=="NogLiteSuccess") return;
    
    if(idOp==23) //grid2point --volume
    {
        iter=appParameters.find("volume");
        if( iter != appParameters.end())
        {
            volFile=iter->second;
            
            if(volFile.substr(0,6) =="lfn://")
            {
                if(volFile.find(".bin") == std::string::npos)
                    volFile.append(".bin");
                localVolFile=getName(volFile);
                
                iter=appParameters.find("out");
                if(iter != appParameters.end())
                {
                    std::string dir=getDir(iter->second);
                    if(dir==""){
                        localVolFile="/"+localVolFile;
                        localVolFile=getenv("PWD")+localVolFile;
                    }
                    else
                        localVolFile=dir+localVolFile;
                } else {
                    localVolFile="/"+localVolFile;
                    localVolFile=getenv("PWD")+localVolFile;
                }
                
                iter=appParameters.find("VO");
                std::string vo=iter->second;
                if(!remoteDownloadFiles(volFile,vo,localVolFile))
                    return;
                
                std::string head1=volFile+".head";
                std::string head2=localVolFile+".head";
                
                if(!remoteDownloadFiles(head1,vo,head2))
                    return;
                
                iter =appParameters.find("volume");
                appParameters.erase(iter);
                appParameters.insert(make_pair("volume",localVolFile));
            }
        }
        
    }
    
    if(localFilename!="")
    {
        if(idOp==3 || idOp==4 || idOp==8)
        {
            iter =appParameters.find("filelist");
            appParameters.erase(iter);
            appParameters.insert(make_pair("filelist",localFilename));
            
        }else if(idOp==17) {
            iter =appParameters.find("infiles");
            appParameters.erase(iter);
            appParameters.insert(make_pair("infiles",localFilename));
        }else if(idOp==-1) {
            //nothing todo
        }else{
            iter =appParameters.find("file");
            appParameters.erase(iter);
            appParameters.insert(make_pair("file",localFilename));
        }
    }
#endif
    std::vector <std::string> valOutFilename;
    std::string filename;
    switch(idOp)
    {
            /*** Randomizer OP **/
        case 1:
        {
            //serial
            if(rank==0)
            {
                iter =appParameters.find("help");
                if( iter != appParameters.end())
                {
                    VSRandomizerTableOp op;
                    op.printHelp();
                    return;
                }
                
                iter =appParameters.find("file");
                if( iter == appParameters.end())
                {
                    if(rank==0) std::cerr <<"No input file table is provided"<<std::endl;
                    return;
                }
                std::stringstream sFilename(iter->second);
                
                sFilename>>filename;
                if(filename.find(".bin") == std::string::npos)
                    filename.append(".bin");
                VSTable table(filename);
                
                if(!table.tableExist())
                {
                    std::cerr <<"No valid input file table is provided"<<std::endl;
                    return;
                }
                
                VSRandomizerTableOp op;
                op.setParameters(appParameters);
                op.addInput(&table);
                op.execute();
                valOutFilename=op.realOutFilename();
                
                if(m_historyEnabled)
                {
                    const char* m_historyFile ="hist.xml";
                    
                    iter =appParameters.find("historyfile");
                    if( iter != appParameters.end())
                    {
                        m_historyFile = appParameters.find("historyfile")->second.c_str();
                    }
                    
                    op.writeHistory(m_historyFile,m_opName,appParameters,valOutFilename);
                }
                break;
            }
            
            
            
        }
            /***END Randomizer OP **/
            /*** Select Columns OP **/
        case 2:
        {
            //serial
            if(rank==0)  //serialized
            {
                iter =appParameters.find("help");
                if( iter != appParameters.end())
                {
                    VSSelectColumnsOp op;
                    op.printHelp();
                    return;
                }
                iter =appParameters.find("file");
                if( iter == appParameters.end())
                {
                    std::cerr <<"No input file table is provided"<<std::endl;
                    return;
                }
                std::stringstream sFilename(iter->second);
                //appParameters.erase(iter);
                sFilename>>filename;
                if(filename.find(".bin") == std::string::npos)
                    filename.append(".bin");
                VSTable table(filename);
                if(!table.tableExist())
                {
                    std::cerr <<"No valid input file table is provided"<<std::endl;
                    return;
                }
                
                
                VSSelectColumnsOp op;
                op.setParameters(appParameters);
                op.addInput(&table);
                op.execute();
                valOutFilename=op.realOutFilename();
                
                if(m_historyEnabled)
                {
                    const char* m_historyFile ="hist.xml";
                    
                    iter =appParameters.find("historyfile");
                    if( iter != appParameters.end())
                    {
                        m_historyFile = appParameters.find("historyfile")->second.c_str();
                    }
                    
                    op.writeHistory(m_historyFile,m_opName,appParameters,valOutFilename);
                }

            }
            break;
        }
            /*** END Select Columns OP **/
            /*** Merge Tables OP **/
        case 3:
        {
            //serial
            if(rank==0)  //serialized
            {
                
                iter =appParameters.find("help");
                if( iter != appParameters.end())
                {
                    VSMergeOp op;
                    op.printHelp();
                    return;
                }
                VSMergeOp op;
                op.setParameters(appParameters);
                op.execute();
                valOutFilename=op.realOutFilename();
                
                if(m_historyEnabled)
                {
                    const char* m_historyFile ="hist.xml";
                    
                    iter =appParameters.find("historyfile");
                    if( iter != appParameters.end())
                    {
                        m_historyFile = appParameters.find("historyfile")->second.c_str();
                    }
                    
                    op.writeHistory(m_historyFile,m_opName,appParameters,valOutFilename);
                }
            }
            break;
        }
            /*** END Merge Tables OP **/
            /*** Append Tables OP **/
        case 4:
        {
            //serial
            if(rank==0)  //serialized
            {
                iter =appParameters.find("help");
                if( iter != appParameters.end())
                {
                    VSAppendOp op;
                    op.printHelp();
                    return;
                }
                VSAppendOp op;
                op.setParameters(appParameters);
                op.execute();
                valOutFilename=op.realOutFilename();
                
                if(m_historyEnabled)
                {
                    const char* m_historyFile ="hist.xml";
                    
                    iter =appParameters.find("historyfile");
                    if( iter != appParameters.end())
                    {
                        m_historyFile = appParameters.find("historyfile")->second.c_str();
                    }
                    
                    op.writeHistory(m_historyFile,m_opName,appParameters,valOutFilename);
                }

            }
            break;
        }
            /*** END Append Tables OP **/
            /*** SelectField  OP **/
        case 5:
        {
            //serial
            if(rank==0)  //serialized
            {
                
                iter =appParameters.find("help");
                if( iter != appParameters.end())
                {
                    VSSelectFieldOp op;
                    op.printHelp();
                    return;
                }
                
                iter =appParameters.find("file");
                if( iter == appParameters.end())
                {
                    std::cerr <<"No input file table is provided"<<std::endl;
                    return;
                }
                std::stringstream sFilename(iter->second);
                //appParameters.erase(iter);
                sFilename>>filename;
                if(filename.find(".bin") == std::string::npos)
                    filename.append(".bin");
                VSTable table(filename);
                if(!table.tableExist())
                {
                    std::cerr <<"No valid input file table is provided"<<std::endl;
                    return;
                }
                
                VSSelectFieldOp op;
                op.setParameters(appParameters);
                op.addInput(&table);
                op.execute();
                valOutFilename=op.realOutFilename();
                if(m_historyEnabled)
                {
                    const char* m_historyFile ="hist.xml";
                    
                    iter =appParameters.find("historyfile");
                    if( iter != appParameters.end())
                    {
                        m_historyFile = appParameters.find("historyfile")->second.c_str();
                    }
                    
                    op.writeHistory(m_historyFile,m_opName,appParameters,valOutFilename);
                }

            }
            break;
        }
            /*** END SelectField  OP **/
            /*** Math  OP **/
        case 6:
        {
            //serial
            if(rank==0)  //serialized
            {
                
                iter =appParameters.find("help");
                if( iter != appParameters.end())
                {
                    VSMathOp op;
                    op.printHelp();
                    return;
                }
                
                iter =appParameters.find("file");
                if( iter == appParameters.end())
                {
                    std::cerr <<"No input file table is provided"<<std::endl;
                    return;
                }
                std::stringstream sFilename(iter->second);
                //appParameters.erase(iter);
                sFilename>>filename;
                if(filename.find(".bin") == std::string::npos)
                    filename.append(".bin");
                VSTable table(filename);
                if(!table.tableExist())
                {
                    std::cerr <<"No valid input file table is provided"<<std::endl;
                    return;
                }
                
                VSMathOp op;
                op.setParameters(appParameters);
                op.addInput(&table);
                op.execute();
                valOutFilename=op.realOutFilename();
                
                if(m_historyEnabled)
                {
                    const char* m_historyFile ="hist.xml";
                    
                    iter =appParameters.find("historyfile");
                    if( iter != appParameters.end())
                    {
                        m_historyFile = appParameters.find("historyfile")->second.c_str();
                    }
                    
                    op.writeHistory(m_historyFile,m_opName,appParameters,valOutFilename);
                }

            }
            break;
        }
            /*** END Math  OP **/
            /*** Decimator OP **/
        case 7:
        {
            //serial
            if(rank==0)  //serialized
            {
                
                iter =appParameters.find("help");
                if( iter != appParameters.end())
                {
                    VSDecimatorTableOp op;
                    op.printHelp();
                    return;
                }
                
                iter =appParameters.find("file");
                if( iter == appParameters.end())
                {
                    std::cerr <<"No input file table is provided"<<std::endl;
                    return;
                }
                
                
                std::stringstream sFilename(iter->second);
                //appParameters.erase(iter);
                
                sFilename>>filename;
                if(filename.find(".bin") == std::string::npos)
                    filename.append(".bin");
                VSTable table(filename);
                if(!table.tableExist())
                {
                    std::cerr <<"No valid input file table is provided"<<std::endl;
                    return;
                }
                
                VSDecimatorTableOp op;
                op.setParameters(appParameters);
                op.addInput(&table);
                op.execute();
                valOutFilename=op.realOutFilename();
                
                if(m_historyEnabled)
                {
                    const char* m_historyFile ="hist.xml";
                    
                    iter =appParameters.find("historyfile");
                    if( iter != appParameters.end())
                    {
                        m_historyFile = appParameters.find("historyfile")->second.c_str();
                    }
                    
                    op.writeHistory(m_historyFile,m_opName,appParameters,valOutFilename);
                }

            }
            break;
        }
            /***END Randomizer OP **/
            /*** Visual OP **/
        case 8:
        {
            //serial
            if(rank==0)  //serialized
            {
                
                iter =appParameters.find("help");
                if( iter != appParameters.end())
                {
                    VSVisualOp op;
                    op.printHelp();
                    return;
                }
                VSVisualOp op;
                op.setParameters(appParameters);
                op.execute();
                valOutFilename=op.realOutFilename();
                
                if(m_historyEnabled)
                {
                    const char* m_historyFile ="hist.xml";
                    
                    iter =appParameters.find("historyfile");
                    if( iter != appParameters.end())
                    {
                        m_historyFile = appParameters.find("historyfile")->second.c_str();
                    }
                    
                    op.writeHistory(m_historyFile,m_opName,appParameters,valOutFilename);
                }

            }
            break;
        }
            /*** END Visual OP **/
            /*** Sphere  OP **/
        case 9:
        {
            //serial
            if(rank==0)  //serialized
            {
                
                iter =appParameters.find("help");
                if( iter != appParameters.end())
                {
                    VSSphereOp op;
                    op.printHelp();
                    return;
                }
                
                iter =appParameters.find("file");
                if( iter == appParameters.end())
                {
                    std::cerr <<"No input file table is provided"<<std::endl;
                    return;
                }
                std::stringstream sFilename(iter->second);
                //appParameters.erase(iter);
                sFilename>>filename;
                if(filename.find(".bin") == std::string::npos)
                    filename.append(".bin");
                VSTable table(filename);
                if(!table.tableExist())
                {
                    std::cerr <<"No valid input file table is provided"<<std::endl;
                    return;
                }
                
                VSSphereOp op;
                op.setParameters(appParameters);
                op.addInput(&table);
                op.execute();
                valOutFilename=op.realOutFilename();
                if(m_historyEnabled)
                {
                    const char* m_historyFile ="hist.xml";
                    
                    iter =appParameters.find("historyfile");
                    if( iter != appParameters.end())
                    {
                        m_historyFile = appParameters.find("historyfile")->second.c_str();
                    }
                    
                    op.writeHistory(m_historyFile,m_opName,appParameters,valOutFilename);
                }

            }
            break;
        }
            /*** END SelectField  OP **/
            /*** ShowTable  OP **/
        case 10:
        {
            //serial
            if(rank==0)  //serialized
            {
                
                iter =appParameters.find("help");
                if( iter != appParameters.end())
                {
                    VSShowTableOp op;
                    op.printHelp();
                    return;
                }
                
                iter =appParameters.find("file");
                
                
                if( iter == appParameters.end())
                {
                    std::cerr <<"No input file table is provided"<<std::endl;
                    return;
                }
                std::stringstream sFilename(iter->second);
                //appParameters.erase(iter);
                sFilename>>filename;
                if(filename.find(".bin") == std::string::npos)
                    filename.append(".bin");
                VSTable table(filename);
                if(!table.tableExist())
                {
                    std::cerr <<"No valid input file table is provided"<<std::endl;
                    return;
                }
                
                VSShowTableOp op;
                op.setParameters(appParameters);
                op.addInput(&table);
                op.execute();
                valOutFilename=op.realOutFilename();
                if(m_historyEnabled)
                {
                    const char* m_historyFile ="hist.xml";
                    
                    iter =appParameters.find("historyfile");
                    if( iter != appParameters.end())
                    {
                        m_historyFile = appParameters.find("historyfile")->second.c_str();
                    }
                    
                    op.writeHistory(m_historyFile,m_opName,appParameters,valOutFilename);
                }

            }
            break;
        }
            /*** END ShowTable  OP **/
            /*** Statistic  OP **/
        case 11:
        {
            //serial
            if(rank==0)  //serialized
            {
                
                iter =appParameters.find("help");
                if(iter != appParameters.end())
                {
                    VSStatisticOp op;
                    op.printHelp();
                    return;
                }
                
                iter = appParameters.find("file");
                

                if(iter == appParameters.end())
                {
                    std::cerr <<"No input file table is provided"<<std::endl;
                    return;
                }
                std::stringstream sFilename(iter->second);
               // appParameters.erase(iter);
                sFilename>>filename;
                if(filename.find(".bin") == std::string::npos)
                    filename.append(".bin");
                VSTable table(filename);
                if(!table.tableExist())
                {
                    std::cerr <<"No valid input file table is provided"<<std::endl;
                    return;
                }
                
                iter = appParameters.find("file");
                
                
                VSStatisticOp op;
                op.setParameters(appParameters);
                op.addInput(&table);
                op.execute();
                valOutFilename=op.realOutFilename();
                
                if(m_historyEnabled)
                {
                    const char* m_historyFile ="hist.xml";
                    
                    iter =appParameters.find("historyfile");
                    if( iter != appParameters.end())
                    {
                        m_historyFile = appParameters.find("historyfile")->second.c_str();
                    }
                    
                    op.writeHistory(m_historyFile,m_opName,appParameters,valOutFilename);
                }

            }
            break;
        }
            /*** END Statistic OP **/
            /*** PointDistribute  OP **/
        case 12:
        {
            //serial
            if(rank==0)  //serialized
            {
                
                iter =appParameters.find("help");
                if( iter != appParameters.end())
                {
                    VSPointDistributeOp op;
                    op.printHelp();
                    return;
                }
                
                iter =appParameters.find("file");
                if( iter == appParameters.end())
                {
                    std::cerr <<"No input file table is provided"<<std::endl;
                    return;
                }
                std::stringstream sFilename(iter->second);
               // appParameters.erase(iter);
                sFilename>>filename;
                if(filename.find(".bin") == std::string::npos)
                    filename.append(".bin");
                VSTable table(filename);
                if(!table.tableExist())
                {
                    std::cerr <<"No valid input file table is provided"<<std::endl;
                    return;
                }
                
                VSPointDistributeOp op;
                op.setParameters(appParameters);
                op.addInput(&table);
                op.execute();
                valOutFilename=op.realOutFilename();
                if(m_historyEnabled)
                {
                    const char* m_historyFile ="hist.xml";
                    
                    iter =appParameters.find("historyfile");
                    if( iter != appParameters.end())
                    {
                        m_historyFile = appParameters.find("historyfile")->second.c_str();
                    }
                    
                    op.writeHistory(m_historyFile,m_opName,appParameters,valOutFilename);
                }

            }
            break;
        }
            /*** END PointDistribute OP **/
            /*** CoarseVolume  OP **/
        case 13:
        {
            //serial
            if(rank==0)  //serialized
            {
                
                iter =appParameters.find("help");
                if( iter != appParameters.end())
                {
                    VSCoarseVolumeOp op;
                    op.printHelp();
                    return;
                }
                
                iter =appParameters.find("file");
                if( iter == appParameters.end())
                {
                    std::cerr <<"No input file table is provided"<<std::endl;
                    return;
                }
                std::stringstream sFilename(iter->second);
                //appParameters.erase(iter);
                sFilename>>filename;
                if(filename.find(".bin") == std::string::npos)
                    filename.append(".bin");
                VSTable table(filename);
                if(!table.tableExist())
                {
                    std::cerr <<"No valid input file table is provided"<<std::endl;
                    return;
                }
                
                VSCoarseVolumeOp op;
                op.setParameters(appParameters);
                op.addInput(&table);
                op.execute();
                valOutFilename=op.realOutFilename();
                if(m_historyEnabled)
                {
                    const char* m_historyFile ="hist.xml";
                    
                    iter =appParameters.find("historyfile");
                    if( iter != appParameters.end())
                    {
                        m_historyFile = appParameters.find("historyfile")->second.c_str();
                    }
                    
                    op.writeHistory(m_historyFile,m_opName,appParameters,valOutFilename);
                }

            }
            break;
        }
            /*** END CoarseVolume OP **/
            /*** ExtractSubVolume  OP **/
        case 14:
        {
            //serial
            if(rank==0)  //serialized
            {
                
                iter =appParameters.find("help");
                if( iter != appParameters.end())
                {
                    VSExtractSubVolumeOp op;
                    op.printHelp();
                    return;
                }
                
                iter =appParameters.find("file");
                if( iter == appParameters.end())
                {
                    std::cerr <<"No input file table is provided"<<std::endl;
                    return;
                }
                std::stringstream sFilename(iter->second);
                //appParameters.erase(iter);
                sFilename>>filename;
                if(filename.find(".bin") == std::string::npos)
                    filename.append(".bin");
                VSTable table(filename);
                if(!table.tableExist())
                {
                    std::cerr <<"No valid input file table is provided"<<std::endl;
                    return;
                }
                
                VSExtractSubVolumeOp op;
                op.setParameters(appParameters);
                op.addInput(&table);
                op.execute();
                valOutFilename=op.realOutFilename();
                if(m_historyEnabled)
                {
                    const char* m_historyFile ="hist.xml";
                    
                    iter =appParameters.find("historyfile");
                    if( iter != appParameters.end())
                    {
                        m_historyFile = appParameters.find("historyfile")->second.c_str();
                    }
                    
                    op.writeHistory(m_historyFile,m_opName,appParameters,valOutFilename);
                }
                
            }
            break;
        }
            /*** END ExtracSubVolume OP **/
            /*** Example  OP Used for developers only **/
        case 15:
        {
            //serial
            if(rank==0)  //serialized
            {
                
                iter =appParameters.find("help");  //! if --help is given, print help
                if( iter != appParameters.end())
                {
                    VSExampleOp op;
                    op.printHelp();
                    return;
                }
                
                iter =appParameters.find("file");  //!check for --file parameter
                if( iter == appParameters.end())
                {
                    std::cerr <<"No input file table is provided"<<std::endl;
                    return;
                }
                std::stringstream sFilename(iter->second);
                sFilename>>filename; //! name of the file binary data format
                if(filename.find(".bin") == std::string::npos)  //! table name must contain .bin extension
                    filename.append(".bin");
                VSTable table(filename); //! create an object VSTable
                if(!table.tableExist()) //! check that the file exist
                {
                    std::cerr <<"No valid input file table is provided"<<std::endl;
                    return;
                }
                
                VSExampleOp op; //! specific operation to be executed (defined in a cpp class. All classes must be derived from VSTableOp)
                
                op.setParameters(appParameters); //!parser of input line command: fill a map(key=string option, value=string value) defined in VSTableOP. The map is a private member  of the VSTableOP. Contained values can be obtained with the method getParameterAsString, getParameterAsFloat, getParameterAsInt public functions of the class. Parameters can be also set with addParameter method.
                
                op.addInput(&table); //! add the created table to the operation (defined in VSTableOp)    			std::vector<VSTable *> m_tables; m_tables[i] is a pointer to a 			VSTable object.
                
                op.execute(); //! execute the action specific of the operation
                valOutFilename=op.realOutFilename();
                
            }
            break;
        }
            /*** END Example OP **/
            /*** PointProperty  OP **/
        case 16:
        {
            //serial
            if(rank==0)  //serialized
            {
                
                iter =appParameters.find("help");
                if( iter != appParameters.end())
                {
                    VSPointPropertyOp op;
                    op.printHelp();
                    return;
                }
                
                iter =appParameters.find("file");
                if( iter == appParameters.end())
                {
                    std::cerr <<"No input file table is provided"<<std::endl;
                    return;
                }
                std::stringstream sFilename(iter->second);
                //appParameters.erase(iter);
                sFilename>>filename;
                if(filename.find(".bin") == std::string::npos)
                    filename.append(".bin");
                VSTable table(filename);
                if(!table.tableExist())
                {
                    std::cerr <<"No valid input file table is provided"<<std::endl;
                    return;
                }
                
                VSPointPropertyOp op;
                op.setParameters(appParameters);
                op.addInput(&table);
                op.execute();
                valOutFilename=op.realOutFilename();
                if(m_historyEnabled)
                {
                    const char* m_historyFile ="hist.xml";
                    
                    iter =appParameters.find("historyfile");
                    if( iter != appParameters.end())
                    {
                        m_historyFile = appParameters.find("historyfile")->second.c_str();
                    }
                    
                    op.writeHistory(m_historyFile,m_opName,appParameters,valOutFilename);
                }

            }
            break;
        }
            /*** END PointProperty OP **/
            /*** Interpolate OP **/
        case 17:
        {
            //serial
            if(rank==0)  //serialized
            {
                
                iter =appParameters.find("help");
                if( iter != appParameters.end())
                {
                    VSInterpolateOp op;
                    op.printHelp();
                    return;
                }
                VSInterpolateOp op;
                op.setParameters(appParameters);
                op.execute();
                valOutFilename=op.realOutFilename();
                if(m_historyEnabled)
                {
                    const char* m_historyFile ="hist.xml";
                    
                    iter =appParameters.find("historyfile");
                    if( iter != appParameters.end())
                    {
                        m_historyFile = appParameters.find("historyfile")->second.c_str();
                    }
                    
                    op.writeHistory(m_historyFile,m_opName,appParameters,valOutFilename);
                }

                
                break;
            }
            /*** END Merge Tables OP **/
        case 18:
            {
                iter =appParameters.find("help");
                if( iter != appParameters.end())
                {
                    VSModuleOp op;
                    op.printHelp();
                    return;
                }
                iter =appParameters.find("file");
                if( iter == appParameters.end())
                {
                    std::cerr <<"No input file table is provided"<<std::endl;
                    return;
                }
                std::stringstream sFilename(iter->second);
                //appParameters.erase(iter);
                sFilename>>filename;
                if(filename.find(".bin") == std::string::npos)
                    filename.append(".bin");
                VSTable table(filename);
                if(!table.tableExist())
                {
                    std::cerr <<"No valid input file table is provided"<<std::endl;
                    return;
                }
                VSModuleOp op;
                op.setParameters(appParameters);
                op.addInput(&table);
                op.execute();
                valOutFilename=op.realOutFilename();
                if(m_historyEnabled)
                {
                    const char* m_historyFile ="hist.xml";
                    
                    iter =appParameters.find("historyfile");
                    if( iter != appParameters.end())
                    {
                        m_historyFile = appParameters.find("historyfile")->second.c_str();
                    }
                    
                    op.writeHistory(m_historyFile,m_opName,appParameters,valOutFilename);
                }

            }
            break;
        }
        case 19:
        {
            //serial
            if(rank==0)  //serialized
            {
                
                iter =appParameters.find("help");
                if( iter != appParameters.end()){
                    VSSigmaContoursOp op;
                    op.printHelp();
                    return;
                }
                
                
                VSSigmaContoursOp op;
                op.setParameters(appParameters);
                iter= appParameters.find("file");
                if(iter==appParameters.end()){
                    std::cerr<<"No input file table is provided!"<<std::endl;
                    return;
                }//close if
                std::stringstream sFilename(iter->second);
                //appParameters.erase(iter);
                sFilename>>filename;
                if(filename.find(".bin") == std::string::npos) filename.append(".bin");
                VSTable table(filename);
                if(!table.tableExist()){
                    std::cerr<<"No valid input file table is provided"<<std::endl;
                    return;
                }
                
                op.addInput(&table);
                op.execute();
                valOutFilename=op.realOutFilename();
                if(m_historyEnabled)
                {
                    const char* m_historyFile ="hist.xml";
                    
                    iter =appParameters.find("historyfile");
                    if( iter != appParameters.end())
                    {
                        m_historyFile = appParameters.find("historyfile")->second.c_str();
                    }
                    
                    op.writeHistory(m_historyFile,m_opName,appParameters,valOutFilename);
                }

            }
            break;
        }
            /*** END SIGMA CONTOURS OP **/
            /*** END SIGMA CONTOURS OP **/
            /*** END SIGMA CONTOURS OP **/
        case 20:
        {
            //serial
            if(rank==0)  //serialized
            {
                iter =appParameters.find("help");
                if( iter != appParameters.end()){
                    VSPolarOp op;
                    op.printHelp();
                    return;
                }
                
                
                VSPolarOp op;
                op.setParameters(appParameters);
                iter= appParameters.find("file");
                if(iter==appParameters.end()){
                    std::cerr<<"No input file table is provided!"<<std::endl;
                    return;
                }//close if
                std::stringstream sFilename(iter->second);
                //appParameters.erase(iter);
                sFilename>>filename;
                if(filename.find(".bin") == std::string::npos) filename.append(".bin");
                VSTable table(filename);
                if(!table.tableExist()){
                    std::cerr<<"No valid input file table is provided"<<std::endl;
                    return;
                }
                
                op.addInput(&table);
                op.execute();
                valOutFilename=op.realOutFilename();
                if(m_historyEnabled)
                {
                    const char* m_historyFile ="hist.xml";
                    
                    iter =appParameters.find("historyfile");
                    if( iter != appParameters.end())
                    {
                        m_historyFile = appParameters.find("historyfile")->second.c_str();
                    }
                    
                    op.writeHistory(m_historyFile,m_opName,appParameters,valOutFilename);
                }

            }
            break;
        }
        case 21:
        {
#ifdef AHF
#ifdef WIN32
            if(rank==0) std::cerr<<"Operation not allowed on this platform";
            return;
#else
            //serial
            if(rank==0)  //serialized
            {
                
                iter =appParameters.find("help");
                if( iter != appParameters.end())
                {
                    std::cout<<"VisIVOFilters --op AHFstep --input AmigaParameterFile"<<std::endl;
                    return;
                }
                iter= appParameters.find("input");
                if(iter==appParameters.end()){
                    std::cerr<<"No input Amiga Parameter File is provided!"<<std::endl;
                    return;
                }//close if
                std::stringstream sFilename(iter->second);
                
                std::string fileInput;
                sFilename>>fileInput;
                char *parFile=new char[fileInput.size()+1];
                parFile[fileInput.size()]=0;
                memcpy(parFile,fileInput.c_str(),fileInput.size());
                main_AHFstep(parFile,&argc,&argv);
                delete [] parFile;
            }
#endif
#else
            if(rank==0) std::cerr<<"Operation not allowed. Compile VisIVO with AHF option";
            return;
#endif
            break;
        }
        case 22:
        {
#ifdef AHF
            //serial
            if(rank==0)  //serialized
            {
                
                iter =appParameters.find("help");
                if( iter != appParameters.end()){
                    VSvbt2ahf op;
                    op.printHelp();
                    return;
                }
                
                
                VSvbt2ahf op;
                op.setParameters(appParameters);
                iter= appParameters.find("file");
                if(iter==appParameters.end()){
                    std::cerr<<"No input file table is provided!"<<std::endl;
                    return;
                }//close if
                std::stringstream sFilename(iter->second);
                //appParameters.erase(iter);
                sFilename>>filename;
                if(filename.find(".bin") == std::string::npos) filename.append(".bin");
                VSTable table(filename);
                if(!table.tableExist()){
                    std::cerr<<"No valid input file table is provided"<<std::endl;
                    return;
                }
                
                op.addInput(&table);
                op.execute();
                valOutFilename=op.realOutFilename();
                if(m_historyEnabled)
                {
                    const char* m_historyFile ="hist.xml";
                    
                    iter =appParameters.find("historyfile");
                    if( iter != appParameters.end())
                    {
                        m_historyFile = appParameters.find("historyfile")->second.c_str();
                    }
                    
                    op.writeHistory(m_historyFile,m_opName,appParameters,valOutFilename);
                }

            }
#else
            if(rank==0) std::cerr<<"Operation not allowed. Compile VisIVO with AHF option";
            return;
#endif
            
            break;
        }
        case 23:
        {
            //serial
            if(rank==0)  //serialized
            {
                
                iter =appParameters.find("help");
                
                if( iter != appParameters.end()){
                    VSGrid2PointDistr op;
                    op.printHelp();
                    return;
                }
                
                VSGrid2PointDistr op;
                op.setParameters(appParameters);
                iter= appParameters.find("file");
                
                if(iter==appParameters.end()){
                    std::cerr<<"No input file table is provided!"<<std::endl;
                    return;
                }//close if
                
                std::stringstream sFilename(iter->second);
                //appParameters.erase(iter);
                sFilename>>filename;
                if(filename.find(".bin") == std::string::npos) filename.append(".bin");
                VSTable table(filename);
                
                if(!table.tableExist()){
                    std::cerr<<"No valid input file table is provided"<<std::endl;
                    return;
                }
                
                op.addInput(&table);
                op.execute();
                valOutFilename=op.realOutFilename();
                if(m_historyEnabled)
                {
                    const char* m_historyFile ="hist.xml";
                    
                    iter =appParameters.find("historyfile");
                    if( iter != appParameters.end())
                    {
                        m_historyFile = appParameters.find("historyfile")->second.c_str();
                    }
                    
                    op.writeHistory(m_historyFile,m_opName,appParameters,valOutFilename);
                }

            }
            break;
        }
            /******* Swap */
        case 24:
        {
            //serial
            if(rank==0)  //serialized
            {
                
                iter =appParameters.find("help");
                if( iter != appParameters.end()){
                    VSSwapOp op;
                    op.printHelp();
                    return;
                }
                
                
                VSSwapOp op;
                op.setParameters(appParameters);
                iter= appParameters.find("file");
                if(iter==appParameters.end()){
                    std::cerr<<"No input file table is provided!"<<std::endl;
                    return;
                }//close if
                std::stringstream sFilename(iter->second);
                //appParameters.erase(iter);
                sFilename>>filename;
                if(filename.find(".bin") == std::string::npos) filename.append(".bin");
                VSTable table(filename);
                if(!table.tableExist()){
                    std::cerr<<"No valid input file table is provided"<<std::endl;
                    return;
                }
                
                op.addInput(&table);
                op.execute();
                valOutFilename=op.realOutFilename();
                if(m_historyEnabled)
                {
                    const char* m_historyFile ="hist.xml";
                    
                    iter =appParameters.find("historyfile");
                    if( iter != appParameters.end())
                    {
                        m_historyFile = appParameters.find("historyfile")->second.c_str();
                    }
                    
                    op.writeHistory(m_historyFile,m_opName,appParameters,valOutFilename);
                }

            }
            break;
        }
            /******* ExtractList */
        case 25:
        {
            //serial
            if(rank==0)  //serialized
            {
                
                iter =appParameters.find("help");
                if( iter != appParameters.end()){
                    VSExtractListRowsOp op;
                    op.printHelp();
                    return;
                }
                
                
                VSExtractListRowsOp op;
                op.setParameters(appParameters);
                iter= appParameters.find("file");
                if(iter==appParameters.end()){
                    std::cerr<<"No input file table is provided!"<<std::endl;
                    return;
                }//close if
                std::stringstream sFilename(iter->second);
                //appParameters.erase(iter);
                sFilename>>filename;
                if(filename.find(".bin") == std::string::npos) filename.append(".bin");
                VSTable table(filename);
                if(!table.tableExist()){
                    std::cerr<<"No valid input file table is provided"<<std::endl;
                    return;
                }
                
                op.addInput(&table);
                op.execute();
                valOutFilename=op.realOutFilename();
                if(m_historyEnabled)
                {
                    const char* m_historyFile ="hist.xml";
                    
                    iter =appParameters.find("historyfile");
                    if( iter != appParameters.end())
                    {
                        m_historyFile = appParameters.find("historyfile")->second.c_str();
                    }
                    
                    op.writeHistory(m_historyFile,m_opName,appParameters,valOutFilename);
                }

            }
            break;
        }
            /******* AddId */
        case 26:
        {
            //serial
            if(rank==0)  //serialized
            {
                
                iter =appParameters.find("help");
                if( iter != appParameters.end()){
                    VSAddIdOp op;
                    op.printHelp();
                    return;
                }
                
                
                VSAddIdOp op;
                op.setParameters(appParameters);
                iter= appParameters.find("file");
                if(iter==appParameters.end()){
                    std::cerr<<"No input file table is provided!"<<std::endl;
                    return;
                }//close if
                std::stringstream sFilename(iter->second);
                //appParameters.erase(iter);
                sFilename>>filename;
                if(filename.find(".bin") == std::string::npos) filename.append(".bin");
                VSTable table(filename);
                if(!table.tableExist()){
                    std::cerr<<"No valid input file table is provided"<<std::endl;
                    return;
                }
                
                op.addInput(&table);
                op.execute();
                valOutFilename=op.realOutFilename();
                if(m_historyEnabled)
                {
                    const char* m_historyFile ="hist.xml";
                    
                    iter =appParameters.find("historyfile");
                    if( iter != appParameters.end())
                    {
                        m_historyFile = appParameters.find("historyfile")->second.c_str();
                    }
                    
                    op.writeHistory(m_historyFile,m_opName,appParameters,valOutFilename);
                }

            }
            break;
        }
            /******* AhfHaloId */
        case 27:
        {
#ifdef AHF
            //serial
            if(rank==0)  //serialized
            {
                
                iter =appParameters.find("help");
                if( iter != appParameters.end()){
                    VSAhfHaloListOp op;
                    op.printHelp();
                    return;
                }
                
                
                VSAhfHaloListOp op;
                op.setParameters(appParameters);
                iter= appParameters.find("file");
                if(iter==appParameters.end()){
                    std::cerr<<"No input file table is provided!"<<std::endl;
                    return;
                }//close if
                std::stringstream sFilename(iter->second);
                //appParameters.erase(iter);
                sFilename>>filename;
                if(filename.find(".bin") == std::string::npos) filename.append(".bin");
                VSTable table(filename);
                if(!table.tableExist()){
                    std::cerr<<"No valid input file table is provided"<<std::endl;
                    return;
                }
                
                op.addInput(&table);
                op.execute();
                valOutFilename=op.realOutFilename();
                if(m_historyEnabled)
                {
                    const char* m_historyFile ="hist.xml";
                    
                    iter =appParameters.find("historyfile");
                    if( iter != appParameters.end())
                    {
                        m_historyFile = appParameters.find("historyfile")->second.c_str();
                    }
                    
                    op.writeHistory(m_historyFile,m_opName,appParameters,valOutFilename);
                }

            }
#else
            if(rank==0) std::cerr<<"Operation not allowed. Compile VisIVO with AHF option";
            return;
#endif
            
            break;
        }
            /******* ChangeColName*/
        case 28:
        {
            //serial
            if(rank==0)  //serialized
            {
                
                iter =appParameters.find("help");
                if( iter != appParameters.end()){
                    VSChangeColNameop op;
                    op.printHelp();
                    return;
                }
                
                
                VSChangeColNameop op;
                op.setParameters(appParameters);
                iter= appParameters.find("file");
                if(iter==appParameters.end()){
                    std::cerr<<"No input file table is provided!"<<std::endl;
                    return;
                }//close if
                std::stringstream sFilename(iter->second);
                //appParameters.erase(iter);
                sFilename>>filename;
                if(filename.find(".bin") == std::string::npos) filename.append(".bin");
                VSTable table(filename);
                if(!table.tableExist()){
                    std::cerr<<"No valid input file table is provided"<<std::endl;
                    return;
                }
                
                op.addInput(&table);
                op.execute();
                valOutFilename=op.realOutFilename();
                if(m_historyEnabled)
                {
                    const char* m_historyFile ="hist.xml";
                    
                    iter =appParameters.find("historyfile");
                    if( iter != appParameters.end())
                    {
                        m_historyFile = appParameters.find("historyfile")->second.c_str();
                    }
                    
                    op.writeHistory(m_historyFile,m_opName,appParameters,valOutFilename);
                }

            }
            break;
        }
            /******* AhfHaloGalaxyExt */
        case 29:
        {
#ifdef AHF
            //serial
            if(rank==0)  //serialized
            {
                
                iter =appParameters.find("help");
                if( iter != appParameters.end()){
                    VSAhfHaloGalaxyExtOp op;
                    op.printHelp();
                    return;
                }
                
                
                VSAhfHaloGalaxyExtOp op;
                op.setParameters(appParameters);
                iter= appParameters.find("file");
                if(iter==appParameters.end()){
                    std::cerr<<"No input file table is provided!"<<std::endl;
                    return;
                }//close if
                std::stringstream sFilename(iter->second);
               // appParameters.erase(iter);
                sFilename>>filename;
                if(filename.find(".bin") == std::string::npos) filename.append(".bin");
                VSTable table(filename);
                if(!table.tableExist()){
                    std::cerr<<"No valid input file table is provided"<<std::endl;
                    return;
                }
                
                op.addInput(&table);
                op.execute();
                valOutFilename=op.realOutFilename();
                if(m_historyEnabled)
                {
                    const char* m_historyFile ="hist.xml";
                    
                    iter =appParameters.find("historyfile");
                    if( iter != appParameters.end())
                    {
                        m_historyFile = appParameters.find("historyfile")->second.c_str();
                    }
                    
                    op.writeHistory(m_historyFile,m_opName,appParameters,valOutFilename);
                }

            }
#else
            if(rank==0) std::cerr<<"Operation not allowed. Compile VisIVO with AHF option";
            return;
#endif
            
            break;
        }
            /******* SplitTable */
        case 30:
        {
            //serial
            if(rank==0)  //serialized
            {
                
                iter =appParameters.find("help");
                if( iter != appParameters.end()){
                    VSSplitTableOp op;
                    op.printHelp();
                    return;
                }
                
                
                VSSplitTableOp op;
                op.setParameters(appParameters);
                iter= appParameters.find("file");
                if(iter==appParameters.end()){
                    std::cerr<<"No input file table is provided!"<<std::endl;
                    return;
                }//close if
                std::stringstream sFilename(iter->second);
                //appParameters.erase(iter);
                sFilename>>filename;
                if(filename.find(".bin") == std::string::npos) filename.append(".bin");
                VSTable table(filename);
                if(!table.tableExist()){
                    std::cerr<<"No valid input file table is provided"<<std::endl;
                    return;
                }
                
                op.addInput(&table);
                op.execute();
                valOutFilename=op.realOutFilename();
                if(m_historyEnabled)
                {
                    const char* m_historyFile ="hist.xml";
                    
                    iter =appParameters.find("historyfile");
                    if( iter != appParameters.end())
                    {
                        m_historyFile = appParameters.find("historyfile")->second.c_str();
                    }
                    
                    op.writeHistory(m_historyFile,m_opName,appParameters,valOutFilename);
                }

            }
            break;
        }
            /******* wrvotrable */
        case 31:
        {
            //serial
            if(rank==0)  //serialized
            {
                
                iter =appParameters.find("help");
                if( iter != appParameters.end()){
                    VSWriteVotableOp op;
                    op.printHelp();
                    return;
                }
                
                
                VSWriteVotableOp op;
                op.setParameters(appParameters);
                iter= appParameters.find("file");
                if(iter==appParameters.end()){
                    std::cerr<<"No input file table is provided!"<<std::endl;
                    return;
                }//close if
                std::stringstream sFilename(iter->second);
                //appParameters.erase(iter);
                sFilename>>filename;
                if(filename.find(".bin") == std::string::npos) filename.append(".bin");
                VSTable table(filename);
                if(!table.tableExist()){
                    std::cerr<<"No valid input file table is provided"<<std::endl;
                    return;
                }
                
                op.addInput(&table);
                op.execute();
                valOutFilename=op.realOutFilename();
                if(m_historyEnabled)
                {
                    const char* m_historyFile ="hist.xml";
                    
                    iter =appParameters.find("historyfile");
                    if( iter != appParameters.end())
                    {
                        m_historyFile = appParameters.find("historyfile")->second.c_str();
                    }
                    
                    op.writeHistory(m_historyFile,m_opName,appParameters,valOutFilename);
                }

            }
            break;
        }
            /******* cut */
        case 32:
        {
            //serial
            if(rank==0)  //serialized
            {
                
                iter =appParameters.find("help");
                if( iter != appParameters.end()){
                    VSCutOp op;
                    op.printHelp();
                    return;
                }
                
                
                VSCutOp op;
                op.setParameters(appParameters);
                iter= appParameters.find("file");
                if(iter==appParameters.end()){
                    std::cerr<<"No input file table is provided!"<<std::endl;
                    return;
                }//close if
                std::stringstream sFilename(iter->second);
                //appParameters.erase(iter);
                sFilename>>filename;
                if(filename.find(".bin") == std::string::npos) filename.append(".bin");
                VSTable table(filename);
                if(!table.tableExist()){
                    std::cerr<<"No valid input file table is provided"<<std::endl;
                    return;
                }
                
                op.addInput(&table);
                op.execute();
                valOutFilename=op.realOutFilename();
                if(m_historyEnabled)
                {
                    const char* m_historyFile ="hist.xml";
                    
                    iter =appParameters.find("historyfile");
                    if( iter != appParameters.end())
                    {
                        m_historyFile = appParameters.find("historyfile")->second.c_str();
                    }
                    
                    op.writeHistory(m_historyFile,m_opName,appParameters,valOutFilename);
                }

            }
            break;
        }
            
            /******* vollimit */
        case 33:
        {
            //serial
            if(rank==0)  //serialized
            {
                
                iter =appParameters.find("help");
                if( iter != appParameters.end()){
                    VSVolLimitOp op;
                    op.printHelp();
                    return;
                }
                
                
                VSVolLimitOp op;
                op.setParameters(appParameters);
                iter= appParameters.find("file");
                if(iter==appParameters.end()){
                    std::cerr<<"No input file table is provided!"<<std::endl;
                    return;
                }//close if
                std::stringstream sFilename(iter->second);
                //appParameters.erase(iter);
                sFilename>>filename;
                if(filename.find(".bin") == std::string::npos) filename.append(".bin");
                VSTable table(filename);
                if(!table.tableExist()){
                    std::cerr<<"No valid input file table is provided"<<std::endl;
                    return;
                }
                
                op.addInput(&table);
                op.execute();
                valOutFilename=op.realOutFilename();
                if(m_historyEnabled)
                {
                    const char* m_historyFile ="hist.xml";
                    
                    iter =appParameters.find("historyfile");
                    if( iter != appParameters.end())
                    {
                        m_historyFile = appParameters.find("historyfile")->second.c_str();
                    }
                    
                    op.writeHistory(m_historyFile,m_opName,appParameters,valOutFilename);
                }
            }
            break;
        }
            /******* include */
        case 34:
        {
            //serial
            if(rank==0)  //serialized
            {
                
                iter =appParameters.find("help");
                if( iter != appParameters.end()){
                    VSIncludeOp op;
                    op.printHelp();
                    return;
                }
                VSIncludeOp op;
                op.setParameters(appParameters);
                iter=appParameters.find("file");
                if(iter==appParameters.end())
                {
                    std::cerr<<"No input file table is provided!"<<std::endl;
                    return;
                }
                std::stringstream sFilename(iter->second);
                sFilename>>filename;
                if(filename.find(".bin") == std::string::npos) filename.append(".bin");
                VSTable table(filename);
                if(!table.tableExist())
                {
                    std::cerr<<"No valid input file table is provided"<<std::endl;
                    return;
                }
                op.addInput(&table);
                op.execute();
                valOutFilename=op.realOutFilename();
                
                if(m_historyEnabled)
                {
                    const char* m_historyFile ="hist.xml";
                    
                    iter =appParameters.find("historyfile");
                    if( iter != appParameters.end())
                    {
                        m_historyFile = appParameters.find("historyfile")->second.c_str();
                    }
                    
                    op.writeHistory(m_historyFile,m_opName,appParameters,valOutFilename);
                }

                
            }
            break;
        }
            /*** MRcampos **/
        case 35:
        {
            //serial
            if(rank==0)  //serialized
            {
                
                iter =appParameters.find("help");
                if( iter != appParameters.end())
                {
                    VSMRCamPos op;
                    op.printHelp();
                    return;
                }
                
                iter =appParameters.find("file");
                if( iter == appParameters.end())
                {
                    std::cerr <<"No input file table is provided"<<std::endl;
                    return;
                }
                std::stringstream sFilename(iter->second);
                //appParameters.erase(iter);
                sFilename>>filename;
                if(filename.find(".bin") == std::string::npos)
                    filename.append(".bin");
                VSTable table(filename);
                if(!table.tableExist())
                {
                    std::cerr <<"No valid input file table is provided"<<std::endl;
                    return;
                }
                
                VSMRCamPos op;
                op.setParameters(appParameters);
                op.addInput(&table);
                op.execute();
                valOutFilename=op.realOutFilename();
                if(m_historyEnabled)
                {
                    const char* m_historyFile ="hist.xml";
                    
                    iter =appParameters.find("historyfile");
                    if( iter != appParameters.end())
                    {
                        m_historyFile = appParameters.find("historyfile")->second.c_str();
                    }
                    
                    op.writeHistory(m_historyFile,m_opName,appParameters,valOutFilename);
                }

            }
            break;
        }
            /*** clusterCell ** **QUI da fare!!!
             case 36:
             {
             //serial
             if(rank==0)  //serialized
             {
             iter =appParameters.find("help");
             if( iter != appParameters.end())
             {
             VSClusterCellOp op;
             op.printHelp();
             return;
             }
             
             iter =appParameters.find("file");
             if( iter == appParameters.end())
             {
             std::cerr <<"No input file table is provided"<<std::endl;
             return;
             }
             std::stringstream sFilename(iter->second);
             sFilename>>filename;
             if(filename.find(".bin") == std::string::npos)
             filename.append(".bin");
             VSTable table(filename);
             if(!table.tableExist())
             {
             std::cerr <<"No valid input file table is provided"<<std::endl;
             return;
             }
             
             
             VSClusterCellOp op;
             op.setParameters(appParameters);
             op.addInput(&table);
             op.execute();
             valOutFilename=op.realOutFilename();
             }
             break;
             }
             */
            /*** muportal **/
        case 37:
        {
            //serial
            if(rank==0)  //serialized
            {
                iter =appParameters.find("help");
                if( iter != appParameters.end())
                {
                    VSMuPortalOp op;
                    op.printHelp();
                    return;
                }
                
                iter =appParameters.find("file");
                if( iter == appParameters.end())
                {
                    std::cerr <<"No input file table is provided"<<std::endl;
                    return;
                }
                std::stringstream sFilename(iter->second);
                sFilename>>filename;
                if(filename.find(".bin") == std::string::npos)
                    filename.append(".bin");
                VSTable table(filename);
                if(!table.tableExist())
                {
                    std::cerr <<"No valid input file table is provided"<<std::endl;
                    return;
                }
                
                
                VSMuPortalOp op;
                op.setParameters(appParameters);
                op.addInput(&table);
                op.execute();
                               
                valOutFilename=op.realOutFilename();
                if(m_historyEnabled)
                {
                    const char* m_historyFile ="hist.xml";
                    
                    iter =appParameters.find("historyfile");
                    if( iter != appParameters.end())
                    {
                        m_historyFile = appParameters.find("historyfile")->second.c_str();
                    }
                    
                    op.writeHistory(m_historyFile,m_opName,appParameters,valOutFilename);
                }

            }
            break;
        }
            
            /*** Default **/
        default:
        {
            //serial
            if(rank==0)  //serialized
            {
                
                std::cerr <<"No valid operation was given"<<std::endl;
                std::cerr<<"VisIVOFilters version Version 2.1.1 June 28th 2013 "<<std::endl<<std::endl;
                std::cerr <<"Syntax1: VisIVOFilters --op operation  [PARAMETERS] [--help]"<<std::endl;
                std::cerr <<"Syntax2: VisIVOFilters parameterFile"<<std::endl;
                std::cerr <<"An operation code is expected: randomizer selcolumns merge append selfield mathop decimator extraction visualop showtable statistic pointdistribute pointproperty coarsevolume extractsubvolume interpolate module sigmacontours cartesian2polar AHFstep vbt2ahf grid2point extractlist addId ahfhalolist ahfhalogalaxyext changecolname splittable wrvotable include mres poca"<<std::endl<<std::endl;
                           
            }
            return;
            
            
        }
            
            
            /*** END Default  OP **/
    }
    

    
    
#ifdef GLITE
    // --append
    if(rank==0)  //QUI VERIFICA
    {
        switch(idOp)
        {
            case 20:
            case 23:
            case 34:
            case 6:
            case 18:
            case 16:
            {
                std::map<std::string,std::string>::iterator iterappend;  
                iterappend=appParameters.find("append");
                if(iterappend != appParameters.end() && !fileIsLocal) //lfn file + --append
                {
                    iterappend=appParameters.find("lfnout");
                    if(iterappend!=appParameters.end()) appParameters.erase(iterappend);
                    appParameters.insert(make_pair("lfnout", inputlfn));
                }
            }
                
        } // end switch
        
        iter = appParameters.find("lfnout");
        std::string lfnOutRoot;
        
        if(iter != appParameters.end())
            lfnOutRoot=iter->second;
        
        if(localVolFile!=""){ //grid2point --volume
            remove(localVolFile.c_str());
            localVolFile += ".head";
            remove(localVolFile.c_str());
        }
        for(int i=0;i<valOutFilename.size();i++)
        {
            if(iter != appParameters.end())
            {
                appParameters.insert(make_pair("realOutFilename", valOutFilename[i]));
                
                std::map<std::string,std::string>::iterator iter2;
                iter2=appParameters.find("realOutFilename");
                
                if(idOp==17){ 
                    size_t pos1 = valOutFilename[i].find("interpolate_") + 12;
                    size_t pos2 = valOutFilename[i].find(".bin");
                    size_t size = pos2 - pos1; 
                    std::string value = valOutFilename[i].substr(pos1, size); 
                    std::string lfnName = lfnOutRoot + "_" + value + ".bin";
                    //appParameters.erase(iter);
                    appParameters.insert(make_pair("lfnout", lfnName));
                }
                
                if(idOp==30){
                    size_t pos1 = valOutFilename[i].find("_split_") + 7;
                    size_t pos2 = valOutFilename[i].find(".bin");
                    size_t size = pos2 - pos1;
                    std::string value = valOutFilename[i].substr(pos1, size);
                    std::string lfnName = lfnOutRoot + "_split_" + value + ".bin";
                    //appParameters.erase(iter);
                    appParameters.insert(make_pair("lfnout", lfnName));
                }
                
                gLInterface.writeVF(appParameters, idOp);
                appParameters.erase(iter2);
                
                remove(valOutFilename[i].c_str());
                
                if(idOp !=10 && idOp!=11 && idOp!=31 && idOp!=33) // --out is a VBT and not ascii
                {
                    std::string tmpfile=valOutFilename[i]+".head";
                    remove(tmpfile.c_str());
                }
            }
        }  
        
        if(idOp==17) // infiles option
        {
            std::stringstream ssInfileparameters;
            ssInfileparameters.str(localFilename);
            std::string f1,f2;
            ssInfileparameters>>f1;
            ssInfileparameters>>f2;
            if(file1.substr(0,6) =="lfn://") remove(f1.c_str());
            if(file2.substr(0,6) =="lfn://") remove(f2.c_str());
            
        } else{
            if(localFilename!="" && !fileIsLocal){  
                remove(localFilename.c_str());
                if(idOp!=33) // --out is a VBT and not ascii
                {
                    std::string tmpfile=localFilename+".head";
                    remove(tmpfile.c_str());
                }
            }
        }
        gLInterface.rmv();
    } //if(rank==0)
#endif
    
    return;
}

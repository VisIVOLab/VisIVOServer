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

#include "commandline.h"
#include <cstdlib>
#include <cstring>

#include <iostream>
#include <sstream>
#include <fstream>

#include "visivoutils.h"
#include "asciisource.h"
#include "csvsource.h"
#include "binsource.h"
#include "flysource.h"
#ifdef HAS_CHANGA_IMPORTER
#include "changasource.h"
#endif
#include "vosourcenew.h"
#include "gadgetsource.h"

#ifndef LIGHT
  #include "xmlsource.h"
  #include "vosource.h"
#endif

#include "fitstablesource.h"
#include "fitsimagesource.h"
#include "hdf5source.h"
#include "rawgridsource.h"
#include "rawpointssource.h"
#include "muportalsource.h"
#include "ramsessource.h"
#ifdef WIN32
	#include <io.h>
#else
	#include <unistd.h>
#endif
//---------------------------------------------------------------------
CommandLine::CommandLine ( )
//---------------------------------------------------------------------
{
   
	m_type="notype";
	m_currentPath="nopath";
	m_binaryDir="./";
	m_binaryName="VisIVOServerBinary.bin";
	m_binaryPath="";
	m_endian="little";
	m_dataType="float";
	m_out="nopath";
	m_login="nouser";
  	m_remoteFile="noremote";
	m_npoints=0;
	m_file="table";
  	m_binaryHeader="noheader";
    m_historyEnabled=false;
    m_historyFile="hist.xml";

	m_size[0]=m_size[1]=m_size[2]=1;
	m_comput[0]=m_comput[1]=m_comput[2]=0;
    m_missing=-1.0918273645e+23;
 	m_text=-1.4536271809e+15;
	m_lfn="";
	m_VO="";
	m_gLiteOut=false;
	m_se="";
	m_outPath="";
	m_datasetList.clear();
    m_fitshdunum=-1;


}
//---------------------------------------------------------------------
CommandLine::~CommandLine ( )
//---------------------------------------------------------------------
{
 

}
//---------------------------------------------------------------------
int CommandLine::parseOption (const std::vector<std::string>  arguments )
//---------------------------------------------------------------------
{
  
  
  
	std::stringstream ss;
	int i;
	for (i=0; i<arguments.size(); i++)
	{
		if (arguments[i]=="--missingvalue")
		{
/*      		  std::string ckInput=arguments[i+1];
      		  if(ckInput.rfind('--')==0)
      		  {
        	    std::cerr<<"Error on "<<arguments[i]<< "argument: "<<ckInput<<std::endl;
		    return -1;
      		  }	*/
		std::stringstream tmp; 
		tmp<<arguments[++i];
		tmp>>m_missing;
		}
		else if (arguments[i]=="--textvalue")
		{
/*      		  std::string ckInput=arguments[i+1];
      		  if(ckInput.rfind('--')==0)
      		  {
        	    std::cerr<<"Error on "<<arguments[i]<< "argument: "<<ckInput<<std::endl;
		    return -1;
      		  }	*/
		std::stringstream tmp; 
		tmp<<arguments[++i];
		tmp>>m_text;
		}
		else if (arguments[i]=="--userpwd")
		{
      		  std::string ckInput=arguments[i+1];
     		  if(ckInput.find_first_of('-')==0)
      		  {
        	    std::cerr<<"Error on "<<arguments[i]<< " argument: "<<ckInput<<std::endl;
		    return -1;
      		  }	
		m_login=arguments[++i];
		}
  		else if (arguments[i]=="--volume")
		{
			m_file="volume";
		}
        
		else if (arguments[i]=="--fformat")
		{
      		  std::string ckInput=arguments[i+1];
      		  if(ckInput.find_first_of('-')==0)
      		  {
        	    std::cerr<<"Error on "<<arguments[i]<< " argument: "<<ckInput<<std::endl;
		    return -1;
      		  }	
			m_type=arguments[++i];
		} 		
		else if (arguments[i]=="--fields")
		{
      		  std::string ckInput=arguments[i+1];
      		  if(ckInput.find_first_of('-')==0)
      		  {
        	    std::cerr<<"Error on "<<arguments[i]<< " argument: "<<ckInput<<std::endl;
		    	return -1;
      		  }	
			  else{
				while(arguments[i+1].find_first_of('-')!=0){
					m_fields.push_back(arguments[i+1]);
					i++;
				}
			  }
		} 
		else if (arguments[i]=="--binaryheader")
		{
      		  std::string ckInput=arguments[i+1];
      		  if(ckInput.find_first_of('-')==0)
      		  {
        	    std::cerr<<"Error on "<<arguments[i]<< " argument: "<<ckInput<<std::endl;
		    return -1;
      		  }	
			m_binaryHeader=arguments[++i];
		}
        else if (arguments[i]=="--fitshdunum")
		{
            std::string ckInput=arguments[i+1];
            if(ckInput.find_first_of('-')==0)
            {
        	    std::cerr<<"Error on "<<arguments[i]<< " argument: "<<ckInput<<std::endl;
                return -1;
            }
			std::stringstream ss1;
            
			ss1<<arguments[++i];
			ss1>>m_fitshdunum ;
		}

    
		else  if (arguments[i]=="--bigendian")
		{
			m_endian="big";
		}
 
    
		else if(arguments[i]=="--out")
		{
		  if(i==arguments.size()-2)
		  {
				std::cerr<<"Output file and input file cannot be have the same filename"<<std::endl;
				return -1;
		  }
      		  std::string ckInput=arguments[i+1];
      		  if(ckInput.find_first_of('-')==0)
      		  {
        	    std::cerr<<"Error on "<<arguments[i]<< " argument: "<<ckInput<<std::endl;
		    return -1;
      		  }	
		  m_out=arguments[++i];
		}

		else if (arguments[i]=="--sizex")
		{
      		  std::string ckInput=arguments[i+1];
      		  if(ckInput.find_first_of('-')==0)
      		  {
        	    std::cerr<<"Error on "<<arguments[i]<< " argument: "<<ckInput<<std::endl;
		    return -1;
      		  }	
			std::stringstream ss1;
      
			ss1<<arguments[++i];
			ss1>>m_size[0];
		}
    
		else if (arguments[i]=="--sizey")
		{
      		  std::string ckInput=arguments[i+1];
      		  if(ckInput.find_first_of('-')==0)
      		  {
        	    std::cerr<<"Error on "<<arguments[i]<< " argument: "<<ckInput<<std::endl;
		    return -1;
      		  }	
			std::stringstream ss2;
      
			ss2<<arguments[++i];
			ss2>>m_size[1];
		}
    
		else if (arguments[i]=="--sizez")
		{
      		  std::string ckInput=arguments[i+1];
      		  if(ckInput.find_first_of('-')==0)
      		  {
        	    std::cerr<<"Error on "<<arguments[i]<< " argument: "<<ckInput<<std::endl;
		    return -1;
      		  }	
			std::stringstream ss3;
      
			ss3<<arguments[++i];
			ss3>>m_size[2];
		}
    
		else if (arguments[i]=="--compx")
		{
      		  std::string ckInput=arguments[i+1];
      		  if(ckInput.find_first_of('-')==0)
      		  {
        	    std::cerr<<"Error on "<<arguments[i]<< " argument: "<<ckInput<<std::endl;
		    return -1;
      		  }	
			std::stringstream ss4;
      
			ss4<<arguments[++i];
			ss4>>m_comput[0];
		}
    
		else if (arguments[i]=="--compy")
		{
      		  std::string ckInput=arguments[i+1];
      		  if(ckInput.find_first_of('-')==0)
      		  {
        	    std::cerr<<"Error on "<<arguments[i]<< " argument: "<<ckInput<<std::endl;
		    return -1;
      		  }	
			std::stringstream ss5;
      
			ss5<<arguments[++i];
			ss5>>m_comput[1];
		}
    
		else if (arguments[i]=="--compz")
		{
      		  std::string ckInput=arguments[i+1];
      		  if(ckInput.find_first_of('-')==0)
      		  {
        	    std::cerr<<"Error on "<<arguments[i]<< " argument: "<<ckInput<<std::endl;
		    return -1;
      		  }	
			std::stringstream ss6;
      
			ss6<<arguments[++i];
			ss6>>m_comput[2];
		}
    
		else if (arguments[i]=="--npoints")
		{
      		  std::string ckInput=arguments[i+1];
      		  if(ckInput.find_first_of('-')==0)
      		  {
        	    std::cerr<<"Error on "<<arguments[i]<< " argument: "<<ckInput<<std::endl;
		    return -1;
      		  }	
			std::stringstream ss6;
      
			ss6<<arguments[++i];
			ss6>>m_npoints;
		}
		else if (arguments[i]=="--datasetlist")
		{
      		  std::string ckInput=arguments[i+1];
      		  if(ckInput.find_first_of('-')==0)
      		  {
        	    std::cerr<<"Error on "<<arguments[i]<< " argument: "<<ckInput<<std::endl;
		    return -1;
      		  }	
			while(true)
			{  
			    m_datasetList=m_datasetList+arguments[++i]+" ";
			    ckInput=arguments[i+1];
			    if(ckInput.find_first_of('-')==0 || i==arguments.size()-2)
			      break;			    
			}
		}
		else if (arguments[i]=="--hyperslab")
		{
      		  std::string ckInput=arguments[i+1];
      		  if(ckInput.find_first_of('-')==0)
      		  {
        	    std::cerr<<"Error on "<<arguments[i]<< " argument: "<<ckInput<<std::endl;
		    return -1;
      		  }	
			int tmpCount=0;
			std::string hyperString;
			while(true)
			{  
			    hyperString.append(arguments[++i]);
			    hyperString.append(" ");
			    tmpCount++;
			    if(tmpCount==3)
			    {
			      tmpCount=0;
			       m_hyperslab.push_back(hyperString); 
			       hyperString.erase();
			    }
			    ckInput=arguments[i+1];
			    if(ckInput.find_first_of('-')==0 || i==arguments.size()-2)
			    {
			       if(tmpCount!=0) 
				 m_hyperslab.push_back(hyperString);     
			      break;
			    }
			}
		}
        
		else  if (arguments[i]=="--double")
		{
			m_dataType="double";
		}   
		else  if (arguments[i]=="--VO")
		{
      		  std::string ckInput=arguments[i+1];
      		  if(ckInput.find_first_of('-')==0)
      		  {
        	    std::cerr<<"Error on "<<arguments[i]<< " argument: "<<ckInput<<std::endl;
		    return -1;
      		  }	
			m_VO=arguments[++i];		  
		}   
		else  if (arguments[i]=="--lfnout")
		{
      		  std::string ckInput=arguments[i+1];
      		  if(ckInput.find_first_of('-')==0)
      		  {
        	    std::cerr<<"Error on "<<arguments[i]<< " argument: "<<ckInput<<std::endl;
		    return -1;
      		  }	
			m_outlfn=arguments[++i];
		        m_gLiteOut=true;
		}   
		else  if (arguments[i]=="--se")
		{
      		  std::string ckInput=arguments[i+1];
      		  if(ckInput.find_first_of('-')==0)
      		  {
        	    std::cerr<<"Error on "<<arguments[i]<< " argument: "<<ckInput<<std::endl;
		    return -1;
      		  }	
			m_se=arguments[++i];        
        
        }
        
        else  if (arguments[i]=="--history")
        {
            m_historyEnabled=true;
        }
        else if (arguments[i]=="--historyfile")
		{
           
            m_historyEnabled=true;
            std::string ckInput=arguments[i+1];
            if(ckInput.find_first_of('-')==0)
            {
                std::cerr<<"Error on "<<arguments[i]<< " argument: "<<ckInput<<std::endl;
                return -1;
            }
            m_historyFile=arguments[++i];
            
		}
 
	}
	m_currentPath=arguments[(arguments.size()-1)]; //!filename including path!!
//check for gLite lfn
	if(m_gLiteOut && m_outlfn.substr(0,6) != "lfn://")
	{
	   std::cerr<<"Invalid output logical file name: "<<m_outlfn<<std::endl;
	   return -1;
	} 

	m_datasetList=trim(m_datasetList);  //trim for multiple options
        if(m_out !="nopath")
	{   
		std::ofstream file;
		std::string outDirTest=getDir(m_out);
		outDirTest.append("testOutputVisIVOServerO_testfile");
		file.open( outDirTest.c_str(), std::ios_base::out );
		if (!file.is_open())
		{
				std::cerr<<"Output file directory does not exist"<<std::endl;
				return -1;
		}
		remove(outDirTest.c_str());
	}


	if(m_currentPath.substr(0,7) == "http://" ||m_currentPath.substr(0,7) == "sftp://" || m_currentPath.substr(0,6) == "ftp://")
	{
		bool remoteOp;
		std::string downloadedPath, remoteFile;
		remoteFile=m_currentPath;
		if(m_out !="nopath") downloadedPath=getDir(m_out);
		downloadedPath.append(getName(m_currentPath));
		remoteOp=remoteDownloadFiles(m_currentPath,m_login,downloadedPath);
		if(!remoteOp)
		{
			std::cerr<<"Remote download failed"<<std::endl;
			return -1;
		}
		m_currentPath=downloadedPath;
		m_remoteFile=downloadedPath;
		if((m_type=="binary" && m_binaryHeader=="noheader") ||(m_type=="binary" && (m_binaryHeader.substr(0,7) == "http://" ||m_binaryHeader.substr(0,7) == "sftp://" || m_binaryHeader.substr(0,6) == "ftp://") ))
		{
			std::string headerFile, downloadedHeader;
			if(m_binaryHeader=="noheader")
			{
				headerFile=remoteFile;
				headerFile=headerFile.append(".head");
				downloadedHeader=m_remoteFile;
				downloadedHeader=downloadedHeader.append(".head");
			} else
			{
				headerFile=m_binaryHeader;
				if(m_out !="nopath") downloadedHeader=getDir(m_out);
				downloadedHeader=downloadedHeader.append(getName(m_binaryHeader));
				m_binaryHeader=downloadedHeader;
			}
			remoteOp=remoteDownloadFiles(headerFile,m_login,downloadedHeader);
			if(!remoteOp)
			{
				std::cerr<<"Remote header download failed"<<std::endl;
				return -1;
			}
		}
	}
	if(m_binaryHeader.substr(0,7) == "http://" ||m_binaryHeader.substr(0,7) == "sftp://" || m_binaryHeader.substr(0,6) == "ftp://")
	{
		bool remoteOp;
		std::string headerFile, downloadedHeader;
		headerFile=m_binaryHeader;
		if(m_out !="nopath") downloadedHeader=getDir(m_out);
		downloadedHeader=downloadedHeader.append(getName(m_binaryHeader));
		m_binaryHeader=downloadedHeader;
		remoteOp=remoteDownloadFiles(headerFile,m_login,downloadedHeader);
		if(!remoteOp)
		{
			std::cout<<"Remote header download failed"<<std::endl;
			return -1;
		}
	
	}

//download gLite    
	if(m_currentPath.substr(0,6) == "lfn://")
	{
		bool remoteOp;
		std::string downloadedPath, remoteFile;
		remoteFile=m_currentPath;
		if(m_out !="nopath") downloadedPath=getDir(m_out);
		downloadedPath.append(getName(m_currentPath));
		remoteOp=remoteDownloadFiles(remoteFile,m_VO,downloadedPath);
		if(!remoteOp)
		{
			std::cerr<<"Remote download failed"<<std::endl;
			return -1;
		}
		m_currentPath=downloadedPath;
		m_remoteFile=downloadedPath;
		if(m_type=="binary")
		{
		    downloadedPath.erase();
		    remoteFile.erase();
		    remoteFile=m_currentPath+".head";
		    if(m_out !="nopath") downloadedPath=getDir(m_out);
		    downloadedPath.append(getName(m_currentPath));
		    downloadedPath=downloadedPath+".head";
		    remoteOp=remoteDownloadFiles(remoteFile,m_VO,downloadedPath);
		    if(!remoteOp)
		    {
			std::cerr<<"Remote header download failed"<<std::endl;
			return -1;
		    }

		}  
	}
//
	if(m_file=="volume")
		if((m_comput[0]<=0 || m_comput[1]<=0 || m_comput[2]<=0) && (m_type!="binary"))
		{
			if(m_type!="hdf5" && m_type!="fitsimage")
			{
			  std::cerr<<"A volume must  have valid x y and z mesh point dimensions"<<std::endl;
			  return -1;
			}

		}
  
	if(m_type!="fitstable" && m_type!="ramses" && m_type!="fitsimage")
	{ 
		std::ifstream inFile;
		inFile.open(m_currentPath.c_str());

		if(!inFile)
		{
			std::cerr<<"the path is incorrect, please try again  or --help for help \n";
			return -1;
		}

		else
			inFile.close();
	}
    
    
    
	return 0;
}


//---------------------------------------------------------------------
int CommandLine::loadFile ()
//---------------------------------------------------------------------
{
	if(m_currentPath=="nopath")
	{
		std::cerr<<"Please give the path \n";
		return -1;
	}
  
	if (m_out!="nopath")
	{
		if (m_out.find("/")!=std::string::npos)
		{

			m_binaryDir=getDir(m_out);
			if (m_binaryDir.at(0)!='/')
				m_binaryDir="./"+ m_binaryDir;

			m_binaryName=getName(m_out);
		}

		else
			m_binaryName=m_out;
	}

	if(m_binaryName.find(".bin") != std::string::npos)
		m_binaryPath=m_binaryDir+m_binaryName;
	else
		m_binaryPath=m_binaryDir+m_binaryName+".bin";
 

	m_ext=getExt(m_currentPath.c_str());
  
	if (m_type=="notype")
	{
		std::cerr <<"please give the file format --fformat [file format] \n";
		return -1;
	}
	else
	{ 
		AbstractSource* pSource;
  
		if ( m_type=="ascii")
			pSource = new AsciiSource();

   
		else if ( m_type=="binary")
			pSource = new BinSource();
 
		else if (m_type=="fitstable")
			pSource = new FitsTableSource();
   
		else if(m_type=="fitsimage")
        {
			pSource = new FitsImageSource();
        }
		else if(m_type=="csv" )
			pSource = new CSVSource();

   
		else if(m_type=="fly" )
			pSource = new FlySource();
 
    
		else if(m_type=="votablefast")
			pSource = new VOSourcenew();

//		else if(m_type=="votable")
//			pSource = new VOSource();
		else if(m_type=="votable")
			pSource = new VOSourcenew();
 
  
		else if(m_type=="gadget")
			pSource = new GadgetSource();

		else if(m_type=="changa"){
#ifdef HAS_CHANGA_IMPORTER
			pSource = new ChangaSource();
#else
    		std::clog << "VisIVOServer Changa Importer not enabled at compile time" << std::endl;
			return -1;
#endif
		}

#ifndef LIGHT
		else if(m_type=="xml")
			pSource = new XmlSource();
#endif
		else if(m_type=="hdf5")
			pSource = new HDF5Source();
     
    
		else if(m_type=="rawgrids")
		{
			pSource = new RawGridSource();
			m_file="volume";
		}
		else if(m_type=="rawpoints")
			pSource = new RawPointsSource();
        
		else if(m_type=="muportal")
			pSource = new MuPortalSource();
		else if(m_type=="ramses"){
		  pSource = new RamsesSource();

		}
		else
		{ 
			std::cerr<<"the format given '"<<m_type<<"' is incorrect, please try again  or --help for help"<<std::endl;
			return -1;
		}
		pSource->setPointsFileName(m_currentPath.c_str(),m_binaryPath.c_str(),m_file.c_str(),
					   m_size,m_comput,m_file.c_str(),m_endian.c_str(),
					   m_dataType.c_str(),m_npoints,m_login.c_str(),
					   m_binaryHeader.c_str(),m_missing,m_text,m_datasetList,
					   m_hyperslab,m_fitshdunum, m_fields);
		if(pSource->readHeader()==0)
			pSource->readData();
        if(m_historyEnabled)
        {
            pSource->writeHistory(m_historyFile.c_str(),m_type.c_str(), m_out.c_str(), m_file.c_str(), m_comput, m_size, m_login.c_str(), m_binaryHeader.c_str(), m_missing, m_text, m_endian.c_str(), m_dataType.c_str(), m_npoints, m_VO.c_str(), m_se.c_str(), m_outlfn.c_str(),m_currentPath.c_str());
        }
        
        
		delete pSource;
		if(m_gLiteOut)
		{
		  bool isvbt=true;
		  if(m_outlfn.find(".bin") == std::string::npos)
		    m_outlfn=m_outlfn+".bin";
		  bool oplfn=remoteLfn(m_binaryPath,m_se,m_VO,m_outlfn,isvbt);
		  remove(m_binaryPath.c_str());
		  if(!oplfn)  return -1;
		}
	}
  
	return 0;
}
 
 //---------------------------------------------------------------------
void CommandLine::showHelp ()
//---------------------------------------------------------------------
 
{
	std::cout<<std::endl;
   	std::cout<<"VisIVOImporter Version 3.5.0 March 19th 2025 "<<std::endl<<std::endl;
  
	std::cout<<" --fformat   [typefile]  (mandatory) Select file type: ascii, csv, votable, binary, fly, gadget, xml, rawpoints, rawgrids, fitstable, fitsimage, hdf5, muportal, ramses, changa"<<std::endl<<std::endl;

	std::cout<<"[pathfile] (mandatory) Absolute path file. Path must be the last command( /home/user/myfile.ascii)"<<std::endl<<std::endl;

	std::cout<<"--out [filename]   (optional) Change default file name  and or directory ( --out /home/user/myfile.bin )  "<<std::endl<<std::endl;

	std::cout<<"--volume  (optional)  if you want a table for volume. It is mandatory if you want select cell size and/or computational cell size "<<std::endl<<std::endl; 
  
	std::cout<<"--sizex [double] --sizey [double]  --sizez [double]  (optional)  if you want  select cell size(--sizex 1 --sizey 1  --sizez 1 ). If you use this command --volume is mandatory .If don't use this commands default size is 1 1 1 "<<std::endl<<std::endl; 
  
	std::cout<<"--compx [double] --compy [double] --compz [double]  (optional)  if you want  select computational cell size. If the mathematical product of this tree values is different from field size the created output will be a table."<<std::endl<<std::endl; 
  
	std::cout<<"--bigendian (optional) use this command only if 'gadget' and 'fly' format file are big endian"<<std::endl<<std::endl;

	std::cout<<"--fields [field_1] [field_2] ... (optional) use this command only when using 'gadget' format to specify which fields you need to import"<<std::endl<<std::endl;
   
	std::cout<<"--double (optional) use this command only if the 'fly' format file has double data type"<<std::endl<<std::endl; 
  
	std::cout<<"--npoints (optional) use this command only to set the number of points in the 'fly' format file "<<std::endl<<std::endl; 

	std::cout<<"--datasetlist  use this command only to select datasets in hdf5 file"<<std::endl<<std::endl; 

	std::cout<<"--hyperslab  use this command only to select dataset hyperslab in hdf5 file"<<std::endl<<std::endl;

    std::cout<<"--history (optional) create an XML file which contains the history of operations performed (default create hist.xml file)"<<std::endl<<std::endl;
    
    std::cout<<"--historyfile [filename]   (optional) Change default history file name  and or directory "<<std::endl<<std::endl; 
}

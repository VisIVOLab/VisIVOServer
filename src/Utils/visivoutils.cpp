/***************************************************************************
 *   Copyright (C) 2008 by Gabriella Caniglia, Marco Comparato *
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

#define M_PI 3.14159265358979323846f
#define M_PI_2 (M_PI/2.F)
#define DEG_TO_RAD_FACTOR (M_PI/180.F)
#include <cstdlib>
#include <cstring>

#include "visivoutils.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cassert>
#include <math.h>
#include "curl/curl.h"

#include <cstring>
#include <cstdlib>
#include <time.h>

#ifdef GLITE
extern "C"{

  #include "lcg_util.h"
}
#endif



static size_t writeData(void *buffer, size_t size, size_t nMemb, void *userP);

//---------------------------------------------------------------------
double masToRad(double mas)
//---------------------------------------------------------------------
{
  double seconds = mas/1000;
  double degs = seconds/3600;

  return degs*DEG_TO_RAD_FACTOR;
}

//---------------------------------------------------------------------
void masToRad(double *mas, int n)
//---------------------------------------------------------------------
{
  int i = 0;

  double seconds = 0;
  double degs = 0;

  for(i = 0; i < n; i++)
  {
    seconds = mas[i]/1000;
    degs = seconds/3600;
    mas[i] = degs * DEG_TO_RAD_FACTOR;
  }

  return;
}

//---------------------------------------------------------------------
void masToRad(float *mas, int n)
//---------------------------------------------------------------------
{
  int i = 0;

  double seconds = 0;
  double degs = 0;

  for(i = 0; i < n; i++)
  {
    seconds = mas[i]/1000;
    degs = seconds/3600;
    mas[i] = degs * DEG_TO_RAD_FACTOR;
  }

  return;
}


//---------------------------------------------------------------------
double degToRad(double deg)
//---------------------------------------------------------------------
{
  return deg*DEG_TO_RAD_FACTOR;
}

//---------------------------------------------------------------------
void degToRad(double *deg, int n)
//---------------------------------------------------------------------
{
  for(int i = 0; i < n; i++)
    deg[i] *= DEG_TO_RAD_FACTOR;

  return;
}
//---------------------------------------------------------------------
void degToRad(float *deg, int n)
//---------------------------------------------------------------------
{
  for(int i = 0; i < n; i++)
    deg[i] *= DEG_TO_RAD_FACTOR;

  return;
}

//---------------------------------------------------------------------
double hmsToRad(char *hms)
//---------------------------------------------------------------------
{
  double hours = 0;
  double minutes = 0;
  double seconds = 0;

  std::stringstream tmpStream;

  tmpStream << hms;

  tmpStream >> hours;
  tmpStream >> minutes;
  tmpStream >> seconds;

  double outcome = hours + (minutes/60) + (seconds/3600);
  outcome = (outcome/24)*360;
  outcome = outcome*DEG_TO_RAD_FACTOR;

  return outcome;
}

//---------------------------------------------------------------------
void hmsToRad(char **hmsIn, double *radOut, int n)
//---------------------------------------------------------------------
{
  double hours = 0;
  double minutes = 0;
  double seconds = 0;

  for(int i = 0; i < n; i++)
  {
    std::stringstream tmpStream;

    tmpStream << hmsIn[i];

    tmpStream >> hours;
    tmpStream >> minutes;
    tmpStream >> seconds;

    radOut[i] = hours + (minutes/60) + (seconds/3600);
    radOut[i] = (radOut[i]/24)*360;
    radOut[i] = radOut[i]*DEG_TO_RAD_FACTOR;
  }

  return;
}
//---------------------------------------------------------------------
void hmsToRad(char **hmsIn, float *radOut, int n)
//---------------------------------------------------------------------
{
  double hours = 0;
  double minutes = 0;
  double seconds = 0;

  for(int i = 0; i < n; i++)
  {
    std::stringstream tmpStream;

    tmpStream << hmsIn[i];

    tmpStream >> hours;
    tmpStream >> minutes;
    tmpStream >> seconds;

    radOut[i] = hours + (minutes/60) + (seconds/3600);
    radOut[i] = (radOut[i]/24)*360;
    radOut[i] = radOut[i]*DEG_TO_RAD_FACTOR;
  }

  return;
}

//---------------------------------------------------------------------
double dmsToRad(char *dms)
//---------------------------------------------------------------------
{
  double degs = 0;
  double minutes = 0;
  double seconds = 0;

  std::stringstream tmpStream;

  tmpStream << dms;

  tmpStream >> degs;
  tmpStream >> minutes;
  tmpStream >> seconds;

  if(degs < 0)
  {
    minutes = -minutes;
    seconds = -seconds;
  }

  double outcome = degs + (minutes/60) + (seconds/3600);
  outcome = outcome*DEG_TO_RAD_FACTOR;

  return outcome;
}


//---------------------------------------------------------------------
void dmsToRad(char **dmsIn, double *radOut, int n)
//---------------------------------------------------------------------
{
  double degs = 0;
  double minutes = 0;
  double seconds = 0;

  for(int i = 0; i < n; i++)
  {
    std::stringstream tmpStream;

    tmpStream << dmsIn[i];

    tmpStream >> degs;
    tmpStream >> minutes;
    tmpStream >> seconds;

    if(degs < 0)
    {
      minutes = -minutes;
      seconds = -seconds;
    }

    radOut[i] = degs + (minutes/60) + (seconds/3600);
    radOut[i] = radOut[i]*DEG_TO_RAD_FACTOR;
  }

  return;
}
//---------------------------------------------------------------------
void dmsToRad(char **dmsIn, float *radOut, int n)
//---------------------------------------------------------------------
{
  double degs = 0;
  double minutes = 0;
  double seconds = 0;

  for(int i = 0; i < n; i++)
  {
    std::stringstream tmpStream;

    tmpStream << dmsIn[i];

    tmpStream >> degs;
    tmpStream >> minutes;
    tmpStream >> seconds;

    if(degs < 0)
    {
      minutes = -minutes;
      seconds = -seconds;
    }

    radOut[i] = degs + (minutes/60) + (seconds/3600);
    radOut[i] = radOut[i]*DEG_TO_RAD_FACTOR;
  }

  return;
}

//---------------------------------------------------------------------
int iCompare(std::string str1, std::string str2)
//---------------------------------------------------------------------
{
#ifdef WIN32
  return strcmpi(str1.c_str(), str2.c_str());
#else
  return strcasecmp(str1.c_str(), str2.c_str());
#endif
}

//---------------------------------------------------------------------
std::string getTempFileName(std::string suffix)
//---------------------------------------------------------------------
{
  std::stringstream tmpFileName;

#ifdef WIN32
  char *tempDir = getenv("TEMP");

  if(!tempDir)
    tempDir = getenv("TMP");

  if(tempDir)
  {
    tmpFileName << tempDir;
    tmpFileName << "\\";
  }
#else

  tmpFileName << "/tmp/";    

#endif

  tmpFileName << "VisIVOTmp" << rand() << suffix;

  return tmpFileName.str();
}

//---------------------------------------------------------------------
std::string getDir(std::string path)
//---------------------------------------------------------------------
{
#ifdef _WIN32
  int idx = path.rfind('\\');
#else
  int idx = path.rfind('/');
#endif
  if(idx >= 0)
    path.erase(idx + 1, path.length() - (idx + 1));
  else
    path = "";

  return path;
}

//---------------------------------------------------------------------
std::string getExt(std::string path)
//---------------------------------------------------------------------
{
  int idx = path.rfind('.');

  if(idx < 0)
    return "";

  return path.erase(0, idx + 1);
}

//---------------------------------------------------------------------
std::string trimRight(const std::string & source, const std::string & t /*= " "*/)
//---------------------------------------------------------------------
{
  std::string str = source;
  return str.erase(str.find_last_not_of(t) + 1);
}

//---------------------------------------------------------------------
std::string trimLeft(const std::string & source, const std::string & t /*= " "*/)
//---------------------------------------------------------------------
{
  std::string str = source;
  return str.erase(0, str.find_first_not_of(t));
}

//---------------------------------------------------------------------
std::string trim(const std::string & source, const std::string & t /*= " "*/)
//---------------------------------------------------------------------
{
  std::string str = source;
  return trimLeft(trimRight(str, t), t);
}

//---------------------------------------------------------------------
void findAndReplace(std::string &str, char find, char replace)
//---------------------------------------------------------------------
{
  size_t i;

  for(; (i = str.find(find)) != std::string::npos;)
    str.replace(i, 1, 1, replace);

  return;
}

//---------------------------------------------------------------------
static size_t writeData(void *buffer, size_t size, size_t nMemb, void *userP)
//---------------------------------------------------------------------
{
  unsigned int byteToWrite = size*nMemb;

  std::ofstream *out = (std::ofstream *)userP;
  out->write((char *)buffer, byteToWrite);

  return byteToWrite;
}

//---------------------------------------------------------------------
std::string getName(std::string path)
//---------------------------------------------------------------------
{
#ifdef _WIN32
  int idx = path.rfind('\\');
#else
  int idx = path.rfind('/');
#endif
  if(idx >= 0)
    path.erase(0, (idx + 1));
//  else
//    path = "";

  return path;
}
//---------------------------------------------------------------------
void makeHeader( unsigned long long  int rows,std::string path,const std::vector<std::string>fields,double size[],double comp[],std::string file)
//---------------------------------------------------------------------
{
    
  int cellComp[3], cellSize[3];
  std::string headerBinaryName = path + ".head";
  //std::clog <<"header="<<headerFileName<<endl;

  std::ofstream outHeader(headerBinaryName.c_str());
  if(!outHeader)  return;

  int cols=fields.size();

  outHeader << "float" << std::endl;
  outHeader << cols << std::endl;
  outHeader << rows ;
  
  //std::clog<<"file="<<file<<std::endl;
  if (file=="volume")
  {
    if (comp[0]==0 && comp[1]==0 && comp[2]==0)
      comp[0]=comp[1]=comp[2]=tryToSetDimension(rows);
    
    if (comp[0]*comp[1]*comp[2]==rows)
    {
      outHeader <<' '<<comp[0]<<' '<<comp[1]<<' '<<comp[2];
      //std::clog<<"comput "<<comp[0]<<' '<<comp[1]<<' '<<comp[2]<<std::endl;
    
      outHeader <<' '<<size[0]<<' '<<size[1]<<' '<<size[2];
     // std::clog<<"size "<<size[0]<<' '<<size[1]<<' '<<size[2]<<std::endl;
    }
  
   
  }
  outHeader<< std::endl;
#ifdef VSBIGENDIAN
  outHeader << "big" << std::endl;
#else
  outHeader << "little" << std::endl;
#endif  

  //std::clog<<"cols="<<cols<<std::endl;
  //std::clog<<"rows="<<rows<<std::endl;
  for(int i = 0; i < cols; i++)
  { 
//     std::clog<<fields[i]<<endl;
    outHeader << fields[i]<< std::endl;
  
  }
  
  
  outHeader.close();
  
  if (file=="volume" && comp[0]*comp[1]*comp[2]!=rows)
  {
    std::cout<< "Grid size and field size don't mach:"<<std::endl<<std::endl;
    std::cout<< "Grid size=compx*compy*compz="<<comp[0]*comp[1]*comp[2]<<std::endl<<std::endl;
    std::cout<<"field size=Rows="<<rows<<std::endl<<std::endl;
    std::cout<< "Your binary file is table";
  }
  
  return;

}

//---------------------------------------------------------------------
void makeHeader( int rows,std::string path,const std::vector<std::string>fields,double size[],double comp[],std::string file)
//---------------------------------------------------------------------
{
    
  int cellComp[3], cellSize[3];
  std::string headerBinaryName = path + ".head";
  //std::clog <<"header="<<headerFileName<<endl;

  std::ofstream outHeader(headerBinaryName.c_str());
  if(!outHeader)  return;

  int cols=fields.size();

  outHeader << "float" << std::endl;
  outHeader << cols << std::endl;
  outHeader << rows ;
  
  //std::clog<<"file="<<file<<std::endl;
  if (file=="volume")
  {
    if (comp[0]==0 && comp[1]==0 && comp[2]==0)
      comp[0]=comp[1]=comp[2]=tryToSetDimension(rows);
    
    if (comp[0]*comp[1]*comp[2]==rows)
    {
      outHeader <<' '<<comp[0]<<' '<<comp[1]<<' '<<comp[2];
      //std::clog<<"comput "<<comp[0]<<' '<<comp[1]<<' '<<comp[2]<<std::endl;
    
      outHeader <<' '<<size[0]<<' '<<size[1]<<' '<<size[2];
     // std::clog<<"size "<<size[0]<<' '<<size[1]<<' '<<size[2]<<std::endl;
    }
  
   
  }
  outHeader<< std::endl;
#ifdef VSBIGENDIAN
  outHeader << "big" << std::endl;
#else
  outHeader << "little" << std::endl;
#endif  

  //std::clog<<"cols="<<cols<<std::endl;
  //std::clog<<"rows="<<rows<<std::endl;
  for(int i = 0; i < cols; i++)
  { 
//     std::clog<<fields[i]<<endl;
    outHeader << fields[i]<< std::endl;
  
  }
  
  
  outHeader.close();
  
  if (file=="volume" && comp[0]*comp[1]*comp[2]!=rows)
  {
    std::cout<< "Grid size and field size don't mach:"<<std::endl<<std::endl;
    std::cout<< "Grid size=compx*compy*compz="<<comp[0]*comp[1]*comp[2]<<std::endl<<std::endl;
    std::cout<<"field size=Rows="<<rows<<std::endl<<std::endl;
    std::cout<< "Your binary file is table";
  }
  
  return;

}

//----------------------------
double tryToSetDimension(int nRows)

//---------------------------
{
  double size=pow((double)nRows,1.0/3);
 
  double dimensionUp=ceil(size);
  //std::clog<<"dimensionUp="<<dimensionUp<<std::endl;
 
  double dimensionDown= floor(size);
  //std::clog<<"dimensionDown="<<dimensionDown<<std::endl;
  
  if(nRows==pow(dimensionUp,3))
    return dimensionUp;
 
  if (nRows==pow(dimensionDown,3))
    return dimensionDown;
 
  else 
    return 0.0;
  
}

//----------------------------
float floatSwap(char *value)

//---------------------------
{

  int size =sizeof(float);
  float swapped;
  char *buffer;
  buffer = new char [sizeof(float)];

  for (int i=0; i<size; i++)
    buffer[ i ] = value[ size-1-i ];

  swapped= *( (float *) buffer );
  delete [] buffer;
  return swapped;

}

//----------------------------
double doubleSwap(char *value)

//---------------------------
{
  int size =sizeof(double);
  char *buffer;
  double swapped;
  buffer = new char[sizeof(double)];

  for (int i=0; i<size; i++)
    buffer[ i ] = value[ size-1-i ];
  

  swapped= *( (double *) buffer );
  delete [] buffer;
  return swapped;

}
//----------------------------
long double longdoubleSwap(char *value)

//---------------------------
{
  int size =sizeof(long double);
  char *buffer;
  long double swapped;
  buffer = new char[sizeof(long double)]; 

  for (int i=0; i<size; i++)
    buffer[ i ] = value[ size-1-i ];
   

  swapped= *( (long double *) buffer );
  delete [] buffer;
  return swapped;
}
//----------------------------
int intSwap(char *value)

//---------------------------
{
  int size =sizeof(int);
  char *buffer;
  int swapped;
  buffer = new char[sizeof(int)];

  for (int i=0; i<size; i++)
    buffer[ i ] = value[ size-1-i ];
  

  swapped= *( (int *) buffer );
  delete [] buffer;
  return swapped;
}
//----------------------------
long int longintSwap(char *value)

//---------------------------
{
  int size =sizeof(long int);
  char *buffer;
  long int swapped;
  buffer = new char[sizeof(long int)];

  for (int i=0; i<size; i++)
    buffer[ i ] = value[ size-1-i ];
  

  swapped= *( (long int *) buffer );
  delete [] buffer;
  return swapped;
}
//----------------------------
long long int longlongintSwap(char *value)

//---------------------------
{
  int size =sizeof(long long int);
  char *buffer;
 long long int swapped;
  buffer = new char[sizeof(long long int)];

  for (int i=0; i<size; i++)
    buffer[ i ] = value[ size-1-i ];
  

  swapped= *( (long long int *) &buffer );
  delete [] buffer;
  return swapped;
}
//----------------------------

 void sortarray(int *vect, int elements)

//---------------------------
{	
for(int i=elements-1;i>=0;i=i-1) 
{
	for(int j=1;j<=i;j++) 
	{
     	    if(vect[j-1] > vect[j])  
	    {   
		int temp=vect[j];
              	vect[j]=vect[j-1];
	       	vect[j-1]=temp;	      

               }
      
	    }
}
	return;
}

//---------------------------------------------------------------------
/*int debug_cb ( CURL *handle, curl_infotype type,
               char *data, size_t size, void *userp) 
//---------------------------------------------------------------------
{
  switch ( type ) 
  {
    case CURLINFO_TEXT: {
//      printf("INFO: %s\n", data);
      std::cerr<<"INFO: "<<data<<std::endl;
      break;
    }
    case CURLINFO_HEADER_IN: {
//      printf("RESPONSE: %s", data);
      std::cerr<<"RESPONSE: "<<data<<std::endl;
      break;
    }
    case CURLINFO_HEADER_OUT: { /* There is no null-terminator on this one ! 
      size_t i;
//      printf("REQUEST: \n");
      std::cerr<<"REQUEST: "<<std::endl;
      for ( i = 0; i < size; i ++) std::cerr<<data[i] <<std::endl;
      break;
    }
    case CURLINFO_DATA_IN: {
//      printf("RECIEVED: %d bytes\n", size);
      std::cerr<<"RECIEVED: "<<size<<" bytes"<<std::endl;
      break;
    }
    case CURLINFO_DATA_OUT: {
//       printf("TRANSMIT: %d bytes\n", size);
      std::cerr<<"TRANSMIT: "<<size<<" bytes"<<std::endl;
       break;
    }
    case CURLINFO_END: {
//       printf("This should never happen!");
       std::cerr<<"This should never happen!"<<std::endl;
       break;
    }
  }
  return 0;
} */
//---------------------------------------------------------------------
size_t testWriteCurl(void *ptr, size_t size, size_t nmemb,
                                 FILE *stream)
//---------------------------------------------------------------------

{
  fwrite(ptr, size, nmemb, stream);
  return(nmemb*size);
}
//---------------------------------------------------------------------
bool remoteDownloadFiles(std::string remoteFile,std::string login,std::string downloadedFile)
//---------------------------------------------------------------------
//!Metodo per scaricare TUTTI i file remoti in locale
// downloadedFile must be complete of PATH
{
  std::string fileName;
  CURL *curl;  
  CURLcode result;  
  static char errorBuffer[CURL_ERROR_SIZE];  
   FILE * buffer;

#ifdef LIGHT
	if(remoteFile.substr(0,7) =="sftp://")
	{                
		std::cerr << "Invalid sft::// request. VisIVO LIGHT model cannot access file ";
		std::cerr <<remoteFile<<std::endl;
    		return false;

	}
#endif
	if(remoteFile.substr(0,6) =="lfn://")
	{                
#ifndef GLITE
		  std::cerr << "Invalid lfn file ";
		  std::cerr <<remoteFile<<" VisIVO is not compiled with GLITE option."<<std::endl;
		  return false;	
#else
		if(login.empty())
		{
		  std::cerr << "Empty VO Operation in logical filename aborted."<<std::endl;
		  return false;	
		} 
	        int i=0;
		std::string downLoadFile;
		if(downloadedFile.substr(0,1)=="/")
		{
		  downLoadFile="file:/"+downloadedFile;
		}else{
		  std::string tmp;
		  tmp=getenv("PWD");
		  downLoadFile=tmp+"/"+downloadedFile;
		}
		if(downloadedFile.substr(0,7)=="file://") downLoadFile=downloadedFile;
		
	
		
                char p1[512], p2[512],p3[512];
		
                strcpy(p1,remoteFile.c_str());
                strcpy(p2,downLoadFile.c_str());
                strcpy(p3,login.c_str());
//		std::clog<<"lgc_cp: "<<p1<<" "<<p2<<" "<<p3;
                i=lcg_cp(p1, p2, p3, 1, NULL, 0, 0);
		if(i!=0)
		{  
		  std::cerr <<"Invalid copy of lfn file on ";
		  std::cerr <<downloadedFile<<" Operation failed."<<std::endl;
		  return false;
		}
    		return true;
#endif
	}

     curl = curl_easy_init();
     if (curl)  
     {  
      buffer=fopen(downloadedFile.c_str(),"wb");
      
       //! Now set up all of the curl options  
       curl_easy_setopt(curl, CURLOPT_VERBOSE, 1);  
       curl_easy_setopt(curl, CURLOPT_ERRORBUFFER, errorBuffer);  
       curl_easy_setopt(curl, CURLOPT_URL,remoteFile.c_str());  
//       curl_easy_setopt(curl, CURLOPT_DEBUGFUNCTION, debug_cb); //QUI eleiminare?

//       curl_easy_setopt(curl, CURLOPT_URL,"sftp://astrct.oact.inaf.it/home/ube/samp_out_fl_0.5000.bin");  
//       curl_easy_setopt(curl, CURLOPT_URL,"sftp://inaf-node-84.ct.trigrid.it/~/samp_out_fl_0.3000.bin");  
//sftp://user:pass@hostname/       
//COMMENT
	curl_easy_setopt(curl, CURLOPT_HEADER, 0);  
//       curl_easy_setopt(curl, CURLOPT_USERPWD, "user:passwd");
       if(login !="nouser") curl_easy_setopt(curl, CURLOPT_USERPWD, login.c_str());
//       curl_easy_setopt(curl, CURLOPT_HEADER, 1);  
       curl_easy_setopt(curl,CURLOPT_FTP_USE_EPSV, 0);  
       curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1);  
//       curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, writer);  // 
       curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, NULL);  
//       curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, *testWriteCurl);  
       curl_easy_setopt(curl, CURLOPT_WRITEDATA, buffer);  
//       long curlBufferSize=10000;
//       curl_easy_setopt(curl, CURLOPT_BUFFERSIZE,curlBufferSize);//QUI eleiminare?
//	curl_off_t curlMaxFileSizeLarge=96000000000;
//        curl_easy_setopt(curl,CURLOPT_MAXFILESIZE_LARGE,curlMaxFileSizeLarge);
///// Test Portsmouth
	curl_easy_setopt(curl, CURLOPT_NOSIGNAL,1L);
	curl_easy_setopt(curl, CURLOPT_TIMEOUT,45000L);
	curl_easy_setopt(curl, CURLOPT_FTP_RESPONSE_TIMEOUT,44000L);
/////
      //! Attempt to retrieve the remote page  
       result = curl_easy_perform(curl);  
       fclose(buffer);
       if(result != CURLE_OK)  //!check if LOCAL file and copy it!!
       {
       		//! Always cleanup  
       		curl_easy_cleanup(curl);  
                std::cerr << "Error: [" << result << "] - " << errorBuffer;  
   		std::cerr<<" Cannot open or access  file "<<remoteFile<<std::endl;
    		return false;
     	}
       //! Always cleanup  
       curl_easy_cleanup(curl);  

       time_t rawtime;
       time ( &rawtime );

       std::string listFilename;
	listFilename=getDir(downloadedFile);
	listFilename.append("downloaded_VisIVO_filelist");
	std::ofstream fileListOutput(listFilename.c_str(), std::ios::out | std::ios::app);
 	fileListOutput<<"Remote file: "<<remoteFile<<" Local file: "<<downloadedFile<<" Date: "<<ctime(&rawtime)<<std::endl;
 	std::cout<<"Downloaded file "<<remoteFile<<" in local file "<<downloadedFile<<" Date: "<<ctime(&rawtime)<<std::endl;
	fileListOutput.close();
	return true;

    }
  return true;
}
//----------------------------------------------------------------------------
bool remoteLfn(std::string localFile,std::string se,std::string vo,std::string lfnFile,bool vbt)
//----------------------------------------------------------------------------
{
// localFile must be complete of PATH
#ifndef GLITE
		  std::cerr << "Invalid lfn file ";
		  std::cerr <<lfnFile<<" VisIVO is not compiled with GLITE option."<<std::endl;
		  return false;	
#else		
/*		std::clog<<"localFile="<<localFile<<std::endl;
		std::clog<<"se="<<se<<std::endl;
		std::clog<<"vo="<<vo<<std::endl;
		std::clog<<"lfnFile="<<lfnFile<<std::endl;*/
		if(localFile.empty())
		{
		  std::cerr << "Error in local filename ";
		  std::cerr <<localFile<<" Copy in logical filename aborted."<<std::endl;
		  return false;	
		} 
		if(lfnFile.empty() || lfnFile.substr(0,6) != "lfn://")
		{
		  std::cerr << "Error in lfn filename ";
		  std::cerr <<lfnFile<<" Copy in logical filename aborted."<<std::endl;
		  return false;	
		} 
		if(vo.empty())
		{
		  std::cerr << "Error in se or vo options se: "<<se<<" vo:"<<vo;
		  std::cerr<<" Copy in logical filename aborted."<<std::endl;
		  return false;	
		} 
		if(se.empty())
		  se=getenv("DPM_HOST");

		std::string localFileDir=getDir(localFile);
//		std::clog<<"localFileDir="<<localFileDir<<std::endl;
		if(localFileDir.substr(0,2)=="./")
		{ 
		  localFile.erase(0,2);
		  std::string tmp;
		  tmp=getenv("PWD");
		  localFile=tmp+"/"+localFile;
		  localFileDir=getDir(localFile);
		}
		if(localFileDir.empty())
		{
		  std::string tmp;
		  tmp=getenv("PWD");
		  localFile=tmp+"/"+localFile;
		}
		
		std::string localFilename;
		if(localFile.substr(0,1)=="/") localFilename="file:/"+localFile;
		if(localFile.substr(0,7)=="file://") localFilename=localFile;
		
		int i=0;
                char p1[512], p2[512],p3[512],p4[512];

                strcpy(p1,localFilename.c_str());
                strcpy(p2,se.c_str());
                strcpy(p3,lfnFile.c_str());
                strcpy(p4,vo.c_str());
//		std::clog<<"p1="<<p1<<" p2="<<p2<<" p3="<<p3<<" p4="<<p4<<std::endl;
                lcg_del(p3, 1, NULL, p4, NULL, 0, 0);
                i=lcg_cr(p1, p2, NULL, p3, p4, NULL, 1, NULL, 0, 0, NULL);
		
		remove(localFile.c_str());
		if(i!=0)
		{  
		  std::cerr << "Invalid output lfn file ";
		  std::cerr <<lfnFile<<" Operation failed."<<std::endl;
		  return false;
		}
		if(vbt)
		{
		  localFilename=localFilename+".head";
		  lfnFile=lfnFile+".head";
                  strcpy(p1,localFilename.c_str());
                  strcpy(p2,se.c_str());
                  strcpy(p3,lfnFile.c_str());
                  strcpy(p4,vo.c_str());
		  lcg_del(p3, 1, NULL, p4, NULL, 0, 0);
                  i=lcg_cr(p1, p2, NULL, p3, p4, NULL, 1, NULL, 0, 0, NULL);
		  localFile=localFile+".head";
		  remove(localFile.c_str());
		  if(i!=0)
		 {  
		   std::cerr << "Invalid header output lfn file ";
		   std::cerr <<lfnFile<<" Operation failed."<<std::endl;
		   return false;
		  }
		}
    		return true;
#endif  
}
//----------------------------------------------------------------------------
// rgb,hsv in range [0.0-1.0], 
//----------------------------------------------------------------------------
void HSVtoRGB( double *r, double *g, double *b, double h, double s, double v )
{
	int i;
	double f, p, q, t;

	if( s == 0 ) {
		// achromatic (grey)
		*r = *g = *b = v;
		return;
	}

	h *= 6;			// sector 0 to 5
	i = floor( h );
	f = h - i;			// factorial part of h
	p = v * ( 1 - s );
	q = v * ( 1 - s * f );
	t = v * ( 1 - s * ( 1 - f ) );

	switch( i ) {
		case 0:
			*r = v;
			*g = t;
			*b = p;
			break;
		case 1:
			*r = q;
			*g = v;
			*b = p;
			break;
		case 2:
			*r = p;
			*g = v;
			*b = t;
			break;
		case 3:
			*r = p;
			*g = q;
			*b = v;
			break;
		case 4:
			*r = t;
			*g = p;
			*b = v;
			break;
		default:		// case 5:
			*r = v;
			*g = p;
			*b = q;
			break;
	}

}

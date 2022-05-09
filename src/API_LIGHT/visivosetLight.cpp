#include "visivodef.h"
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <iostream>
#include <fstream>
#include <stdio.h>  
#include <vector>  

#include "vuparametersparser.h"
#include "vscreatepath.h"
#include "vscreateslices.h"
#include "vscreategenericslices.h"
#ifdef WIN32
#include <time.h>
#endif
#include "vstable.h"
#include "vsstatisticop.h"

std::vector<std::string> VFextFile;
std::vector<std::string> VIextFile;

// extern "C"{
  
int VI_Import(VisIVOImporter *env);
int VF_Filter(VisIVOFilter *env);

    
    
//----------------------------
int VI_Init(VisIVOImporter *env)
//---------------------------
{
  for(int i=0;i<NPAR;i++)
    env->setatt[i]=0;  
  return noError;
}  
//----------------------------
int VI_Clean(VisIVOImporter *env)
//---------------------------
{
  for(int i=0;i<NPAR;i++)
    env->setatt[i]=0;
  
    for(int i=0;i<VIextFile.size();i++)
      remove(VIextFile[i].c_str());
    
  VIextFile.erase(VIextFile.begin(),VIextFile.end());

  return noError;
}  
//----------------------------
int VI_SetAtt(VisIVOImporter *env, int code, char *value)
//---------------------------
{
std::string sValue(value);
if(code<0 || code>=NPAR)
{
    std::cerr<<"Invalid VI_SetAtt parameter code: "<<code<<std::endl;
    std::cerr<<"SetAtt ignored"<<code<<std::endl;
    return invalidParCode;
} else
  env->setatt[code]=1;
if(code==VI_SET_FFORMAT)
{  
  if(sValue!="ascii" && sValue!="csv" && sValue!="binary" 
	  && sValue!="fly" && sValue!="fitstable"  && sValue!="gadget"
	  && sValue!="hdf5" && sValue!="rawpoints" && sValue!="rawgrids"
	  && sValue!="xml" && sValue!="votable" && sValue!="muportal"
	  && sValue !="ramses")
  {
    std::cerr<<"Invalid value VI_SET_FFORMAT "<<sValue<<std::endl;
    env->setatt[code]=0;
    return invalidParCode;
  } else
    strcpy(env->fformat ,sValue.c_str());
}
if(code==VI_SET_FILEPATH)
	strcpy(env->infile,sValue.c_str());
if(code==VI_SET_OUTFILEVBT)
{
	strcpy(env->outfile,sValue.c_str());
	VIextFile.push_back(sValue);
	VIextFile.push_back(sValue+".head");
}
if(code==VI_SET_VOLUME)
{
	std::stringstream sstmp(sValue);
	sstmp>>env->comp[0];
	sstmp>>env->comp[1];
	sstmp>>env->comp[2];
	sstmp>>env->size[0];
	sstmp>>env->size[1];
	sstmp>>env->size[2];
	if(env->comp[0]<=0 || env->comp[1]<=0 || env->comp[2]<=0 ||
	     env->size[0]<=0.0 || env->size[1]<=0.0 || env->size[2]<=0.0)
	{
	    std::cerr<<"Invalid value VI_SET_VOLUME "<<sValue<<std::endl;
	    env->setatt[code]=0;
	    return invalidParCode;
	} 
}
if(code==VI_SET_USERPWD)
	strcpy(env->userpwd,sValue.c_str());
if(code==VI_SET_BINARYHEADER)
	strcpy(env->binaryheader,sValue.c_str());
if(code==VI_SET_MISSINGVALUE)
{
	std::stringstream sstmp(sValue);
	sstmp>>env->missing;
}
if(code==VI_SET_TEXTVALUE)
{
	std::stringstream sstmp(sValue);
	sstmp>>env->text;
}
if(code==VI_SET_NPOINTS)
{
	std::stringstream sstmp(sValue);
	sstmp>>env->npoints;
}
if(code==VI_SET_DATASETLIST)
	strcpy(env->datasetList,sValue.c_str());
if(code==VI_SET_HYPERSLAB)
	strcpy(env->hyperslab,sValue.c_str());
if(code==VI_SET_VO)
	strcpy(env->VO,sValue.c_str());
if(code==VI_SET_LFNOUT)
	strcpy(env->lfnout,sValue.c_str());
if(code==VI_SET_SE)
	strcpy(env->se,sValue.c_str());

return noError;

} //end VI_SetAtt

//----------------------------
int VF_Init(VisIVOFilter *env)
//---------------------------
{
  for(int i=0;i<NPAR;i++)
    env->setatt[i]=0;  
  strcpy(env->appendList,"NoFile");
  strcpy(env->limitsList,"NoFile");
  strcpy(env->geometryList,"NoFile");
  strcpy(env->mergeList,"NoFile");  
  strcpy(env->visualList,"NoFile");  
  strcpy(env->MRgeometry,"NoFile");  
  for(int i=0;i<VFextFile.size();i++)
      remove(VFextFile[i].c_str());
  VFextFile.erase(VFextFile.begin(),VFextFile.end());
 
  return noError;
}  
//----------------------------
int VF_Clean(VisIVOFilter *env)
//---------------------------
{
  for(int i=0;i<NPAR;i++)
    env->setatt[i]=0;  
  strcpy(env->appendList,"NoFile");
  strcpy(env->limitsList,"NoFile");
  strcpy(env->geometryList,"NoFile");
  strcpy(env->mergeList,"NoFile");
  strcpy(env->visualList,"NoFile");
  strcpy(env->MRgeometry,"NoFile");
  for(int i=0;i<VFextFile.size();i++)
      remove(VFextFile[i].c_str());
  VFextFile.erase(VFextFile.begin(),VFextFile.end());
  
  return noError;
}  
//----------------------------
int VF_SetAtt(VisIVOFilter *env, int code, char *value)
//---------------------------
{
std::string sValue(value);

if(code<VF_PARAM || code>=VF_PARAM+NPAR*2)
{
    std::cerr<<"Invalid VF_SetAtt parameter code: "<<code<<std::endl;
    std::cerr<<"SetAtt ignored"<<code<<std::endl;
    return invalidParCode;
} else
  env->setatt[code-VF_PARAM]=1;

if(code==VF_SET_OPERATION)
    strcpy(env->op ,sValue.c_str());
if(code==VF_SET_FILEVBT)
    strcpy(env->vbt ,sValue.c_str());
if(code==VF_SET_RANDOMPERC)
    strcpy(env->perc ,sValue.c_str());
if(code==VF_SET_RANDOMSEED)
    strcpy(env->seed ,sValue.c_str());
if(code==VF_SET_FIELD)
    strcpy(env->field ,sValue.c_str());
if(code==VF_SET_OUTCOL)
    strcpy(env->outcol ,sValue.c_str());
if(code==VF_SET_OUTVBT)
{
    strcpy(env->out ,sValue.c_str());
    VFextFile.push_back(sValue);
    VFextFile.push_back(sValue+".head");
}
if(code==VF_SET_ADDIDSTART)
    strcpy(env->start ,sValue.c_str());
if(code==VF_SET_APPENDLISTPURGE)
{
  strcpy(env->appendList,"NoFile");
  env->setatt[VF_SET_APPENDLIST-VF_PARAM]=0;
}
if(code==VF_SET_APPENDLIST)
{
    std::stringstream sstmp(sValue);
    std::ofstream appListFile;
    std::string outFile, stCompare;
    stCompare="NoFile";
    outFile=env->appendList;
    if(outFile.compare(stCompare) ==0)
    {
      std::stringstream ssoutFile;
      time_t rawtime;
      struct tm * timeinfo;
      char buffer [80];
      time ( &rawtime );
      timeinfo = localtime ( &rawtime );
      strftime (buffer,80,"%Y%m%d%H%M",timeinfo);
      ssoutFile<<".VS_"<<rand()<<buffer<<"_appendList.txt";  
      outFile=ssoutFile.str();
      strcpy(env->appendList, outFile.c_str());    
      VFextFile.push_back(outFile); 
    } 
    appListFile.open(outFile.c_str(), std::ios::app);
    while(!sstmp.eof())
    {
      std::string stmp;
      sstmp >> stmp;
      appListFile<<stmp<<std::endl;
    }
    appListFile.close();
}
if(code==VF_SET_NEWCOLNAMES)
  strcpy(env->newcolnames,sValue.c_str());
if(code==VF_SET_VOLUMEPERC)
  strcpy(env->volperc,sValue.c_str());
if(code==VF_SET_NEWRES)
  strcpy(env->newres,sValue.c_str());
if(code==VF_SET_LIMITSPURGE)
{
  strcpy(env->limitsList,"NoFile");
  env->setatt[VF_SET_LIMITS-VF_PARAM]=0;
}
if(code==VF_SET_NOAPPEND)
  env->setatt[VF_SET_APPEND-VF_PARAM]=0;

if(code==VF_SET_LIMITS)
{
    std::stringstream sstmp(sValue);
    std::ofstream limitsListFile;
    std::string outFile, stCompare;
    stCompare="NoFile";
    outFile=env->limitsList;
    if(outFile.compare(stCompare) ==0)
    {
      std::stringstream ssoutFile;
      time_t rawtime;
      struct tm * timeinfo;
      char buffer [80];
      time ( &rawtime );
      timeinfo = localtime ( &rawtime );
      strftime (buffer,80,"%Y%m%d%H%M",timeinfo);
      ssoutFile<<".VS_"<<rand()<<buffer<<"_limitsList.txt";  
      outFile=ssoutFile.str();
      strcpy(env->limitsList, outFile.c_str());    
      VFextFile.push_back(outFile); 
    } 
    limitsListFile.open(outFile.c_str(), std::ios::app);
    while(!sstmp.eof())
    {
      std::string stmp;
      sstmp >> stmp;
      limitsListFile<<stmp<<" ";
      sstmp >> stmp;
      limitsListFile<<stmp<<" ";
      sstmp >> stmp;
      limitsListFile<<stmp<<std::endl;
    }
    limitsListFile.close();
}
if(code==VF_SET_CUTTHRESHOLD)
  strcpy(env->threshold,sValue.c_str());
if(code==VF_SET_DECIMATORSKIP)
  strcpy(env->skip,sValue.c_str());
if(code==VF_SET_EXTRACTIONGEOMETRYPURGE)
{
  strcpy(env->geometryList,"NoFile");  
  env->setatt[VF_SET_EXTRACTIONGEOMETRY-VF_PARAM]=0;
}
if(code==VF_SET_EXTRACTIONGEOMETRY)
{
    std::stringstream sstmp(sValue);
    std::ofstream geometryListFile;
    std::string outFile, stCompare;
    stCompare="NoFile";
    outFile=env->geometryList;
    if(outFile.compare(stCompare) ==0)
    {
      std::stringstream ssoutFile;
      time_t rawtime;
      struct tm * timeinfo;
      char buffer [80];
      time ( &rawtime );
      timeinfo = localtime ( &rawtime );
      strftime (buffer,80,"%Y%m%d%H%M",timeinfo);
      ssoutFile<<".VS_"<<rand()<<buffer<<"_geometry.txt";  
      outFile=ssoutFile.str();
      strcpy(env->geometryList, outFile.c_str());    
      VFextFile.push_back(outFile); 
    } 
    geometryListFile.open(outFile.c_str(), std::ios::app);
    while(!sstmp.eof())
    {
      std::string stmp;
      sstmp >> stmp;
      geometryListFile<<stmp<<" ";
      sstmp >> stmp;
      geometryListFile<<stmp<<std::endl;
    }
    geometryListFile.close();
}

if(code==VF_SET_MRGEOMETRYPURGE)
{
  strcpy(env->MRgeometry,"NoFile");  
  env->setatt[VF_SET_MRGEOMETRY-VF_PARAM]=0;
}
if(code==VF_SET_MRGEOMETRY)
{
    std::stringstream sstmp(sValue);
    std::ofstream geometryListFile;
    std::string outFile, stCompare;
    stCompare="NoFile";
    outFile=env->MRgeometry;
    if(outFile.compare(stCompare) ==0)
    {
      std::stringstream ssoutFile;
      time_t rawtime;
      struct tm * timeinfo;
      char buffer [80];
      time ( &rawtime );
      timeinfo = localtime ( &rawtime );
      strftime (buffer,80,"%Y%m%d%H%M",timeinfo);
      ssoutFile<<".VS_"<<rand()<<buffer<<"_mrgeometry.txt";  
      outFile=ssoutFile.str();
      strcpy(env->MRgeometry, outFile.c_str());    
      VFextFile.push_back(outFile); 
    } 
    geometryListFile.open(outFile.c_str(), std::ios::app);
    while(!sstmp.eof())
    {
      std::string stmp;
      sstmp >> stmp;
      geometryListFile<<stmp<<" ";
      sstmp >> stmp;
      geometryListFile<<stmp<<std::endl;
    }
    geometryListFile.close();
}


if(code==VF_SET_STARTINGCELL)
  strcpy(env->startingcell,sValue.c_str());
if(code==VF_SET_RESOLUTION)
  strcpy(env->resolution,sValue.c_str());
if(code==VF_SET_POINTCOLUMNS)
  strcpy(env->pointcolumns,sValue.c_str());
if(code==VF_SET_GRIDORIGIN)
  strcpy(env->gridorigin,sValue.c_str());
if(code==VF_SET_GRIDSPACING)
  strcpy(env->gridspacing,sValue.c_str());
if(code==VF_SET_VOLUME)
  strcpy(env->volume,sValue.c_str());
if(code==VF_SET_BOX)
  strcpy(env->box,sValue.c_str());
if(code==VF_SET_NOPERIODIC)
  env->setatt[VF_SET_PERIODIC-VF_PARAM]=0;
if(code==VF_SET_NUMBIN)
  strcpy(env->numbin,sValue.c_str());
if(code==VF_SET_INTERVAL)
{  
  std::stringstream sstmp(sValue);
  std::string stmp;
  sstmp>>stmp;
  strcpy(env->intfrom,stmp.c_str());
  sstmp>>stmp;
  strcpy(env->intto,stmp.c_str());
}
if(code==VF_SET_INFILES)
{  
  std::stringstream sstmp(sValue);
  std::string stmp;
  sstmp>>stmp;
  strcpy(env->infiles1,stmp.c_str());
  sstmp>>stmp;
  strcpy(env->infiles2,stmp.c_str());
}
if(code==VF_SET_CENTER)
  strcpy(env->center,sValue.c_str());
if(code==VF_SET_RADIUS)
  strcpy(env->radius,sValue.c_str());
if(code==VF_SET_INVALUE)
  strcpy(env->invalue,sValue.c_str());
if(code==VF_SET_OUTVALUE)
  strcpy(env->outvalue,sValue.c_str());
if(code==VF_SET_MATEXPRESSION)
  strcpy(env->math,sValue.c_str());
if(code==VF_SET_MERGEPAD)
  strcpy(env->pad,sValue.c_str());
if(code==VF_SET_MERGELISTPURGE)
{  
  strcpy(env->mergeList,"NoFile");
  env->setatt[VF_SET_MERGELIST-VF_PARAM]=0;
}
if(code==VF_SET_MERGELIST)
{
    std::stringstream sstmp(sValue);
    std::ofstream mergeListFile;
    std::string outFile, stCompare;
    stCompare="NoFile";
    outFile=env->mergeList;
    if(outFile.compare(stCompare) ==0)
    {
      std::stringstream ssoutFile;
      time_t rawtime;
      struct tm * timeinfo;
      char buffer [80];
      time ( &rawtime );
      timeinfo = localtime ( &rawtime );
      strftime (buffer,80,"%Y%m%d%H%M",timeinfo);
      ssoutFile<<".VS_"<<rand()<<buffer<<"_merge.txt";  
      outFile=ssoutFile.str();
      strcpy(env->mergeList, outFile.c_str());    
      VFextFile.push_back(outFile); 
    } 
    mergeListFile.open(outFile.c_str(), std::ios::app);
    while(!sstmp.eof())
    {
      std::string stmp;
      sstmp >> stmp;
      mergeListFile<<stmp<<" ";
      sstmp >> stmp;
      mergeListFile<<stmp<<std::endl;
    }
    mergeListFile.close();
}
if(code==VF_SET_NUMROWS)
  strcpy(env->numrows,sValue.c_str());
if(code==VF_SET_RANGEROWS)
{  
  std::stringstream sstmp(sValue);
  std::string stmp;
  sstmp>>stmp;
  strcpy(env->rangerowsfrom,stmp.c_str());
  sstmp>>stmp;
  strcpy(env->rangerowsto,stmp.c_str());
}
if(code==VF_SET_WIDTH)
  strcpy(env->width,sValue.c_str());
if(code==VF_SET_PRECISION)
  strcpy(env->precision,sValue.c_str());
if(code==VF_SET_OUT)
  strcpy(env->outtxt,sValue.c_str());
if(code==VF_SET_NUMCELLS)
  strcpy(env->numcells,sValue.c_str());
if(code==VF_SET_NSIGMA)
  strcpy(env->nsigma,sValue.c_str());
if(code==VF_SET_VOLUMESPLIT)
  strcpy(env->volumesplit,sValue.c_str());
if(code==VF_SET_NUMOFTABLES)
  strcpy(env->numoftables,sValue.c_str());
if(code==VF_SET_MAXSIZETABLE)
  strcpy(env->maxsizetable,sValue.c_str());
if(code==VF_SET_VISUALSIZE)
  strcpy(env->visualsize,sValue.c_str());

if(code==VF_SET_VISUALLISTPURGE)
{  
  strcpy(env->visualList,"NoFile");
  env->setatt[VF_SET_VISUALLIST-VF_PARAM]=0;
}
if(code==VF_SET_VISUALLIST)
{
    std::stringstream sstmp(sValue);
    std::ofstream visualListFile;
    std::string outFile, stCompare;
    stCompare="NoFile";
    outFile=env->visualList;
    if(outFile.compare(stCompare) ==0)
    {
      std::stringstream ssoutFile;
      time_t rawtime;
      struct tm * timeinfo;
      char buffer [80];
      time ( &rawtime );
      timeinfo = localtime ( &rawtime );
      strftime (buffer,80,"%Y%m%d%H%M",timeinfo);
      ssoutFile<<".VS_"<<rand()<<buffer<<"_visual.txt";  
      outFile=ssoutFile.str();
      strcpy(env->visualList, outFile.c_str());    
      VFextFile.push_back(outFile); 
    } 
    visualListFile.open(outFile.c_str(), std::ios::app);
    while(!sstmp.eof())
    {
      std::string stmp;
      sstmp >> stmp;
      visualListFile<<stmp<<" ";
      sstmp >> stmp;
      visualListFile<<stmp<<std::endl;
    }
    visualListFile.close();
}

if(code==VF_SET_HISTOGRAMBIN)
{
  if(sValue.empty())
    strcpy(env->histobin,"-1");
  else
    strcpy(env->histobin,sValue.c_str());
}
if(code==VF_SET_LFNOUT)
  strcpy(env->lfnout,sValue.c_str());
if(code==VF_SET_VO)
  strcpy(env->vo,sValue.c_str());
if(code==VF_SET_SE)
  strcpy(env->se,sValue.c_str());
if(code==VF_SET_MRPOS)
  strcpy(env->MRpos,sValue.c_str());
if(code==VF_SET_MRBACKGROUND)
  strcpy(env->MRbackground,sValue.c_str());
if(code==VF_SET_DIMVOX)
  strcpy(env->dimvox,sValue.c_str());
if(code==VF_SET_TRACKPLANEDIST)
  strcpy(env->trackplanedist,sValue.c_str());
if(code==VF_SET_INNERDIST)
  strcpy(env->innerdist,sValue.c_str());
if(code==VF_SET_OUTPOINTS)
  strcpy(env->outpoints,sValue.c_str());
if(code==VF_SET_OUTVOL)
  strcpy(env->outvol,sValue.c_str());


return noError;

}
//
//----------------------------
int VS_VBTMetaData(VBT *tab, int satistic ,char *value)
//---------------------------
{
std::string filename(value);
if(filename.find(".bin") == std::string::npos)
	    filename.append(".bin");
VSTable table(filename);
if(!table.tableExist())
  {
    std::cerr <<"No valid input file table is provided"<<std::endl;
    return invalidVBT;
   }
strcpy(tab->datatype, table.getType().c_str()); 
strcpy(tab->locator, filename.c_str()); 
tab->nOfFields=table.getNumberOfColumns();
tab->nOfRows=table.getNumberOfRows();
tab->nCells[0]=table.getCellNumber()[0];
tab->nCells[1]=table.getCellNumber()[1];
tab->nCells[2]=table.getCellNumber()[2];
tab->cellSize[0]=table.getCellSize()[0];
tab->cellSize[1]=table.getCellSize()[1];
tab->cellSize[2]=table.getCellSize()[2];
strcpy(tab->endianity, table.getEndiannes().c_str());

tab->field= new char*[tab->nOfFields];
int maxSize=MAXCOLNAMESIZE;
for(unsigned int i=0;i<tab->nOfFields;i++)
{
  tab->field[i]= new char[maxSize];
  strcpy(tab->field[i],table.getColName(i).c_str());
}


if(satistic>0)
{
  VSStatisticOp op;
  bool ret;
  ret=op.addInput(&table);
  if(!ret) return invalidVBT;
  op.addParameter("silent","");
  ret=op.execute();
  if(!ret) return invalidVBTMetadata;
  
  tab->statistic=new float*[tab->nOfFields];
  for(unsigned int i=0;i<tab->nOfFields;i++)
  {  
      tab->statistic[i]=new float[4];
      op.getRange(i,tab->statistic[i][0],tab->statistic[i][1],
		  tab->statistic[i][2],tab->statistic[i][3]);
  }
}

return noError;
}
//----------------------------
int VS_VBTPointers(VBT *tab, unsigned int *colList,unsigned  int nOfCols ,
		   long unsigned int fromRow, long unsigned int toRow, float **fArray)
//---------------------------
{
std::string filename(tab->locator);
VSTable table(filename);
if(!table.tableExist())
  {
    std::cerr <<"No valid input file table is provided"<<std::endl;
    return invalidVBT;
   }
int ret=table.getColumn(colList,nOfCols,fromRow,toRow,fArray);
if(ret<0) return invalidPointers;
return noError;
}
//----------------------------
int VS_VBTAllColumn(VBT *tab, int idCol ,float* oneCol) 
//----------------------------
{
std::string filename(tab->locator);
VSTable table(filename);
if(!table.tableExist())
  {
    std::cerr <<"No valid input file table is provided"<<std::endl;
    return invalidVBT;
   }
int ret=table.getColumn(idCol,oneCol);
if(ret<0) return invalidPointers;
return noError;

}
//----------------------------
int VS_VBTPartialColumn(VBT *tab, int idCol ,float* oneCol,int from, int to) 
//----------------------------
{
std::string filename(tab->locator);
VSTable table(filename);
if(!table.tableExist())
  {
    std::cerr <<"No valid input file table is provided"<<std::endl;
    return invalidVBT;
   }
int ret=table.getColumn(idCol,oneCol,from,to);
if(ret<0) return invalidPointers;
return noError;

}
//****************************
//****************************
//****************************

//------------------------Generic
int VS_DirectImage(VisIVOViewer *VVenv,VisIVOImporter *VIenv,float random)
//----------------------------
{
      char tmp[256];
      std::string ImpOutFile,FilOutFile;
      std::stringstream sstmp;
      time_t rawtime;
      struct tm * timeinfo;
      char buffer [80];
      time ( &rawtime );
      timeinfo = localtime ( &rawtime );
      strftime (buffer,80,"%Y%m%d%H%M",timeinfo);
      sstmp<<".VS_"<<rand()<<buffer;
      ImpOutFile=sstmp.str()+"_ImpDirectImage.bin";  
      FilOutFile=ImpOutFile;  
      strcpy(VIenv->outfile,ImpOutFile.c_str());
    
      int errorCode=VI_SetAtt(VIenv,VI_SET_OUTFILEVBT,VIenv->outfile);
      errorCode=VI_Import(VIenv);
      if(random <100.0)
      {
          std::stringstream sstmp1;
	  sstmp1<<random;      
	  std::string stmp;
	  stmp=sstmp1.str();
	  VisIVOFilter DIFilStr;
	  strcpy(tmp,"randomizer");
	  errorCode=VF_SetAtt(&DIFilStr,VF_SET_OPERATION,tmp);
	  strcpy(tmp,stmp.c_str());
	  errorCode=VF_SetAtt(&DIFilStr,VF_SET_RANDOMPERC,tmp);
          FilOutFile=sstmp.str()+"_FilDirectImage.bin";  
	  
	  strcpy(tmp,FilOutFile.c_str());
	  errorCode=VF_SetAtt(&DIFilStr,VF_SET_OUT,tmp);
	  strcpy(tmp,ImpOutFile.c_str());
	  errorCode=VF_SetAtt(&DIFilStr,VF_SET_FILEVBT,tmp);
          FilOutFile=sstmp.str()+"_FilDirectImage.bin";
	  errorCode=VF_Filter(&DIFilStr);
     }
     strcpy(tmp,FilOutFile.c_str());
//     errorCode=VV_SetAtt(VVenv,VV_SET_FILEVBT,tmp);
//     errorCode=VV_View(VVenv);
     remove(ImpOutFile.c_str());
     ImpOutFile+=".head";
     remove(ImpOutFile.c_str());
     remove(FilOutFile.c_str());
     FilOutFile+=".head";
     remove(FilOutFile.c_str());
     
    return noError;
  
}


//---------------------------
void *VS_DirectImage_Thread(void *id)
//---------------------------
{
	int iret;
	VisIVOAsynchId *idTh;
	idTh=(VisIVOAsynchId *)id;
	iret=VS_DirectImage((VisIVOViewer *)(idTh->env), (VisIVOImporter *)(idTh->envSecond), idTh->random);
	idTh->errorCode=iret;
	idTh->state=successfulEndThread;
	return 0;
}

//---------------------------
int createProcessDirectImage(VisIVOAsynchId *id)
//---------------------------
{
	int iret;
	#ifndef WIN32
	pid_t pid;
	id->state=runningThread;
	pid=fork();
	if(pid==-1)
	{
		id->state=errorThread;
		return 0;
	}
	else
	{
		if(pid!=0)
		{
			id->pid=pid;
			return noError;
		}
		else
		{
			// Child
			iret=VS_DirectImage((VisIVOViewer *)(id->env), (VisIVOImporter *)(id->envSecond), id->random);
			exit(iret);
		}
	}
	#endif
	return 0;	
}


//---------------------------
void *VS_DirectImage_Process(void *id)
//---------------------------
{
	int iret, i;
	VisIVOAsynchId *idTh;
	
	#ifndef WIN32
	pid_t pid;
	idTh=(VisIVOAsynchId *)id;
	// New process.
	pid=fork();
	if(pid==-1)
	{
		idTh->state=errorThread;
		return 0;
	}
	else
	{
		if(pid!=0)
		{
			// Parent
			i=waitpid(pid, &iret, 0);
			idTh->state=successfulEndThread;
			if(WIFEXITED(iret))	idTh->errorCode=WEXITSTATUS(iret);
		}
		else
		{
			// Child
			iret=VS_DirectImage((VisIVOViewer *)(idTh->env), (VisIVOImporter *)(idTh->envSecond), idTh->random);
			exit(iret);
		}
	}
	#endif
	return 0;	
}



//------------------------Direct thread
int VA_DirectImage(VisIVOViewer *VVenv, VisIVOImporter *VIenv, float random, VisIVOAsynchId *id)
//------------------------
{
	int iret;
	#ifdef WIN32
 		iret=VS_DirectImage(VVenv, VIenv, random);
 		id->errorCode=iret;
 		return iret;
	#endif
#ifdef MAC
        id->withThread=1;
#endif
	#ifndef WIN32
	id->env=VVenv;
	id->envSecond=VIenv;
	id->random=random;
	id->state=runningThread;
	if(id->withThread>0)	iret=pthread_create(&(id->threadId), NULL, &VS_DirectImage_Thread, (void *)id);
	else	iret=createProcessDirectImage(id);//pthread_create(&(id->threadId), NULL, &VS_DirectImage_Process, (void *)id);
	if(iret!=0)	
	{
		id->state=undefined;
		iret=errorThread;
	}
	return iret;
	#endif
}

//*********************
//*********************
//*********************
//*********************
//*********************

//----------------------------
void VA_Init(VisIVOAsynchId *id)
//----------------------------
{
	id->threadId=NULL;
	id->state=undefined;
	id->errorCode=undefined;
	id->withThread=0;
	id->random=100;
	id->env=NULL;
	id->envSecond=NULL;
}

/*
 * Aspetta la terminazione del thread identificato da id. Ritorna uno dei seguenti valori:
 * 
 * 
 */
//----------------------------
int VA_Wait(VisIVOAsynchId *id)
//----------------------------
{
	int err=undefined;
	int s=id->state;
	#ifndef WIN32
	pid_t p;
	if(id->withThread>0)	// is thread
	{
		if((id->state==runningThread) && (id->threadId!=NULL))	err=pthread_join(id->threadId, NULL);
		if(err!=0)	err=errorThread;
		else   err=noError;
		if(s==successfulEndThread)	return noError; // Thread terminated
	}
	else    // Process	
	{	
		if(id->state==runningThread)
		{
			p=waitpid(id->pid, &err, 0);
			id->state=successfulEndThread;
			if(WIFEXITED(err))	id->errorCode=WEXITSTATUS(err);
		}
	}
	#endif
	return err;
}

//----------------------------
int VA_GetError(VisIVOAsynchId *id)
//----------------------------
{
	int err;
	err=id->errorCode;
	if(id->state==runningThread)	err=undefined;
	return err;
}

//----------------------------
int VA_GetState(VisIVOAsynchId *id)
//----------------------------
{
	int state=-1;
	state=id->state;
	return state;
}

void VA_SetMultiThread(VisIVOAsynchId *id)
{
	id->withThread=1;
}
void VA_SetMultiProc(VisIVOAsynchId *id)
{
	id->withThread=0;
}

//
// } //extern "C"

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
std::vector<std::string> VVextFile;
extern "C"
{
  
int VI_Import(VisIVOImporter *env);
int VF_Filter(VisIVOFilter *env);
int VV_View(VisIVOViewer *env);


//----------------------------
int VV_Init(VisIVOViewer *env)
//---------------------------
{
  for(int i=0;i<NPAR;i++)
    env->setatt[i]=0;  
  
  env->nRows=0;
  env->nCols=0;
  env->cycleOffset=0;
  std::string cycleFile;
  std::stringstream sscycleFile;
  time_t rawtime;
  struct tm * timeinfo;
  char buffer [80];
  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  strftime (buffer,80,"%Y%m%d%H%M%S",timeinfo);
  sscycleFile<<".VS_"<<rand()<<buffer<<"_cycle.par";  
  cycleFile=sscycleFile.str();
  strcpy(env->cycleFile, cycleFile.c_str());    
  return noError;

}  
//----------------------------
int VV_Clean(VisIVOViewer *env)
//---------------------------
{
  for(int i=0;i<NPAR;i++)
    env->setatt[i]=0;
  for(int i=0;i<VVextFile.size();i++)
      remove(VVextFile[i].c_str());
  VVextFile.erase(VVextFile.begin(),VVextFile.end());
  return noError;

}  
//----------------------------
int VV_SetXYZ(VisIVOViewer *env,float *X,float *Y, float *Z,char *names,int nRows)
//---------------------------
{
   std::string snames(names);
   std::stringstream sstmp(snames);
   std::string stmp;

  env->setatt[VV_SET_FILEVBT-VV_PARAM]=0;
  env->setatt[VV_SET_INTERNAL-VV_PARAM]=1;
  env->setatt[VV_SET_INTERNALFIELD-VV_PARAM]=1;
  env->nRows=nRows;
  env->nCols+=3;
  sstmp >> stmp;
  strcpy(env->xField,stmp.c_str());//set X axis name
  strcpy(env->internalNames[env->nCols-3],stmp.c_str());//set X axis name
  sstmp >> stmp;
  strcpy(env->yField,stmp.c_str());//set Y axis name
  strcpy(env->internalNames[env->nCols-2],stmp.c_str());//set Y axis name
  sstmp >> stmp;
  strcpy(env->zField,stmp.c_str());//set Z axis name
  strcpy(env->internalNames[env->nCols-1],stmp.c_str());//set Z axis name
  env->internalPointers[env->nCols-3]=X;
  env->internalPointers[env->nCols-2]=Y;
  env->internalPointers[env->nCols-1]=Z;
  return noError;

}
//----------------------------
int VV_SetVect(VisIVOViewer *env,float *VX,float *VY, float *VZ,char *names,int nRows)
//---------------------------
{
   std::string snames(names);
   std::stringstream sstmp(snames);
   std::string stmp;
  env->setatt[VV_SET_INTERNAL-VV_PARAM]=1;
  env->setatt[VV_SET_INTERNALVFIELD-VV_PARAM]=1;
  env->setatt[VV_SET_VECTOR-VV_PARAM]=1;
//  std::clog<<"TEST: "<<env->setatt[65]<<std::endl;
  if(env->nRows==0) env->nRows=nRows;
  env->nCols+=3;
  sstmp >> stmp;
  strcpy(env->xVectorField,stmp.c_str());//set VX name
  strcpy(env->internalNames[env->nCols-3],stmp.c_str());//set VX name
  sstmp >> stmp;
  strcpy(env->yVectorField,stmp.c_str());//set VY name
  strcpy(env->internalNames[env->nCols-2],stmp.c_str());//set VY name
  sstmp >> stmp;
  strcpy(env->zVectorField,stmp.c_str());//set VZ  name
  strcpy(env->internalNames[env->nCols-1],stmp.c_str());//set VZ  name
  env->internalPointers[env->nCols-3]=VX;
  env->internalPointers[env->nCols-2]=VY;
  env->internalPointers[env->nCols-1]=VZ;
  return noError;
}
//----------------------------
int VV_SetColorScalar(VisIVOViewer *env,float *color, char *name, int nRows)
//---------------------------
{
  env->setatt[VV_SET_INTERNAL-VV_PARAM]=1;
  env->setatt[VV_SET_INTERNALCOLORSCALAR-VV_PARAM]=1;
  if(env->nRows==0) env->nRows=nRows;
  env->nCols+=1;
  strcpy(env->colorScalar,name);//set
  strcpy(env->internalNames[env->nCols-1],name);//set 
  env->internalPointers[env->nCols-1]=color;
  return noError;
}
//----------------------------
int VV_SetRadiusScalar(VisIVOViewer *env,float *radiusscalar, char *name, int nRows)
//---------------------------
{
  env->setatt[VV_SET_INTERNAL-VV_PARAM]=1;
  env->setatt[VV_SET_INTERNALRADIUSSCALAR-VV_PARAM]=1;
  if(env->nRows==0) env->nRows=nRows;
  env->nCols+=1;
  strcpy(env->radiusscalar,name);//set  
  strcpy(env->internalNames[env->nCols-1],name);//set 
  env->internalPointers[env->nCols-1]=radiusscalar;
  return noError;
}
//----------------------------
int VV_SetHeightScalar(VisIVOViewer *env,float *heightscalar, char *name, int nRows)
//---------------------------
{
  env->setatt[VV_SET_INTERNAL-VV_PARAM]=1;
  env->setatt[VV_SET_INTERNALHEIGHTSCALAR-VV_PARAM]=1;
  if(env->nRows==0) env->nRows=nRows;
  env->nCols+=1;
  strcpy(env->heightscalar,name);//set  
  strcpy(env->internalNames[env->nCols-1],name);//set 
  env->internalPointers[env->nCols-1]=heightscalar;
  return noError;
}
//----------------------------
int VV_SetVolume(VisIVOViewer *env,float *vol,char *name,char *geo,char *size)
//---------------------------
{
  env->setatt[VV_SET_FILEVBT-VV_PARAM]=0;
  env->setatt[VV_SET_INTERNAL-VV_PARAM]=1;
  env->setatt[VV_SET_INTERNALVRENDERFIELD-VV_PARAM]=1;
  env->nCols+=1;
  strcpy(env->vRenderingField,name);//set vol name
  strcpy(env->isosurfaceField,name);//set vol name
  strcpy(env->slicefield,name);//set vol name
  strcpy(env->internalNames[env->nCols-1],name);//set Volume  Field name
  env->internalPointers[env->nCols-1]=vol;
  env->internalPointers[env->nCols-1]=vol;
  std::string sComp(geo); 
  std::stringstream ssComp(sComp);
  ssComp>>env->comp[0];
  ssComp>>env->comp[1];
  ssComp>>env->comp[2];
  env->nRows=env->comp[0]*env->comp[1]*env->comp[2];
  std::string sSize(size); 
  std::stringstream ssSize(sSize);
  ssSize>>env->size[0];
  ssSize>>env->size[1];
  ssSize>>env->size[2];
  return noError;
}
//----------------------------
int VV_SetCameraPath(VisIVOViewer *env,int type,float *camera,int zoomend,float *zsf,int framesec, 
		     int length, float *campos, float *camfp, float *camroll,int fcycle)
//---------------------------
{
  
    if(type!=0 && type!=1 && type!=2 && type!=3) 
    {
	std::cerr<<"Invalid type in CameraPath function. Camera path Ignored"<<std::endl;
	return invalidParCode;
    }
    if(fcycle==2)
    {
      env->setatt[VV_SET_CYCLE-VV_PARAM]=0;
      return 0;
    }
    env->setatt[VV_SET_CYCLE-VV_PARAM]=1;
    VSCreatePathUT op;
    
    std::stringstream commandParametersSStream;
    std::map<std::string, std::string> appParameters;
    std::map<std::string, std::string>::iterator iter;
    
    commandParametersSStream<<" --op createpath";
    commandParametersSStream<<" --type "<<type;
    if(type!=2) commandParametersSStream<<" --azimuth "<<camera[0]<<" "<<camera[1];
    if(type!=2) commandParametersSStream<<" --elevation "<<camera[2]<<" "<<camera[3];
    commandParametersSStream<<" --zoom "<<camera[4]<<" "<<camera[5];
    switch(type)
    {
      case 0:{
	if(camroll[0]!=NOSET_CAM && camroll[1]!=NOSET_CAM)
	  commandParametersSStream<<" --camroll "<<camroll[0]<<" "<<camroll[1];
	break;
      }
      
      case 1:{
	bool fpValid=true;
	for(int i=0;i<6;i++)
	  if(camfp[i]==NOSET_CAM) fpValid=false;
	  
	if(fpValid) 
	{
	  commandParametersSStream<<" --camfp";
	  for(int i=0;i<6;i++)
	    commandParametersSStream<<" "<<camfp[i];
	}
	if(camroll[0]!=NOSET_CAM && camroll[1]!=NOSET_CAM)
	  commandParametersSStream<<" --camroll "<<camroll[0]<<" "<<camroll[1];	
	break;
      }

      case 2:
      case 3:{
	bool fpValid=true;
	bool posValid=true;
	for(int i=0;i<6;i++)
	  if(camfp[i]==NOSET_CAM) fpValid=false;
	for(int i=0;i<6;i++)
	  if(campos[i]==NOSET_CAM) posValid=false;
	  
	if(posValid) 
	{
	  commandParametersSStream<<" --campos";
	  for(int i=0;i<6;i++)
	    commandParametersSStream<<" "<<campos[i];
	}
	if(fpValid) 
	{
	  commandParametersSStream<<" --camfp";
	  for(int i=0;i<6;i++)
	    commandParametersSStream<<" "<<camfp[i];
	}
	if(camroll[0]!=NOSET_CAM && camroll[1]!=NOSET_CAM)
	  commandParametersSStream<<" --camroll "<<camroll[0]<<" "<<camroll[1];	
	break;
      }
      
    }
    if(zoomend==1)
    { 
      commandParametersSStream<<" --zoomend ";
      if(*zsf>0.0) commandParametersSStream<<*zsf;
    }  
    if(framesec>0) commandParametersSStream<<" --framesec "<<framesec;
    if(length>0) commandParametersSStream<<" --length "<<length;
    
    if(fcycle==0)
    {  
      std::string outFile;
      std::stringstream ssoutFile;
      time_t rawtime;
      struct tm * timeinfo;
      char buffer [80];
      time ( &rawtime );
      timeinfo = localtime ( &rawtime );
      strftime (buffer,80,"%Y%m%d%H%M%S",timeinfo);
      ssoutFile<<".VS_"<<rand()<<buffer<<"_cycle.par";  
      outFile=ssoutFile.str();
      strcpy(env->cycleFile, outFile.c_str());    
      commandParametersSStream<<" --out "<<outFile;
      VVextFile.push_back(outFile);
    }
    if(fcycle==1)
    {  
      std::string outFile(env->cycleFile);
      commandParametersSStream<<" --out "<<outFile;
      VVextFile.push_back(outFile);
    }
    
//    std::clog<<"Cycle Utility "<<commandParametersSStream.str();

    VUParametersParser myparser(commandParametersSStream.str());

    appParameters=myparser.getParameters();
    
    op.setParameters(appParameters);
    bool ret=op.execute();
    if(!ret) return invalidPathCreation;
    return noError;
}
//----------------------------
int VV_SetGSliceScan(VisIVOViewer *env,float *fgslice, int *iglsice)
// Generic Slice scan
//---------------------------
{
    env->setatt[VV_SET_CYCLE-VV_PARAM]=1;
    env->setatt[VV_SET_ISOSURFACE-VV_PARAM]=0;
    VSCreateGenericSlicesUT op;
    
    std::stringstream commandParametersSStream;
    std::map<std::string, std::string> appParameters;
    std::map<std::string, std::string>::iterator iter;
    
    commandParametersSStream<<" --op genericslices";
    commandParametersSStream<<" --point "<<fgslice[0]<<" "<<fgslice[1]<<" "<<fgslice[2];
    commandParametersSStream<<" --normal "<<fgslice[3]<<" "<<fgslice[4]<<" "<<fgslice[5];
    commandParametersSStream<<" --size "<<fgslice[6];
    commandParametersSStream<<" --step "<<iglsice[0];
    if(iglsice[1]==1)  commandParametersSStream<<" --movedown";
    
    std::string outFile;
    std::stringstream ssoutFile;
    time_t rawtime;
    struct tm * timeinfo;
    char buffer [80];
    time ( &rawtime );
    timeinfo = localtime ( &rawtime );
    strftime (buffer,80,"%Y%m%d%H%M",timeinfo);
    ssoutFile<<".VS_"<<rand()<<buffer<<"_cycle.par";  
    outFile=ssoutFile.str();
    strcpy(env->cycleFile, outFile.c_str());    
    commandParametersSStream<<" --out "<<outFile;
    VVextFile.push_back(outFile);
    
//    std::clog<<"Cycle Utility "<<commandParametersSStream.str()<<std::endl;

    VUParametersParser myparser(commandParametersSStream.str());

    appParameters=myparser.getParameters();
    
    op.setParameters(appParameters);
    bool ret=op.execute();
    if(!ret) return invalidPathCreation;
    return noError;

}
//----------------------------
int VV_EnableSplotch(VisIVOViewer *env)
//---------------------------
{    
    env->setatt[VV_SET_SPLOTCH-VV_PARAM]=1;
    env->setatt[VV_SET_VOLUME-VV_PARAM]=0; //disable volume
    std::string outFile;
    std::stringstream ssoutFile;
    time_t rawtime;
    struct tm * timeinfo;
    char buffer [80];
    time ( &rawtime );
    timeinfo = localtime ( &rawtime );
    strftime (buffer,80,"%Y%m%d%H%M",timeinfo);
    ssoutFile<<".VS_"<<rand()<<buffer<<"_splotch.par";  
    outFile=ssoutFile.str();
    strcpy(env->splotchpar, outFile.c_str());    
    VVextFile.push_back(outFile);
  
  return noError;

}
//----------------------------
int VV_SetAtt(VisIVOViewer *env, int code, char *value)
//---------------------------
{
std::string sValue(value);

if(code<VV_PARAM || code>=VV_PARAM+NPAR)
{
    std::cerr<<"Invalid VV_SetAtt parameter code: "<<code<<std::endl;
    std::cerr<<"SetAtt ignored"<<code<<std::endl;
    return invalidParCode;
} else
  env->setatt[code-VV_PARAM]=1;

if(code==VV_SET_EXTERNAL)
{
  env->setatt[VV_SET_INTERNAL-VV_PARAM]=0;
  env->setatt[VV_SET_INTERNALFIELD-VV_PARAM]=0;
  env->setatt[VV_SET_INTERNALVFIELD-VV_PARAM]=0;
  env->setatt[VV_SET_INTERNALCOLORSCALAR-VV_PARAM]=0;
  env->setatt[VV_SET_INTERNALRADIUSSCALAR-VV_PARAM]=0;
  env->setatt[VV_SET_INTERNALHEIGHTSCALAR-VV_PARAM]=0;
  env->setatt[VV_SET_INTERNALVRENDERFIELD-VV_PARAM]=0; 
}

if(code==VV_SET_FILEVBT)
{
	if(env->setatt[VV_SET_INTERNAL-VV_PARAM]==1)
	{
	  std::cerr<<"SET_INTERNAL is present. SET_PATH ignored"<<std::endl;
	  env->setatt[code-VV_PARAM]=0;
	} else
          strcpy(env->path, sValue.c_str());
}
if(code==VV_SET_FIELD)
{
	std::stringstream sstmp(sValue);
	std::string stmp;

	sstmp >> stmp;
        strcpy(env->xField, stmp.c_str());
	sstmp >> stmp;
        strcpy(env->yField, stmp.c_str());
	sstmp >> stmp;
        strcpy(env->zField, stmp.c_str());
}
if(code==VV_SET_COLORSCALAR)
	strcpy(env->colorScalar,sValue.c_str());

if(code==VV_SET_VRENDERFIELD)
{
	strcpy(env->vRenderingField,sValue.c_str());
	strcpy(env->colorScalar,sValue.c_str());
}

if(code==VV_SET_ISOSURFIELD)
	strcpy(env->isosurfaceField,sValue.c_str());

if(code==VV_SET_SLICEPLANE)
{
	if(sValue!="x" && sValue!="y" && sValue!="z")
	{
	  std::cerr<<"Invalid value VV_SET_SLICEPLANE "<<sValue<<std::endl;
	  env->setatt[code-VV_PARAM]=0;	  
	} else
	  
	strcpy(env->slicePlane,sValue.c_str());
}

if(code==VV_SET_SLICEPOS)
	strcpy(env->slicePosition,sValue.c_str());

if(code==VV_SET_SLICEFIELD)
	strcpy(env->slicefield,sValue.c_str());

if(code==VV_SET_VFIELD)
{
	std::stringstream sstmp(sValue);
	std::string stmp;

	sstmp >> stmp;
        strcpy(env->xVectorField, stmp.c_str());
	sstmp >> stmp;
        strcpy(env->yVectorField, stmp.c_str());
	sstmp >> stmp;
        strcpy(env->zVectorField, stmp.c_str());
}
if(code==VV_SET_OUT)
	strcpy(env->imageName,sValue.c_str());

if(code==VV_SET_CAMERA)
{
	std::stringstream sstmp(sValue);
	sstmp>>env->azimuth;
	sstmp>>env->elevation;
	sstmp>>env->zoom;
}
if(code==VV_SET_CAMPOS)
{
	std::stringstream sstmp(sValue);
	sstmp>>env->campos[0];
	sstmp>>env->campos[1];
	sstmp>>env->campos[2];
}
if(code==VV_SET_CAMFP)
{
	std::stringstream sstmp(sValue);
	sstmp>>env->camfp[0];
	sstmp>>env->camfp[1];
	sstmp>>env->camfp[2];
}
if(code==VV_SET_AZIMUTH)
{
	std::stringstream sstmp(sValue);
	sstmp>>env->azimuth;
}
if(code==VV_SET_ELEVATION)
{
	std::stringstream sstmp(sValue);
	sstmp>>env->elevation;
}
if(code==VV_SET_ZOOM)
{
	std::stringstream sstmp(sValue);
	sstmp>>env->zoom;
}
if(code==VV_SET_IMAGESIZE)
{
	if(sValue!="small" && sValue!="medium" && sValue!="large")
	{
	  std::cerr<<"Invalid value VV_SET_IMAGESIZE "<<sValue<<std::endl;
	  env->setatt[code-VV_PARAM]=0;
	} else
	  
	strcpy(env->imageSize,sValue.c_str());
}
if(code==VV_SET_BACKCOLOR)
{ 
	if(sValue!="yellow" && sValue!="red" && sValue!="green" 
	  && sValue!="blue" && sValue!="white"  && sValue!="black"
	  && sValue!="cyan" && sValue!="violet" )
	{
	  std::cerr<<"Invalid value VV_SET_BACKCOLOR "<<sValue<<std::endl;
	  env->setatt[code-VV_PARAM]=0;
	} else
	  
	strcpy(env->backColor,sValue.c_str());
}
if(code==VV_SET_ONECOLOR)
{ 
	if(sValue!="yellow" && sValue!="red" && sValue!="green" 
	  && sValue!="blue" && sValue!="white"  && sValue!="black"
	  && sValue!="cyan" && sValue!="violet" )
	{
	  std::cerr<<"Invalid value VV_SET_ONECOLOR "<<sValue<<std::endl;
	  env->setatt[code-VV_PARAM]=0;
	} else
	  
	strcpy(env->oneColor,sValue.c_str());
}
if(code==VV_SET_VO)
	strcpy(env->vo,sValue.c_str());
if(code==VV_SET_LFNOUT)
	strcpy(env->lfnout,sValue.c_str());
if(code==VV_SET_SE)
	strcpy(env->se,sValue.c_str());

if(code==VV_SET_COLORTABLE)
	strcpy(env->colorTable,sValue.c_str());

if(code==VV_SET_COLORRANGE)
{
	std::stringstream sstmp(sValue);
	sstmp>>env->colorRangeFrom;
	sstmp>>env->colorRangeTo;
}
if(code==VV_SET_COLORRANGEFROM)
{
	std::stringstream sstmp(sValue);
	sstmp>>env->colorRangeFrom;
}
if(code==VV_SET_COLORRANGETO)
{
	std::stringstream sstmp(sValue);
	sstmp>>env->colorRangeTo;
}
if(code==VV_SET_STEREO)
{
	if(sValue!="RedBlue" && sValue!="CrystalEyes" && sValue!="Anaglyph")
	{
	  std::cerr<<"Invalid value VV_SET_STEREO "<<sValue<<std::endl;
	  env->setatt[code-VV_PARAM]=0;
	} else
	  
	strcpy(env->stereoMode,sValue.c_str());
}
if(code==VV_SET_ANAGLYPHSAT)
{
	std::stringstream sstmp(sValue);
	sstmp>>env->anaglyphsat;
	if(env->anaglyphsat<0.0)env->anaglyphsat=0.0; 
	if(env->anaglyphsat>1.0)env->anaglyphsat=1.0; 
}
if(code==VV_SET_ANAGLYPHMASK)
	strcpy(env->anaglyphmask,sValue.c_str());

if(code==VV_SET_CYCLEOFFSET)
{
	std::stringstream sstmp(sValue);
	sstmp>>env->cycleOffset;
}
if(code==VV_SET_GLYPHS)
{   
	if(sValue!="pixel" && sValue!="sphere" && sValue!="cone" 
	  && sValue!="cylinder" && sValue!="cube")
	{
	  std::cerr<<"Invalid value VV_SET_GLYPHS "<<sValue<<std::endl;
	  env->setatt[code-VV_PARAM]=0;
	} else
	  
	strcpy(env->glyphs,sValue.c_str());
}
if(code==VV_SET_RADIUS)
{
	std::stringstream sstmp(sValue);
	sstmp>>env->radius;
}
if(code==VV_SET_HEIGHT)
{
	std::stringstream sstmp(sValue);
	sstmp>>env->height;
}
if(code==VV_SET_RADIUSSCALAR)
	strcpy(env->radiusscalar,sValue.c_str());
if(code==VV_SET_HEIGHTSCALAR)
	strcpy(env->heightscalar,sValue.c_str());

if(code==VV_SET_OPACITY)
{
	std::stringstream sstmp(sValue);
	sstmp>>env->opacity;
}
if(code==VV_SET_OPACITYTF)
{	
        std::stringstream sstmp(sValue);
	sstmp>>env->opacityTF[0];
	sstmp>>env->opacityTF[1];
	sstmp>>env->opacityTF[2];
}
if(code==VV_SET_ISOSURVALUE)
{
	std::stringstream sstmp(sValue);
	sstmp>>env->isosurfaceValue;
	if(env->isosurfaceValue>255) env->isosurfaceValue=255;
	if(env->isosurfaceValue<0) env->isosurfaceValue=0;
}
if(code==VV_SET_ISOSMOOTH)
	strcpy(env->isoSmooth,sValue.c_str());
if(code==VV_SET_ISOSURFACE)
{  
      env->setatt[VV_SET_VOLUMERENDERING-VV_PARAM]=0;
      env->setatt[VV_SET_SLICE-VV_PARAM]=0;

}
if(code==VV_SET_VOLUMERENDERING)
{  
      env->setatt[VV_SET_ISOSURFACE-VV_PARAM]=0;
      env->setatt[VV_SET_SLICE-VV_PARAM]=0;

}
if(code==VV_SET_SLICE)
{  
      env->setatt[VV_SET_ISOSURFACE-VV_PARAM]=0;
      env->setatt[VV_SET_VOLUMERENDERING-VV_PARAM]=0;

}
if(code==VV_SET_SLICEPLANEPOINT)
	strcpy(env->slicePlanePoint,sValue.c_str());

if(code==VV_SET_SLICEPLANENORMAL)
	strcpy(env->slicePlaneNormal,sValue.c_str());
if(code==VV_SET_VECTORSCALINGFACTOR)
{
	std::stringstream sstmp(sValue);
	sstmp>>env->vectorScalingFactor;
}
if(code==VV_SET_VECTORSCALE)
{
	std::stringstream sstmp(sValue);
	sstmp>>env->vectorScale;
	if(env->vectorScale!=0 ||env->vectorScale!=1 )
	{ 
	  std::cerr<<"Invalid VV_SET_VECTORSCALE value "<<sValue<<std::endl;
	  env->vectorScale=-1;
	  env->setatt[VV_SET_VECTORSCALE-VV_PARAM]=0;
	}
}
if(code==VV_SET_CLIPRANGE)
{
	std::stringstream sstmp(sValue);
	sstmp>>env->cliprange[0];
	sstmp>>env->cliprange[1];
}

//////////////////// SPLOTCH Section ///////////////



if(code==VV_SET_SPL_INTENSITYLOG0)
{
    std::ofstream splFile;
    splFile.open(env->splotchpar, std::ios::app);
    splFile<<"intensity_log0="<<sValue<<std::endl;
    splFile.close();
}
if(code==VV_SET_SPL_INTENSITYMIN0)
{
    std::ofstream splFile;
    splFile.open(env->splotchpar, std::ios::app);
    splFile<<"intensity_min0="<<sValue<<std::endl;
    splFile.close();
} 
if(code==VV_SET_SPL_INTENSITYMAX0)
{
    std::ofstream splFile;
    splFile.open(env->splotchpar, std::ios::app);
    splFile<<"intensity_max0="<<sValue<<std::endl;
    splFile.close();
} 
if(code==VV_SET_SPL_SIZEFIX0)
{
    std::ofstream splFile;
    splFile.open(env->splotchpar, std::ios::app);
    splFile<<"size_fix0="<<sValue<<std::endl;
    splFile.close();
} 
if(code==VV_SET_SPL_SIZEFAC0)
{
    std::ofstream splFile;
    splFile.open(env->splotchpar, std::ios::app);
    splFile<<"size_fac0="<<sValue<<std::endl;
    splFile.close();
} 
if(code==VV_SET_SPL_COLORISVECTOR0)
{
    std::ofstream splFile;
    splFile.open(env->splotchpar, std::ios::app);
    splFile<<"color_is_vector0="<<sValue<<std::endl;
    splFile.close();
} 
if(code==VV_SET_SPL_XRES)
{
    std::ofstream splFile;
    splFile.open(env->splotchpar, std::ios::app);
    splFile<<"xres="<<sValue<<std::endl;
    splFile.close();
} 
if(code==VV_SET_SPL_YRES)
{
    std::ofstream splFile;
    splFile.open(env->splotchpar, std::ios::app);
    splFile<<"yres="<<sValue<<std::endl;
    splFile.close();
} 
if(code==VV_SET_SPL_FOV)
{
    std::ofstream splFile;
    splFile.open(env->splotchpar, std::ios::app);
    splFile<<"fov="<<sValue<<std::endl;
    splFile.close();
} 
if(code==VV_SET_SPL_GRAYABSORBTION)
{
    std::ofstream splFile;
    splFile.open(env->splotchpar, std::ios::app);
    splFile<<"gray_absorption="<<sValue<<std::endl;
    splFile.close();
} 
if(code==VV_SET_SPL_BRIGHTNESS0)
{
    std::ofstream splFile;
    splFile.open(env->splotchpar, std::ios::app);
    splFile<<"brightness0="<<sValue<<std::endl;
    splFile.close();
} 
if(code==VV_SET_SPL_COLORLOG0)
{
    std::ofstream splFile;
    splFile.open(env->splotchpar, std::ios::app);
    splFile<<"color_log0="<<sValue<<std::endl;
    splFile.close();
} 
if(code==VV_SET_SPL_COLORASINH0)
{
    std::ofstream splFile;
    splFile.open(env->splotchpar, std::ios::app);
    splFile<<"color_asinh0="<<sValue<<std::endl;
    splFile.close();
} 

if(code==VV_SET_SPL_ROLSKY_X)
{
    std::ofstream splFile;
    splFile.open(env->splotchpar, std::ios::app);
    splFile<<"sky_x="<<sValue<<std::endl;
    splFile.close();
} 
if(code==VV_SET_SPL_ROLSKY_Y)
{
    std::ofstream splFile;
    splFile.open(env->splotchpar, std::ios::app);
    splFile<<"sky_y="<<sValue<<std::endl;
    splFile.close();
} 
if(code==VV_SET_SPL_ROLSKY_Z)
{
    std::ofstream splFile;
    splFile.open(env->splotchpar, std::ios::app);
    splFile<<"sky_z="<<sValue<<std::endl;
    splFile.close();
} 

////////////////////////////////////////////////////
if(code==VV_SET_CLEAN)  VV_Clean(env);

return noError;

}
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
	  && sValue !="ramses" && sValue !="changa" && sValue !="hecuba")
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
     errorCode=VV_SetAtt(VVenv,VV_SET_FILEVBT,tmp);
     errorCode=VV_View(VVenv);
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
} //extern "C"

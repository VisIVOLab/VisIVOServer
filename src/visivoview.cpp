#include "visivodef.h"
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <iostream>
#include "optionssetter.h"
#include "vtkGraphicsFactory.h"

extern "C"
{
int VV_View(VisIVOViewer *env);
void* VV_View_Thread(void *id);
void* VV_View_Process(void *id);
int createProcessView(VisIVOAsynchId *id);

#ifndef WIN32
	pthread_mutex_t mutex=PTHREAD_MUTEX_INITIALIZER;
#endif

//----------------------------
int VA_View(VisIVOViewer *e, VisIVOAsynchId *id)
//----------------------------
{
	int iret;
	#ifdef WIN32
 		iret=VV_View(e);
 		return iret;
	#endif
	#ifdef MAC
		id->withThread=1;
	#endif
	#ifndef WIN32
	id->env=e;
	id->state=runningThread;
	if(id->withThread>0)	iret=pthread_create(&(id->threadId), NULL, &VV_View_Thread, (void *)id);
	else	iret=createProcessView(id);//pthread_create(&(id->threadId), NULL, &VV_View_Process, (void *)id);
	if(iret!=0)
	{
		id->state=undefined;
		iret=errorThread;
	}
	return iret;
	#endif
}

//---------------------------
void *VV_View_Thread(void *id)
//---------------------------
{
	int iret;
	VisIVOAsynchId *idTh;
	idTh=(VisIVOAsynchId *)id;
	iret=VV_View((VisIVOViewer *)(idTh->env));
	idTh->errorCode=iret;
	idTh->state=successfulEndThread;
	return 0;
}
//---------------------------//---------------------------
int createProcessView(VisIVOAsynchId *id)
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
			iret=VV_View((VisIVOViewer *)(id->env));
			exit(iret);
		}
	}
	#endif
	return 0;	
}

//---------------------------
void *VV_View_Process(void *id)
//---------------------------
{
	int iret, i;
	VisIVOAsynchId *idTh;
	idTh=(VisIVOAsynchId *)id;
	#ifndef WIN32
	
	pid_t pid;
	// New process.
	pid=fork();
	if(pid==-1)
	{
		idTh->state=errorThread;
	}
	else
		if(pid!=0)
		{
			// Parent
			i=waitpid(pid, &iret, 0);
			idTh->state=successfulEndThread;
			if(WIFEXITED(iret))	idTh->errorCode=(int)WEXITSTATUS(iret);
			return 0;
		}
		else
		{
			// Child
			iret=VV_View((VisIVOViewer *)(idTh->env));
			exit(iret);
		}

	#endif
	return 0;	
}
	
//----------------------------
int VV_SetTable(VisIVOViewer *env, VSTable *table)
//----------------------------
{
    if (env == nullptr || table == nullptr)
        return invalidParCode;

    env->table = table;
    env->enableInMemory = true;

    return noError;
}

int VV_SetTableFromImporter(VisIVOViewer* viewer, VisIVOImporter* importer, size_t index) {
    if (!viewer || !importer) return -1;
    if (index >= importer->memTables->size() || index < 0) {
		std::cerr << "Desired table index not available" << std::endl;
		return -1;
	}
	
    VSTable* importerTable = importer->memTables->at(index);

    if (!importerTable) {
		std::cerr << "Table not defined" << std::endl;
		return -1;
	}

    viewer->table = importerTable;
    viewer->enableInMemory = true;
	return 0;
}

//----------------------------
int VV_View(VisIVOViewer *env)
//---------------------------
{
#ifndef WIN32
 pthread_mutex_lock(&mutex);
#endif
OptionsSetter *pOptSett= new OptionsSetter;
std::vector<std::string> args;
if(env->setatt[VV_SET_INTERNAL-VV_PARAM]==0 && env->setatt[VV_SET_FILEVBT-VV_PARAM]==0 && !env->enableInMemory)
{
  std::cerr<<"VV_View: Invalid VBT. Please set VI_SET_FILEVBT ";
  std::cerr<<"with VI_SettAtt function"<<std::endl;
  return invalidVBT;
}

if(env->setatt[VV_SET_INTERNAL-VV_PARAM]==1)
{
 // Inavlid external   
  env->setatt[VV_SET_FIELD-VV_PARAM]=0;
  env->setatt[VV_SET_COLORSCALAR-VV_PARAM]=0;
  env->setatt[VV_SET_VFIELD-VV_PARAM]=0;
  env->setatt[VV_SET_VRENDERFIELD-VV_PARAM]=0;
  env->setatt[VV_SET_ISOSURFIELD-VV_PARAM]=0;
  env->setatt[VV_SET_SLICEFIELD-VV_PARAM]=0;
  env->setatt[VV_SET_RADIUSSCALAR-VV_PARAM]=0;
  env->setatt[VV_SET_HEIGHTSCALAR-VV_PARAM]=0;
}

for(int idPar=0; idPar<NPAR; idPar++)
{
  if(idPar+VV_PARAM==VV_SET_DEFAULTIMAGES && env->setatt[idPar]==0)
    args.push_back("--nodefault");
 

  if(env->setatt[idPar]==1)
  {
//    std::clog<<" i="<<idPar<<" setatt "<<env->setatt[idPar]<<std::endl;
    switch(idPar+VV_PARAM)
    {
      case VV_SET_FILEVBT:
      {
	args.push_back(env->path);
	break;
      }
      
      case VV_SET_INTERNALFIELD:
      case VV_SET_FIELD:
      {
	args.push_back("--x");
	args.push_back(env->xField);
	args.push_back("--y");
	args.push_back(env->yField);
	args.push_back("--z");
	args.push_back(env->zField);
	break;
      }
      case VV_SET_COLOR:
      {
	args.push_back("--color");
	break;
      }
      case VV_SET_INTERNALCOLORSCALAR:
      case VV_SET_COLORSCALAR:
      {
	args.push_back("--colorscalar");
	args.push_back(env->colorScalar);
	break;
      }
      case VV_SET_BOX:
      {
	args.push_back("--showbox");
	break;
      }
      case VV_SET_AXES:
      {
	args.push_back("--showaxes");
	break;
      }
      case VV_SET_PALETTE:
      {
	args.push_back("--showlut");
	break;
      }
      case VV_SET_VOLUME:
      {
	args.push_back("--volume");
	break;
      }
      case VV_SET_INTERNALVRENDERFIELD:
      case VV_SET_VRENDERFIELD:
      {
	args.push_back("--vrenderingfield");
	args.push_back(env->vRenderingField);
	args.push_back("--isosurfacefield");
	args.push_back(env->isosurfaceField);
	args.push_back("--slicefield");
	args.push_back(env->slicefield);
	break;
      }
       case VV_SET_ISOSURFACE:
      {
	args.push_back("--isosurface");
	break;
      }
      case VV_SET_ISOSURFIELD:
      {
	args.push_back("--isosurfacefield");
	args.push_back(env->isosurfaceField);
	break;
      }
       case VV_SET_SLICE:
      {
	args.push_back("--slice");
	break;
      }
      case VV_SET_SLICEFIELD:
      {
	args.push_back("--slicefield");
	args.push_back(env->slicefield);
	break;
      }
      case VV_SET_SLICEPLANE:
      {
	args.push_back("--sliceplane");
	args.push_back(env->slicePlane);
	break;
      }
      case VV_SET_SLICEPOS:
      {
	args.push_back("--sliceposition");
	args.push_back(env->slicePosition);
	break;
      }
      case VV_SET_VECTOR:
      {
	args.push_back("--vector");
	break;
      }
      case VV_SET_INTERNALVFIELD:
      case VV_SET_VFIELD:
      {
	args.push_back("--vx");
	args.push_back(env->xVectorField);
	args.push_back("--vy");
	args.push_back(env->yVectorField);
	args.push_back("--vz");
	args.push_back(env->zVectorField);
	break;
      }
      case VV_SET_OUT:
      {
	args.push_back("--out");
	args.push_back(env->imageName);
	break;
      }
      case VV_SET_VO:
      {
	args.push_back("--VO");
	args.push_back(env->vo);
	break;
      }
      case VV_SET_LFNOUT:
      {
	args.push_back("--lfnout");
	args.push_back(env->lfnout);
	break;
      }
      case VV_SET_SE:
      {
	args.push_back("--se");
	args.push_back(env->se);
	break;
      }
      case VV_SET_CAMERA:
      {
	std::stringstream sstmp;
	sstmp << env->azimuth;
	std::string stmp;
	sstmp>>stmp;
	args.push_back("--camazim");
	args.push_back(stmp.c_str());

	sstmp.clear();
	stmp.clear();
	sstmp << env->elevation;
	sstmp>>stmp;
    
	args.push_back("--camelev");
	args.push_back(stmp.c_str());

	sstmp.clear();
	stmp.clear();
	sstmp << env->zoom;
	sstmp>>stmp;
	
	args.push_back("--zoom");
	args.push_back(stmp.c_str());
	break;
      }
      case VV_SET_CAMPOS:
      {
	std::stringstream sstmp;
	sstmp << env->campos[0]<<" ";
	sstmp << env->campos[1]<<" ";
	sstmp << env->campos[2];
	std::string stmp;
	sstmp>>stmp;
	args.push_back("--campos");
	args.push_back(stmp.c_str());
	break;
      }
      case VV_SET_CAMFP:
      {
	std::stringstream sstmp;
	sstmp << env->camfp[0]<<" ";
	sstmp << env->camfp[1]<<" ";
	sstmp << env->camfp[2];
	std::string stmp;
	sstmp>>stmp;
	args.push_back("--camfp");
	args.push_back(stmp.c_str());
	break;
      }
      case VV_SET_AZIMUTH:
      {
	std::stringstream sstmp;
	sstmp << env->azimuth;
	std::string stmp;
	sstmp>>stmp;
	args.push_back("--camazim");
	args.push_back(stmp.c_str());
	break;
      }
      case VV_SET_ELEVATION:
      {
	std::stringstream sstmp;
	sstmp << env->elevation;
	std::string stmp;
	sstmp>>stmp;
	args.push_back("--camelev");
	args.push_back(stmp.c_str());
	break;
      }
      case VV_SET_ZOOM:
      {
	std::stringstream sstmp;
	sstmp << env->zoom;
	std::string stmp;
	sstmp>>stmp;
	args.push_back("--zoom");
	args.push_back(stmp.c_str());
	break;
      }
      case VV_SET_IMAGESIZE:
      {
	args.push_back("--imagesize");
	args.push_back(env->imageSize);
	break;
      }
      case VV_SET_BACKCOLOR:
      {
	args.push_back("--backcolor");
	args.push_back(env->backColor);
	break;
      }
      case VV_SET_ONECOLOR:
      {
	args.push_back("--onecolor");
	args.push_back(env->oneColor);
	break;
      }
      case VV_SET_COLORTABLE:
      {
	args.push_back("--colortable");
	args.push_back(env->colorTable);
	break;
      }
      case VV_SET_COLORRANGE:
      {
	std::stringstream sstmp;
	sstmp << env->colorRangeFrom;
	std::string stmp;
	sstmp>>stmp;
	args.push_back("--colorrangefrom");
	args.push_back(stmp.c_str());

	sstmp.clear();
	stmp.clear();
	sstmp << env->colorRangeTo;
	sstmp>>stmp;
    
	args.push_back("--colorrangeto");
	args.push_back(stmp.c_str());
	break;
      }
      case VV_SET_COLORRANGEFROM:
      {
	std::stringstream sstmp;
	sstmp << env->colorRangeFrom;
	std::string stmp;
	sstmp>>stmp;
	args.push_back("--colorrangefrom");
	args.push_back(stmp.c_str());
	break;
      }
      case VV_SET_COLORRANGETO:
      {
	std::stringstream sstmp;
	sstmp << env->colorRangeTo;
	std::string stmp;
	sstmp>>stmp;
	args.push_back("--colorrangeto");
	args.push_back(stmp.c_str());
	break;
      }
      case VV_SET_STEREO:
      {
	args.push_back("--stereo");
	args.push_back(env->stereoMode);
	break;
      }
      case VV_SET_ANAGLYPHSAT:
      {
	std::stringstream sstmp;
	sstmp << env->anaglyphsat;
	std::string stmp;
	sstmp>>stmp;
	args.push_back("--anaglyphsat");
	args.push_back(stmp.c_str());
	break;
      }
      case VV_SET_ANAGLYPHMASK:
      {
	args.push_back("--anaglyphmask");
	args.push_back(env->anaglyphmask);
	break;
      }
      case VV_SET_CYCLE:
      {
	args.push_back("--cycle");
	args.push_back(env->cycleFile);
	if(env->setatt[VV_SET_CYCLEOFFSET-VV_PARAM]==0)
	{ 
	  std::stringstream sstmp;
	  sstmp << env->cycleOffset;
	  std::string stmp;
	  sstmp>>stmp;
	  args.push_back("--cycleoffset");
	  args.push_back(stmp.c_str());
	}
	break;
      }
      case VV_SET_CYCLEOFFSET:
      {
	std::stringstream sstmp;
	sstmp << env->cycleOffset;
	std::string stmp;
	sstmp>>stmp;
	args.push_back("--cycleoffset");
	args.push_back(stmp.c_str());
	break;
      }
      case VV_SET_SCALE:
      {
	args.push_back("--scale");
	break;
      }
      case VV_SET_LOGSCALE:
      {
	args.push_back("--logscale");
	break;
      }
	  case VV_SET_USEMEMORY:
      {
	args.push_back("--usememory");
	break;
      }
      case VV_SET_GLYPHS:
      {
	args.push_back("--glyphs");
	args.push_back(env->glyphs);
	break;
      }
      case VV_SET_RADIUS:
      {
	std::stringstream sstmp;
	sstmp << env->radius;
	std::string stmp;
	sstmp>>stmp;
	args.push_back("--radius");
	args.push_back(stmp.c_str());
	break;
      }
      case VV_SET_HEIGHT:
      {
	std::stringstream sstmp;
	sstmp << env->height;
	std::string stmp;
	sstmp>>stmp;
	args.push_back("--height");
	args.push_back(stmp.c_str());
	break;
      }
      case VV_SET_SCALEGLYPHS:
      {
	args.push_back("--scaleglyphs");
	break;
      }
      case VV_SET_INTERNALRADIUSSCALAR:
      case VV_SET_RADIUSSCALAR:
      {
	args.push_back("--radiusscalar");
	args.push_back(env->radiusscalar);
	break;
      }
      case VV_SET_INTERNALHEIGHTSCALAR:
      case VV_SET_HEIGHTSCALAR:
      {
	args.push_back("--heightscalar");
	args.push_back(env->heightscalar);
	break;
      }
      case VV_SET_OPACITY:
      {
	std::stringstream sstmp;
	sstmp << env->opacity;
	std::string stmp;
	sstmp>>stmp;
	args.push_back("--opacity");
	args.push_back(stmp.c_str());
	break;
      }
      case VV_SET_OPACITYTF:
      {
	args.push_back("--opacityTF");
	std::stringstream sstmp;
	sstmp << env->opacityTF[0];
	std::string stmp;
	sstmp>>stmp;
	args.push_back(stmp.c_str());
	sstmp.clear();
	sstmp << env->opacityTF[1];
	sstmp>>stmp;
	args.push_back(stmp.c_str());
	sstmp.clear();
	sstmp << env->opacityTF[2];
	sstmp>>stmp;
	args.push_back(stmp.c_str());
	break;
      }
      case VV_SET_SHADOW:
      {
	args.push_back("--shadow");
	break;
      }
      case VV_SET_ISOSURVALUE:
      {
	std::stringstream sstmp;
	sstmp << env->isosurfaceValue;
	std::string stmp;
	sstmp>>stmp;
	args.push_back("--isosurfacevalue");
	args.push_back(stmp.c_str());
	break;
      }
       case VV_SET_WIREFRAME:
      {
	args.push_back("--wireframe");
	break;
      }
      case VV_SET_ISOSMOOTH:
      {
	args.push_back("--isosmooth");
	args.push_back(env->isoSmooth);
	break;
      }
      case VV_SET_SLICEPLANEPOINT:
      {
	args.push_back("--sliceplanepoint");
	args.push_back(env->slicePlanePoint);
	break;
      }
      case VV_SET_SLICEPLANENORMAL:
      {
	args.push_back("--sliceplanenormal");
	args.push_back(env->slicePlaneNormal);
	break;
      }
      case VV_SET_VECTORLINE:
      {
	args.push_back("--vectorline");
	break;
      }
      case VV_SET_VECTORSCALINGFACTOR:
      {
	std::stringstream sstmp;
	sstmp << env->vectorScalingFactor;
	std::string stmp;
	sstmp>>stmp;
	args.push_back("--vectorscalefactor");
	args.push_back(stmp.c_str());
	break;
      }
      case VV_SET_VECTORSCALE:
      {
	std::stringstream sstmp;
	sstmp << env->vectorScale;
	std::string stmp;
	sstmp>>stmp;
	args.push_back("--vectorscale");
	args.push_back(stmp.c_str());
	break;
      }
      case VV_SET_SPLOTCH:
      {
	args.push_back("--splotch");
	args.push_back(env->splotchpar);
	break;	
      }
    case VV_SET_CLIPRANGE:
      {
	std::stringstream sstmp;
	sstmp << env->cliprange[0]<<" ";
	sstmp << env->cliprange[1]<<" ";
	std::string stmp;
	sstmp>>stmp;
	args.push_back("--cliprange");
	args.push_back(stmp.c_str());
	break;
      }
    case VV_SET_CLIPLARGE:
      {
	args.push_back("--cliplarge");
	break;
      }

       default:
     {
	break;
      }
    } 
    
  }//end if env->
}  //end for
// for(int i=0;i<args.size();i++)std::clog<<args[i]<<" "; 
//  std::clog<<std::endl<<"test VisIVO View"<<std::endl;
 
  if(env->setatt[VV_SET_INTERNAL-VV_PARAM]==1)
    if(!pOptSett->setInternalData(env))
    {
    	#ifndef WIN32
 			pthread_mutex_unlock(&mutex);
	    #endif
		return invalidInternalData;
    }  
  
  pOptSett->parseOption(args);
  VisIVOServerOptions opt=pOptSett->returnOptions();
  
  if(env->enableInMemory) pOptSett->setTable(env->table);
  //if(!pOptSett->internalData() && !pOptSett->getUseMemory()) pOptSett->readData(); 
  if(pOptSett->images()<0)
  { 
	std::clog << "Error" << std::endl;
      delete pOptSett;
      
      #ifndef WIN32
 		pthread_mutex_unlock(&mutex);
	  #endif
      
      return invalidImage;
  }

 if(env->setatt[VV_SET_CYCLE-VV_PARAM]==1)
  {
    VisIVOServerOptions VVoptions=pOptSett->returnOptions();
//    std::clog<<"Check ==> "<<VVoptions.numImageToLoad<<std::endl;
    env->cycleOffset+=VVoptions.numImageToLoad;
    
  }

  if ( pOptSett!=0)
    delete pOptSett;
#ifndef WIN32
 pthread_mutex_unlock(&mutex);
#endif
  return noError;
}
/////////////////////////////
} //extern "C"

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


#include "optionssetter.h"

#include "visivoutils.h"
#include "pipe.h"
#include "pointspipe.h"

#ifdef SPLVISIVO

#include "splotchpipecamera.h"
#include "splotchpipe.h"

#endif

#include "volumepipe.h"
#include "slicerpipe.h"
#include "isosurfacepipe.h"
#include "vectorpipe.h"
#include "vtkimagepipe.h"
#include "pointssmoothpipe.h"
#include "glitegw.h"

#include <cstdlib>
#include <cstring>

#include <sstream>
#include <iostream>
#include <fstream>
#include <map>
#include "parfile.h"

#include "historyXmlWriter.h"

#include <math.h>


const unsigned int OptionsSetter::MAX_ELEMENT_TO_LOAD = 400000000;
const float OptionsSetter::INVALID_READ = -100000000.0;
const double OptionsSetter::INVALID_DREAD = -100000000.0;
const double OptionsSetter::INVALID_CAM = -123456789.31;


//---------------------------------------------------------------------
OptionsSetter::OptionsSetter ( )
//---------------------------------------------------------------------
{
    
    m_vServer.path="none";
    
    m_vServer.xField=m_vServer.yField=m_vServer.zField="none";
    m_vServer.xVectorField=m_vServer.yVectorField=m_vServer.zVectorField="none";
    m_vServer.hsmlField=m_vServer.spColorField=m_vServer.spIntensityField="none";
    m_vServer.sprField=m_vServer.spIField=m_vServer.spC1Field=m_vServer.spC2Field=m_vServer.spC3Field="none";
    m_vServer.colorScalar="none";
    
    m_vServer.imageName="VisIVOServerImage";
    
    m_vServer.scale="none";
    
    m_vServer.noDefault="none";
    m_vServer.uselogscale="none";
    
    m_vServer.scaleGlyphs="none";
    m_vServer.radiusscalar="none";
    m_vServer.heightscalar="none";
    
    m_vServer.colorTable="none";
    m_vServer.glyphs="pixel";
    
    m_vServer.savepar="none";
    m_vServer.loadpar="none";
    
    m_vServer.opacity=0.666;
    
    m_vServer.elevation=0;
    m_vServer.azimuth =0;
    m_vServer.zoom=1;
    m_vServer.fov=45;
    m_vServer.fovIsGiven=false;
    m_vServer.cameraPos[0]=0;
    m_vServer.cameraPos[1]=0;
    m_vServer.cameraPos[2]=0;
    m_vServer.cameraFocalPoint[0]=0;
    m_vServer.cameraFocalPoint[1]=0;
    m_vServer.cameraFocalPoint[2]=0;
    m_vServer.cameraRoll=0;
    m_vServer.cameraPosPrev=new double[3];
    m_vServer.cameraPosPrev[0]=INVALID_CAM;
    m_vServer.cameraPosPrev[1]=INVALID_CAM;
    m_vServer.cameraPosPrev[2]=INVALID_CAM;
    m_vServer.cameraFocalPointPrev=new double[3];
    m_vServer.cameraFocalPointPrev[0]=INVALID_CAM;
    m_vServer.cameraFocalPointPrev[1]=INVALID_CAM;
    m_vServer.cameraFocalPointPrev[2]=INVALID_CAM;
    m_vServer.cameraRollPrev=new double[1];
    m_vServer.cameraRollPrev[0]=INVALID_CAM;
    m_vServer.setCameraPos=false;
    m_vServer.setCameraFocalPoint=false;
    m_vServer.setCameraRoll=false;
    m_vServer.cycle=false;
    m_vServer.cycleFile="none";
    m_vServer.cycleOffset=0;
    m_vServer.cycleSkipFrom=-100000010;
    m_vServer.cycleSkipTo=-100000010;
    
    m_vServer.radius=1;
    m_vServer.height=1;
    
    
    m_vServer.x=0;
    m_vServer.y=1;
    m_vServer.z=2;
    
    m_vServer.vx=0;
    m_vServer.vy=0;
    m_vServer.vz=0;
    
    m_vServer.hsml=-1;
    m_vServer.spColor=-1;
    m_vServer.spIntensity=-1;
    m_vServer.spr=-1;
    m_vServer.spI=-1;
    m_vServer.spC1=-1;
    m_vServer.spC2=-1;
    m_vServer.spC3=-1;
    
    
    m_vServer.nRows=0;
    m_vServer.nCols=0;
    m_vServer.nColorScalar=0;
    m_vServer.nColorTable=0;
    
    
    m_vServer.endian="little";
    m_vServer.dataType="float";
    
    m_vServer.nHeight=0;
    m_vServer.nRadius=0;
    
    m_vServer.color="none";
    m_vServer.vector="none";
    
    m_vServer.volume="none";
    m_vServer.vRendering="none";
    m_vServer.vShadow=false;
    m_vServer.vRenderingField="none";
    m_vServer.nVRenderingField=0;
    m_vServer.slice="none";
    m_vServer.sliceOrthoNormal=false;
    m_vServer.sliceField="none";
    m_vServer.nSliceField=0;
    m_vServer.slicePlane="none";
    m_vServer.nSlicePlane=-1;
    m_vServer.slicePosition=0;
    m_vServer.slicePlaneNormal="none";
    m_vServer.slicePlanePoint="none";
    m_vServer.slicePlaneNormalValue[0]=1;
    m_vServer.slicePlaneNormalValue[1]=1;
    m_vServer.slicePlaneNormalValue[2]=1;
    m_vServer.slicePlanePointValue[0]=0;
    m_vServer.slicePlanePointValue[1]=0;
    m_vServer.slicePlanePointValue[2]=0;
    
    m_vServer.isosurface="none";
    m_vServer.isosurfaceField="none";
    m_vServer.nIsosurfaceField=0;
    m_vServer.isosurfaceValue= 0.11;
    
    m_vServer.needSwap=false;
    m_vServer.splotch=false;
    m_vServer.systemEndianism="little";
    m_vServer.dataRead=true;
    bool goodAllocation=false;
    m_vServer.tableData=NULL;
    m_vServer.colorRangeTo=-1;
    m_vServer.colorRangeFrom=-1;
    m_vServer.isColorRangeTo=false;
    m_vServer.isColorRangeFrom=false;
    m_vServer.oneColor="white";
    m_vServer.showBox=false;
    m_vServer.showLut=false;
    m_vServer.showAxes=false;
    m_vServer.wireframe=false;
    m_vServer.vectorLine=false;
    m_vServer.backColor="black";
    m_vServer.imageSize="medium";
    m_vServer.isoSmooth="none";
    m_vServer.vectorScalingFactor=INVALID_DREAD;
    m_vServer.vectorScale=-1;
    
    m_vServer.vtkImage="none";
    m_vServer.vtkImSize="none";
    m_vServer.vtkImSizeValue[0]=128;
    m_vServer.vtkImSizeValue[1]=128;
    m_vServer.vtkImSizeValue[2]=128;
    m_vServer.vtkGausExp=-1;
    m_vServer.vtkSpacingFact=1;
    m_vServer.vtkScale=1;
    m_vServer.vtkEcc=1;
    
    m_vServer.mode="none";
    m_vServer.scenario="etna";
    m_vServer.opacityTF[0]=5;
    m_vServer.opacityTF[1]=3;
    m_vServer.opacityTF[2]=2.5;
    
    m_vServer.stereo=false;
    m_vServer.stereoMode="none";
    m_vServer.anaglyphsat=0.65;
    m_vServer.anaglyphmask="4 3";
    m_vServer.internalData=false;
    m_vServer.colorIsFile=false;
    
    m_vServer.VO="none";
    m_vServer.lfnout="none";
    m_vServer.se="";
    m_inputLfnGiven=false;
    
    m_vServer.clipset=false;
    m_vServer.cliprange[0]=0.0;
    m_vServer.cliprange[1]=1.0e+3;
    
    m_vServer.autocalcCam=false;
    
    fileIsAtTheEnd=true;
    
}
//---------------------------------------------------------------------
OptionsSetter::~OptionsSetter ( )
//---------------------------------------------------------------------
{
    delete []  m_vServer.cameraPosPrev;
    delete []  m_vServer.cameraFocalPointPrev;
    delete []  m_vServer.cameraRollPrev;
}
//---------------------------------------------------------------------
bool OptionsSetter::readSplocthColumn()
//---------------------------------------------------------------------
{
    std::ifstream splFile(m_vServer.splotchpar.c_str());
    if(!splFile)
    {
        std::cerr<<"Invalid Splotch parameter file"<<std::endl;
        return -1;
    }
    Parfile params (m_vServer.splotchpar);
    if(!params.valid())
    {
        std::cerr<<"Invalid Splotch parameter filename"<<std::endl;
        return false;
    }
    m_vServer.spIntensityField = params.find<std::string>("label_intensity","none");
    m_vServer.spColorField = params.find<std::string>("label_color","none");
    m_vServer.hsmlField = params.find<std::string>("label_hsml","none");
    m_vServer.sprField= params.find<std::string>("r_col","none");
    m_vServer.spIField= params.find<std::string>("I_col","none");
    m_vServer.spC1Field= params.find<std::string>("C1_col","none");
    m_vServer.spC2Field= params.find<std::string>("C2_col","none");
    m_vServer.spC3Field= params.find<std::string>("C3_col","none");
    
    

    std::istringstream is(params.find<std::string>("autocalc_cam","false"));
    bool b;
    is >> std::boolalpha >> b;
    
    m_vServer.autocalcCam= b;


    if(m_vServer.colorScalar=="none")
        m_vServer.colorScalar=m_vServer.spC1Field;
    
    splFile.close();
    return true;
    
}
//---------------------------------------------------------------------
int OptionsSetter::readData ( )
//---------------------------------------------------------------------
{
	if(!m_vServer.dataRead || m_vServer.useMemory)
		return -1;
  	std::ifstream inFile;
  	inFile.open(m_vServer.path.c_str(), ios::binary);
  	if(!inFile)
	{
		std::cerr<<"Error: binary data does not exist"<<std::endl;
        return -1;
	}
    
	std::map<std::string, int>::iterator p;
	std::map<std::string, int> colNames;
    
	if(m_vServer.xField!="none")
		colNames.insert(make_pair(m_vServer.xField,0));
	if(m_vServer.yField!="none")
		colNames.insert(make_pair(m_vServer.yField,0));
	if(m_vServer.zField!="none")
		colNames.insert(make_pair(m_vServer.zField,0));
	if(m_vServer.xVectorField!="none")
		colNames.insert(make_pair(m_vServer.xVectorField,0));
	if(m_vServer.yVectorField!="none")
		colNames.insert(make_pair(m_vServer.yVectorField,0));
	if(m_vServer.zVectorField!="none")
		colNames.insert(make_pair(m_vServer.zVectorField,0));
	if(m_vServer.colorScalar!="none")
		colNames.insert(make_pair(m_vServer.colorScalar,0));
	if(m_vServer.radiusscalar!="none")
		colNames.insert(make_pair(m_vServer.radiusscalar,0));
	if(m_vServer.heightscalar!="none")
		colNames.insert(make_pair(m_vServer.heightscalar,0));
	if(m_vServer.vRenderingField!="none")
		colNames.insert(make_pair(m_vServer.vRenderingField,0));
	if(m_vServer.sliceField!="none")
		colNames.insert(make_pair(m_vServer.sliceField,0));
	if(m_vServer.isosurfaceField!="none")
		colNames.insert(make_pair(m_vServer.isosurfaceField,0));
	if(m_vServer.hsmlField!="none")
		colNames.insert(make_pair(m_vServer.hsmlField,0));
	if(m_vServer.spColorField!="none")
		colNames.insert(make_pair(m_vServer.spColorField,0));
	if(m_vServer.spIntensityField!="none")
		colNames.insert(make_pair(m_vServer.spIntensityField,0));
	if(m_vServer.sprField!="none")
		colNames.insert(make_pair(m_vServer.sprField,0));
	if(m_vServer.spIField!="none")
		colNames.insert(make_pair(m_vServer.spIField,0));
	if(m_vServer.spC1Field!="none")
		colNames.insert(make_pair(m_vServer.spC1Field,0));
	if(m_vServer.spC2Field!="none")
		colNames.insert(make_pair(m_vServer.spC2Field,0));
	if(m_vServer.spC3Field!="none")
		colNames.insert(make_pair(m_vServer.spC3Field,0));
    
	if(colNames.size()*m_vServer.nRows>MAX_ELEMENT_TO_LOAD)
		return -1;
    
    m_vServer.tableData= new float*[colNames.size()];
	m_vServer.goodAllocation=true;
	int i=-1;
	for(p=colNames.begin();p!=colNames.end();p++)
	{
		i++;
        try
        {
            m_vServer.tableData[i]=new  float[m_vServer.nRows];
        }
        catch(std::bad_alloc e)
        {
            m_vServer.tableData[i]=NULL;
        }
		if(m_vServer.tableData[i]==NULL)
        {
            for(unsigned int j=0;j<i;j++)
                delete [] m_vServer.tableData[j];
            delete [] m_vServer.tableData;
            m_vServer.goodAllocation=false;
            return -1;
        }
		int numOfCol=-1;
		for(int k=0;k<m_vServer.fieldNames.size();k++)
		{
			if(m_vServer.fieldNames[k]==p->first)
				numOfCol=k;
		}
        if(numOfCol>=0)
		{
            
		} else
		{
			std::cerr<<"Severe Warning: wrong read of data, "<< p->first<<" is not in the in table. No data will be pre read. Error can occur."<<std::endl;
			for(unsigned int j=0;j<i;j++)
				delete [] m_vServer.tableData[j];
			delete [] m_vServer.tableData;
			m_vServer.goodAllocation=false;
			return -1;
		}
		inFile.seekg(numOfCol * m_vServer.nRows* sizeof(float),std::ios::beg);
		inFile.read((char *)(m_vServer.tableData[i]),  sizeof(float)* m_vServer.nRows);
		m_vServer.columns.insert(make_pair(p->first,i));
	}
    // 	for(p=m_vServer.columns.begin();p!=m_vServer.columns.end();p++)
    // 		std::clog<<"mappa "<<p->first<<" " <<p->second<<std::endl;
    // swap float!
	if(m_vServer.needSwap)
	{
		for(int i=0;i<colNames.size();i++)
			for(int j=0;j<m_vServer.nRows;j++)
                m_vServer.tableData[i][j]=floatSwap((char *)(&m_vServer.tableData[i][j]));
	}
    
    return 0;
}
//---------------------------------------------------------------------
int OptionsSetter::parseOption (const std::vector<std::string>  arguments )
//---------------------------------------------------------------------
{
    std::stringstream ss;
    int i;
    std::string paramFilename="none";
    for (i=0; i<arguments.size(); i++)
    {
        if(arguments[i]=="--help")
        {
            showHelp();
            return -1;
        }
    }
    
    if(arguments.size()==1)
    {
        paramFilename=arguments[0];
        std::string tmp="none";
        Parfile params(paramFilename);
        if(!params.valid())
        {
            std::cerr<<"Invalid  parameter filename"<<std::endl;
            return -1;
        }
        
        m_vServer.volume = params.find<std::string>("volume","none");
        m_vServer.vector = params.find<std::string>("vector","none");
        tmp="none";
        tmp = params.find<std::string>("splotch","none");
        if(tmp=="yes")
        {
            m_vServer.splotch=true;
            m_vServer.splotchpar=paramFilename;
            m_vServer.vRendering="none";
            m_vServer.volume="none";
            
            viewParameter["splotch"]=m_vServer.splotchpar;
            
        }
        
        if(!params.param_present("input"))
        {
            std::cerr<<"No input filename is given."<<std::endl;
            return -1;
        }
        m_vServer.path=params.find<std::string>("input","none");
        //     std::clog<<m_vServer.path<<std::endl;
        m_vServer.VO=params.find<std::string>("VO","none");
        if(params.param_present("VO"))
            viewParameter["VO"]=m_vServer.VO;
        
        
        m_vServer.imageName=params.find<std::string>("out","VisIVOServerImage");
        
        m_vServer.lfnout=params.find<std::string>("lfnout","none");
        if(params.param_present("lfnout"))
            viewParameter["lfnout"]=m_vServer.lfnout;
        
        m_vServer.se=params.find<std::string>("se","");
        if(params.param_present("se"))
            viewParameter["se"]=m_vServer.se;
        
        m_vServer.cycleFile = params.find<std::string>("cycle","none");
        if(m_vServer.cycleFile!="none")
        {
            viewParameter["cycle"]=m_vServer.cycleFile;
            m_vServer.cycle=true;
            
        }
        m_vServer.cycleOffset=params.find<int>("cycleoffset",0);
        if(params.param_present("cycleoffset"))
        {
            std::stringstream ss;//create a stringstream
            ss << m_vServer.cycleOffset;//add number to the stream

            viewParameter["cycleoffset"]=ss.str();
        }
        //m_vServer.cycleOffset=params.find<int>("cycleoffset",0);
        m_vServer.cycleSkipFrom=params.find<int>("cycle_skip_from",-100000010);
        if(params.param_present("cycle_skip_from"))
        {
            
            std::stringstream ss;//create a stringstream
            ss << m_vServer.cycleSkipFrom;//add number to the stream

            viewParameter["cycle_skip_from"]=ss.str();
        }
        
        m_vServer.cycleSkipTo=params.find<int>("cycle_skip_to",-100000010);
        if(params.param_present("cycle_skip_to"))
        {
            std::stringstream ss;//create a stringstream
            ss << m_vServer.cycleSkipTo;//add number to the stream
            
            viewParameter["cycle_skip_to"]=ss.str();
        }
        
        if(params.param_present("colorrangeto"))
        {
            m_vServer.isColorRangeTo=true;
            m_vServer.colorRangeTo=params.find<float>("colorrangeto",-1);
            
            std::stringstream ss;//create a stringstream
            ss << m_vServer.colorRangeTo;//add number to the stream

            
            viewParameter["colorrangeto"]=ss.str();
        }
        if(params.param_present("colorrangefrom"))
        {
            m_vServer.isColorRangeFrom=true;
            m_vServer.colorRangeFrom=params.find<float>("colorrangefrom",-1);
            
            
            std::stringstream ss;//create a stringstream
            ss << m_vServer.colorRangeFrom;//add number to the stream

            viewParameter["colorrangefrom"]=ss.str();
        }
        
        if(params.param_present("stereo"))
        {
            m_vServer.stereoMode=params.find<std::string>("stereo","CrystalEyes");
            if(m_vServer.stereoMode == "CrystalEyes" ||  m_vServer.stereoMode =="RedBlue"||  m_vServer.stereoMode =="Anaglyph"){
                m_vServer.stereo=true;
                viewParameter["stereo"]=m_vServer.stereoMode;
            }
            else
                std::cerr<<"Invalid stereo option "<<m_vServer.stereoMode<<std::endl;
            
        }
        
        m_vServer.anaglyphsat=params.find<float>("anaglyphsat",0.65);
        if(params.param_present("anaglyphsat"))
        {
            std::stringstream ss;//create a stringstream
            ss << m_vServer.anaglyphsat;//add number to the stream

            viewParameter["anaglyphsat"]=ss.str();
        }
        
        m_vServer.anaglyphmask=params.find<std::string>("anaglyphmask","4 3");
        if(m_vServer.anaglyphmask=="")m_vServer.anaglyphmask="4 3";
        if(params.param_present("anaglyphmask"))
            viewParameter["anaglyphmask"]=m_vServer.anaglyphmask;
        
        
        
        if(params.param_present("x"))
        {
            m_vServer.xField=params.find<std::string>("x","X");
            viewParameter["x"]=m_vServer.xField;
        }
        if(params.param_present("x_coordinate"))
        {
            m_vServer.xField=params.find<std::string>("x_coordinate","X");
            viewParameter["x"]=m_vServer.xField;
            
        }
        if(params.param_present("y"))
        {
            m_vServer.yField=params.find<std::string>("y","Y");
            viewParameter["y"]=m_vServer.yField;
            
        }
        if(params.param_present("y_coordinate"))
        {
            m_vServer.yField=params.find<std::string>("y_coordinate","Y");
            viewParameter["y_coordinate"]=m_vServer.yField;
            
        }
        if(params.param_present("z"))
        {
            m_vServer.zField=params.find<std::string>("z","Z");
            viewParameter["z"]=m_vServer.zField;
            
        }
        
        if(params.param_present("z_coordinate"))
        {
            m_vServer.zField=params.find<std::string>("z_coordinate","Z");
            viewParameter["z_coordinate"]=m_vServer.zField;
            
        }
        
        if(params.param_present("onecolor"))
        {
            m_vServer.oneColor=params.find<std::string>("onecolor","white");
            viewParameter["onecolor"]=m_vServer.oneColor;
            
        }
        if(params.param_present("showbox"))
        {
            tmp =params.find<std::string>("showbox","no");
            if(tmp=="yes"){
                m_vServer.showBox=true;
                viewParameter["showbox"]="yes";
            }
        }
        if(params.param_present("showlut"))
        {
            tmp =params.find<std::string>("showlut","no");
            if(tmp=="yes")
            {
                m_vServer.showLut=true;
                viewParameter["showlut"]="yes";
            }
        }
        if(params.param_present("backcolor"))
        {
            m_vServer.backColor=params.find<std::string>("backcolor","black");
            viewParameter["backcolor"]=m_vServer.backColor;
            
        }
        
        if(params.param_present("imagesize"))
        {
            m_vServer.imageSize =params.find<std::string>("imagesize","medium");
            viewParameter["imagesize"]=m_vServer.backColor;
            
        }
        
        
        if(params.param_present("showaxes"))
        {
            tmp =params.find<std::string>("showaxes","no");
            if(tmp=="yes")
            {
                m_vServer.showAxes=true;
                viewParameter["showaxes"]="yes";
                
            }
        }
        if(params.param_present("wireframe"))
        {
            tmp =params.find<std::string>("wireframe","no");
            if(tmp=="yes")
            {
                m_vServer.wireframe=true;
                viewParameter["wireframe"]="yes";
                
            }
        }
        if(params.param_present("vectorline"))
        {
            tmp =params.find<std::string>("vectorline","no");
            if(tmp=="yes")
            {
                m_vServer.vectorLine=true;
                viewParameter["vectorline"]="yes";
                
            }
        }
        if(params.param_present("isosmooth"))
        {
            m_vServer.isoSmooth =params.find<std::string>("isosmooth","none");
            if ( m_vServer.isoSmooth != "none") {
                viewParameter["isosmooth"]=m_vServer.isoSmooth;
            }
            
        }
        if(params.param_present("vx"))
        {
            m_vServer.xVectorField=params.find<std::string>("vx","VX");
            if ( m_vServer.xVectorField != "VX") {
                viewParameter["vx"]=m_vServer.xVectorField;
            }
            
        }
        
        if(params.param_present("vx_coordinate"))
        {
            m_vServer.xVectorField=params.find<std::string>("vx_coordinate","VX");
            if ( m_vServer.xVectorField != "VX") {
                viewParameter["vx"]=m_vServer.xVectorField;
            }
            
        }
        
        if(params.param_present("vy"))
        {
            m_vServer.yVectorField=params.find<std::string>("vy","VY");
            if ( m_vServer.yVectorField != "VY") {
                viewParameter["vy"]=m_vServer.yVectorField;
            }
            
        }
        
        if(params.param_present("vy_coordinate"))
        {
            m_vServer.yVectorField=params.find<std::string>("vy_coordinate","VY");
            if ( m_vServer.yVectorField != "VY") {
                viewParameter["vy"]=m_vServer.yVectorField;
            }
            
        }
        
        
        if(params.param_present("vz"))
        {
            m_vServer.zVectorField=params.find<std::string>("vz","VZ");
            if ( m_vServer.zVectorField != "VZ") {
                viewParameter["vz"]=m_vServer.zVectorField;
            }
            
            
        }
        if(params.param_present("vz_coordinate"))
        {
            m_vServer.zVectorField=params.find<std::string>("vz_coordinate","VZ");
            if ( m_vServer.zVectorField != "VZ") {
                viewParameter["vz"]=m_vServer.zVectorField;
            }
            
        }
        
        m_vServer.scale=params.find<std::string>("scale","none");
        if(params.param_present("scale"))
            viewParameter["scale"]="yes";
            

        
        m_vServer.vRendering=params.find<std::string>("vrendering","none");
        if(params.param_present("vrendering"))
            viewParameter["vrendering"]="yes";
            
    
        
        m_vServer.slice=params.find<std::string>("slice","none");
        if(m_vServer.slice=="yes")
        {
            m_vServer.vRendering="none";
            viewParameter["slice"]="yes";
        }
        
        m_vServer.isosurface=params.find<std::string>("isosurface","none");
        if(m_vServer.isosurface=="yes")
        {
            viewParameter["isosurface"]=m_vServer.isosurface;
            m_vServer.vRendering="none";
        }
        tmp="none";
        tmp=params.find<std::string>("shadow","no");
        if(tmp=="yes")
        {
            m_vServer.vShadow=true;
            viewParameter["shadow"]="yes";
        }
        
        m_vServer.vRenderingField=params.find<std::string>("vrenderingfield","none");
        if(params.param_present("vrenderingfield"))
            viewParameter["vrenderingfield"]=m_vServer.vRenderingField;
        
        m_vServer.sliceField=params.find<std::string>("slicefield","none");
        if(params.param_present("slicefield"))
            viewParameter["slicefield"]=m_vServer.sliceField;
        
        m_vServer.isosurfaceField=params.find<std::string>("isosurfacefield","none");
        if(params.param_present("isosurfacefield"))
            viewParameter["isosurfacefield"]=m_vServer.isosurfaceField;
        
        
        m_vServer.isosurfaceValue=params.find<double>("isosurfacevalue",127.5);
        if(params.param_present("isosurfacevalue"))
        {
            std::stringstream ss;//create a stringstream
            ss << m_vServer.isosurfaceValue;//add number to the stream

            viewParameter["isosurfacevalue"]=ss.str();
        }
        
        m_vServer.slicePlane=params.find<std::string>("sliceplane","none");
        if(params.param_present("sliceplane"))
            viewParameter["sliceplane"]=m_vServer.slicePlane;
            
        
        
        
        m_vServer.slicePosition=params.find<int>("sliceposition",0);
        if(params.param_present("sliceposition"))
        {
            std::stringstream ss;//create a stringstream
            ss << m_vServer.slicePosition;//add number to the stream

            viewParameter["sliceposition"]= ss.str();
        }
        
        m_vServer.slicePlaneNormal=params.find<std::string>("sliceplanenormal","none");
        if ( m_vServer.slicePlaneNormal != "none")
        {
            viewParameter["sliceplanenormal"]=m_vServer.slicePlaneNormal;
        }
        
        
        
        m_vServer.slicePlanePoint=params.find<std::string>("sliceplanepoint","none");
        if ( m_vServer.slicePlanePoint != "none")
        {
            viewParameter["sliceplanepoint"]=m_vServer.slicePlanePoint;
        }
        
        
        m_vServer.azimuth=params.find<double>("camazim",0.0);
        if (m_vServer.azimuth != 0.0)
        {
            std::stringstream ss;//create a stringstream
            ss << m_vServer.azimuth;//add number to the stream
            
            viewParameter["camazim"]=ss.str();
            
        }
        
        m_vServer.fovIsGiven=params.param_present("camfov");
        if (m_vServer.fovIsGiven) {
            m_vServer.fov=params.find<double>("camfov",45.0);
            
            std::stringstream ss;//create a stringstream
            ss << m_vServer.fov;//add number to the stream

            viewParameter["camfov"]=ss.str();
        }
        
        m_vServer.elevation=params.find<double>("camelev",0.0);
        if (m_vServer.elevation != 0.0)
        {
            std::stringstream ss;//create a stringstream
            ss << m_vServer.elevation;//add number to the stream

            viewParameter["camelev"]= ss.str();
            
        }
        
        m_vServer.zoom=params.find<double>("zoom",1.0);
        if (m_vServer.zoom != 0.0)
        {
            std::stringstream ss;//create a stringstream
            ss << m_vServer.zoom;//add number to the stream

            viewParameter["zoom"]=ss.str();
            
        }
        
        
        m_vServer.setCameraPos=params.param_present("campos");
        m_vServer.setCameraFocalPoint=params.param_present("camfp");
        m_vServer.setCameraRoll=params.param_present("camroll");
        if(m_vServer.setCameraPos)
        {
            std::string cm;
            cm=params.find<std::string>("campos","0 0 0");
            viewParameter["campos"]=cm;
            
            
            std::stringstream sscm(cm);
            int countcm=0;
            while(!sscm.eof())
            {
                sscm>>m_vServer.cameraPos[countcm];
                countcm++;
                if(countcm==3) break;
            }
        }
        
        if(m_vServer.setCameraFocalPoint)
        {
            std::string cm;
            cm=params.find<std::string>("camfp","0 0 0");
            viewParameter["camfp"]=cm;
            
            std::stringstream sscm(cm);
            int countcm=0;
            while(!sscm.eof())
            {
                sscm>>m_vServer.cameraFocalPoint[countcm];
                countcm++;
                if(countcm==3) break;
            }
        }
        if(m_vServer.setCameraRoll)
        {
            m_vServer.cameraRoll=params.find<double>("camroll",0);
            
            std::stringstream ss;//create a stringstream
            ss << m_vServer.cameraRoll;//add number to the stream

            viewParameter["camroll"]=ss.str();
            
        }
        m_vServer.noDefault=params.find<std::string>("nodefault","none");
        viewParameter["nodefault"]="yes";
        
        m_vServer.noDefault=trim(m_vServer.noDefault);
        if(m_vServer.noDefault!="yes")m_vServer.noDefault="none";
        
        if(m_vServer.slice=="yes") m_vServer.noDefault="none";
        
        tmp="none";
        tmp=params.find<std::string>("largeimage","no");
        if(tmp=="yes"){
            m_vServer.dataRead=false;
            viewParameter["largeimage"]="yes";
        }
        
        m_vServer.color=params.find<std::string>("color","none");
        if (m_vServer.color!= "none")
        {
            viewParameter["color"]="yes";
        }
        
        m_vServer.colorScalar=params.find<std::string>("colorscalar","none");
        if (m_vServer.colorScalar!= "none")
        {
            viewParameter["colorscalar"]=m_vServer.colorScalar;
        }
        
        m_vServer.colorTable=params.find<std::string>("colortable","none");
        if (m_vServer.colorTable!= "none")
        {
            viewParameter["colortable"]=m_vServer.colorTable;
        }
        
        
        m_vServer.opacity=params.find<double>("opacity",0.66);
        if (m_vServer.opacity!=0.66)
        {
            std::stringstream ss;//create a stringstream
            ss << m_vServer.opacity;//add number to the stream

            viewParameter["opacity"]=ss.str();
        }
        m_vServer.uselogscale=params.find<std::string>("logscale","none");
        if (m_vServer.uselogscale== "yes")
        {
            viewParameter["logscale"]="yes";
        }
        
        m_vServer.glyphs=params.find<std::string>("glyphs","pixel");
        if(params.param_present("glyphs"))
            viewParameter["glyphs"]=m_vServer.glyphs;
        
        
        m_vServer.scaleGlyphs=params.find<std::string>("scaleglyphs","none");

        if (m_vServer.scaleGlyphs != "yes")
        {
            m_vServer.scaleGlyphs= "none";
        } else
            viewParameter["scaleglyphs"]="yes";
        
        m_vServer.radius= params.find<double>("radius",1.0);
        if(params.param_present("radius"))
        {
            std::stringstream ss;//create a stringstream
            ss << m_vServer.radius;//add number to the stream

            viewParameter["radius"]= ss.str();
        }
        
        m_vServer.radiusscalar=params.find<std::string>("radiusscalar","none");
        if(params.param_present("radiusscalar"))
            viewParameter["radiusscalar"]=m_vServer.radiusscalar;
        
        m_vServer.height= params.find<double>("height",1.0);
        if(params.param_present("height"))
        {
            std::stringstream ss;//create a stringstream
            ss << m_vServer.height;//add number to the stream

            viewParameter["height"]= ss.str();
        }
        
        m_vServer.heightscalar=params.find<std::string>("heightscalar","none");
        if(params.param_present("heightscalar"))
            viewParameter["heightscalar"]=m_vServer.heightscalar;
        
        m_vServer.vectorScalingFactor=params.find<double>("vectorscalefactor",INVALID_DREAD);
        if(params.param_present("vectorscalefactor"))
        {
            std::stringstream ss;//create a stringstream
            ss << m_vServer.vectorScalingFactor;//add number to the stream

            viewParameter["vectorscalefactor"]=m_vServer.vectorScalingFactor;
        }
        m_vServer.vectorScale=params.find<int>("vectorscale",-1);
        if(params.param_present("vectorscale"))
        {
            std::stringstream ss;//create a stringstream
            ss << m_vServer.vectorScale;//add number to the stream
            
            viewParameter["vectorscale"]= ss.str();

        }
        
        m_vServer.vtkImage=params.find<std::string>("vtkimage","none");
        if(m_vServer.vtkImage=="yes")
        {
            m_vServer.vRendering="none";
            m_vServer.isosurface="none";
            m_vServer.slice="none";
            m_vServer.splotch=false;
        }
        m_vServer.vtkImSize=params.find<std::string>("vtkimsize","128 128 128");
        m_vServer.vtkGausExp=params.find<float>("vtkgausexp",-1);
        m_vServer.vtkScale=params.find<float>("vtkscale",1);
        m_vServer.vtkSpacingFact=params.find<float>("vtkspacingfact",1);
        m_vServer.vtkEcc=params.find<float>("vtkecc",1);
        
        
        tmp =params.find<std::string>("history","no");
        if(params.param_present("history"))
        {
            m_vServer.m_historyEnabled=true;
            m_vServer.m_historyFile ="hist.xml";
            
        }
        
        tmp =params.find<std::string>("historyfile","none");
        if(params.param_present("historyfile"))
        {
            m_vServer.m_historyEnabled=true;
            m_vServer.m_historyFile=tmp;
            
        }
        
        
        //END HISTORY
        
        
        
        m_vServer.mode=params.find<std::string>("mode","none");
        if(params.param_present("mode"))
            viewParameter["mode"]=m_vServer.mode;

        m_vServer.scenario=params.find<std::string>("scenario","etna");
        if(params.param_present("scenario"))
            viewParameter["scenario"]=m_vServer.scenario;

        
        std::string strOpacityTF;
        strOpacityTF=params.find<std::string>("opacityTF","5 3 2.5");
        if(params.param_present("opacityTF"))
            viewParameter["opacityTF"]=strOpacityTF;

        
        std::stringstream ssOpacityTT(strOpacityTF);
        int countTF=0;
        while(!ssOpacityTT.eof())
        {
            ssOpacityTT>>m_vServer.opacityTF[countTF];
            countTF++;
            if(countTF==3) break;
        }
        tmp="none";
        tmp=params.find<std::string>("cliplarge","none");
        if(tmp=="yes")
        {
            m_vServer.cliprange[0]=0.0;
            m_vServer.cliprange[1]=1.0e+13;
            m_vServer.clipset=true;
            viewParameter["cliplarge"]="yes";

        }
        tmp="none";
        m_vServer.clipset=params.param_present("cliprange");
        tmp=params.find<std::string>("cliprange","0.0 1.0e+3");
        
        if(params.param_present("cliprange"))
            viewParameter["cliprange"]=tmp;

        std::stringstream sstmp(tmp);
        int countcr=0;
        while(!sstmp.eof())
        {
            sstmp>>m_vServer.cliprange[countcr];
            countcr++;
            if(countcr==2) break;
        }
        
        
        
    }
    else //if(arguments.size()==1)
    {
        
        for (i=0; i<arguments.size(); i++)
        {
            
            if(arguments[i]=="-x" || arguments[i]=="--x")
            {
                std::string ckInput=arguments[i+1];
                if(ckInput.find_first_of('-')==0)
                {
                    std::cerr<<"Error on --x argument: "<<ckInput<<std::endl;
                    return -1;
                }
                m_vServer.xField=arguments[++i];
                viewParameter["x"]=m_vServer.xField;
                
            }
            
            else if(arguments[i]=="-y"|| arguments[i]=="--y")
            {
                std::string ckInput=arguments[i+1];
                if(ckInput.find_first_of('-')==0)
                {
                    std::cerr<<"Error on --y argument: "<<ckInput<<std::endl;
                    return -1;
                }
                m_vServer.yField=arguments[++i];
                viewParameter["y"]=m_vServer.yField;
                
            }
            
            else if(arguments[i]=="-z"|| arguments[i]=="--z")
            {
                std::string ckInput=arguments[i+1];
                if(ckInput.find_first_of('-')==0)
                {
                    std::cerr<<"Error on --z argument: "<<ckInput<<std::endl;
                    return -1;
                }
                m_vServer.zField=arguments[++i];
                viewParameter["z"]=m_vServer.zField;
                
            }
            
            else if(arguments[i]=="-vx"|| arguments[i]=="--vx")
            {
                std::string ckInput=arguments[i+1];
                if(ckInput.find_first_of('-')==0)
                {
                    std::cerr<<"Error on --vx argument: "<<ckInput<<std::endl;
                    return -1;
                }
                m_vServer.xVectorField=arguments[++i];
                viewParameter["vx"]=m_vServer.xVectorField;
                
            }
            
            else if(arguments[i]=="-vy"|| arguments[i]=="--vy")
            {
                std::string ckInput=arguments[i+1];
                if(ckInput.find_first_of('-')==0)
                {
                    std::cerr<<"Error on --vy argument: "<<ckInput<<std::endl;
                    return -1;
                }
                
                m_vServer.yVectorField=arguments[++i];
                viewParameter["vy"]=m_vServer.yVectorField;
                
            }
            
            else if(arguments[i]=="-vz"|| arguments[i]=="--vz")
            {
                std::string ckInput=arguments[i+1];
                if(ckInput.find_first_of('-')==0)
                {
                    std::cerr<<"Error on --vz argument: "<<ckInput<<std::endl;
                    return -1;
                }
                m_vServer.zVectorField=arguments[++i];
                viewParameter["vz"]=m_vServer.zVectorField;
                
            }
            
            else if(arguments[i]=="--scale")
            {
                m_vServer.scale="yes";
                viewParameter["scale"]="yes";
                
                // std::clog<<"scale="<<m_vServer.scale<<std::endl;
            }
            
            else if(arguments[i]=="--vector")
            {
                m_vServer.vector="yes";
                viewParameter["vector"]="yes";
                
            }
            
            else if(arguments[i]=="--volume")
            {
                m_vServer.volume="yes";
                m_vServer.vRendering="yes";
                // std::clog<<"volume="<<m_vServer.volume<<std::endl;
                viewParameter["volume"]="yes";
                
            }
            else if(arguments[i]=="--splotch")
            {
                std::string ckInput=arguments[i+1];
                if(ckInput.find_first_of('-')==0)
                {
                    std::cerr<<"Error on --splotch argument: "<<ckInput<<std::endl;
                    return -1;
                }
                
                m_vServer.splotch=true;
                m_vServer.vRendering="none";
                m_vServer.splotchpar=arguments[++i];
                viewParameter["splotch"]=m_vServer.splotchpar;
                
                
            }
            
            else if(arguments[i]=="--vrendering")
            {
                m_vServer.vRendering="yes";
                viewParameter["vrendering"]="yes";
                
            }
            else if(arguments[i]=="--shadow")
            {
                m_vServer.vShadow=true;
                viewParameter["shadow"]="yes";
                
            }
            
            else if(arguments[i]=="--vrenderingfield")
            {
                std::string ckInput=arguments[i+1];
                if(ckInput.find_first_of('-')==0)
                {
                    std::cerr<<"Error on --vrenderingfield argument: "<<ckInput<<std::endl;
                    return -1;
                }
                m_vServer.vRenderingField=arguments[++i];
                viewParameter["vrenderingfield"]=m_vServer.vRenderingField;
                
            }
            
            else if(arguments[i]=="--slice")
            {
                m_vServer.slice="yes";
                m_vServer.vRendering="none";
                viewParameter["slice"]="yes";
                
            }
            
            else if(arguments[i]=="--slicefield")
            {
                std::string ckInput=arguments[i+1];
                if(ckInput.find_first_of('-')==0)
                {
                    std::cerr<<"Error on --slicefield argument: "<<ckInput<<std::endl;
                    return -1;
                }
                m_vServer.sliceField=arguments[++i];
                viewParameter["slicefield"]=m_vServer.sliceField;
                
            }
            
            else if(arguments[i]=="--sliceplane")
            {
                std::string ckInput=arguments[i+1];
                if(ckInput.find_first_of('-')==0)
                {
                    std::cerr<<"Error on --sliceplane argument: "<<ckInput<<std::endl;
                    return -1;
                }
                m_vServer.slicePlane=arguments[++i];
                viewParameter["sliceplane"]=m_vServer.slicePlane;
                
            }
            
            else if (arguments[i]=="--sliceposition")
            {
                std::string ckInput=arguments[i+1];
                if(ckInput.substr(0,2).compare("--")==0)
                {
                    std::cerr<<"Error on --sliceposition argument: "<<ckInput<<std::endl;
                    return -1;
                }
                std::stringstream ss1;
                
                ss1<<arguments[++i];
                ss1>>m_vServer.slicePosition;
                viewParameter["sliceposition"]=m_vServer.slicePosition;
                
            }
            else if (arguments[i]=="--sliceplanenormal")
            {
                std::string ck1Input=arguments[i+1];
                std::string ck2Input=arguments[i+2];
                std::string ck3Input=arguments[i+3];
                if(ck1Input.substr(0,2).compare("--")==0 || ck2Input.substr(0,2).compare("--")==0 ||ck3Input.substr(0,2).compare("--")==0)
                {
                    std::cerr<<"Error on --sliceplanenormal argument: "<<ck1Input<<" "<<ck2Input <<" "<< ck3Input<<std::endl;
                    return -1;
                }
                std::stringstream total;
                total<<arguments[i+1]<<" "<<arguments[i+2]<<" "<<arguments[i+3];
                m_vServer.slicePlaneNormal=total.str();
                i=i+3;
                
                viewParameter["sliceplanenormal"]=m_vServer.slicePlaneNormal;
                
            }
            else if (arguments[i]=="--sliceplanepoint")
            {
                std::string ck1Input=arguments[i+1];
                std::string ck2Input=arguments[i+2];
                std::string ck3Input=arguments[i+3];
                if(ck1Input.substr(0,2).compare("--")==0 || ck2Input.substr(0,2).compare("--")==0 ||ck3Input.substr(0,2).compare("--")==0)
                {
                    std::cerr<<"Error on --sliceplanepoint argument: "<<ck1Input<<" "<<ck2Input <<" "<< ck3Input<<std::endl;
                    return -1;
                }
                std::stringstream total;
                total<<arguments[i+1]<<" "<<arguments[i+2]<<" "<<arguments[i+3];
                m_vServer.slicePlanePoint=total.str();
                i=i+3;
                viewParameter["sliceplanepoint"]=m_vServer.slicePlanePoint;
                
            }
            
            else if(arguments[i]=="--isosurface")
            {
                m_vServer.isosurface="yes";
                m_vServer.vRendering="none";
                viewParameter["isosurface"]=m_vServer.isosurface;
                
            }
            
            else if(arguments[i]=="--isosurfacefield")
            {
                std::string ckInput=arguments[i+1];
                if(ckInput.find_first_of('-')==0)
                {
                    std::cerr<<"Error on --isosurfacefield argument: "<<ckInput<<std::endl;
                    return -1;
                }
                
                m_vServer.isosurfaceField=arguments[++i];
                viewParameter["isosurfacefield"]=m_vServer.isosurfaceField;
                
            }
            
            else if(arguments[i]=="--onecolor")
            {
                std::string ckInput=arguments[i+1];
                if(ckInput.find_first_of('-')==0)
                {
                    std::cerr<<"Error on --onecolor argument: "<<ckInput<< "Ignored"<<std::endl;
                    continue;
                }
                m_vServer.oneColor=arguments[++i];
                viewParameter["onecolor"]=m_vServer.oneColor;
                
            }
            else if(arguments[i]=="--showbox")
            {
                m_vServer.showBox=true;
                viewParameter["showbox"]="yes";
                
            }
            else if(arguments[i]=="--showlut")
            {
                m_vServer.showLut=true;
                viewParameter["showlut"]="yes";
                
            }
            else if(arguments[i]=="--showaxes")
            {
                m_vServer.showAxes=true;
                viewParameter["showaxes"]="yes";
                
            }
            else if(arguments[i]=="--wireframe")
            {
                m_vServer.wireframe=true;
                viewParameter["wireframe"]="yes";
                
            }
            else if(arguments[i]=="--vectorline")
            {
                m_vServer.vectorLine=true;
                viewParameter["vectorline"]="yes";
                
            }
            else if(arguments[i]=="--backcolor")
            {
                std::string ckInput=arguments[i+1];
                if(ckInput.find_first_of('-')==0)
                {
                    std::cerr<<"Error on --backcolor argument: "<<ckInput<<" Ignored."<<std::endl;
                    continue;
                }
                m_vServer.backColor=arguments[++i];
                viewParameter["backcolor"]=m_vServer.backColor;
                
            }
            else if(arguments[i]=="--imagesize")
            {
                std::string ckInput=arguments[i+1];
                if(ckInput.find_first_of('-')==0)
                {
                    std::cerr<<"Error on --imagesize argument: "<<ckInput<<" Ignored."<<std::endl;
                    continue;
                }
                m_vServer.imageSize=arguments[++i];
                viewParameter["imagesize"]=m_vServer.backColor;
                
            }
            else if(arguments[i]=="--stereo")
            {
                std::string ckInput=arguments[i+1];
                if(ckInput.find_first_of('-')==0)
                {
                    std::cerr<<"Error on --stereo argument: "<<ckInput<<" Ignored."<<std::endl;
                    continue;
                }
                m_vServer.stereoMode=arguments[++i];
                if(m_vServer.stereoMode == "CrystalEyes" ||  m_vServer.stereoMode =="RedBlue" ||  m_vServer.stereoMode =="Anaglyph"){
                    m_vServer.stereo=true;
                    viewParameter["stereo"]=m_vServer.stereoMode;
                    
                }
                else
                    std::cerr<<"Invalid stereo option "<<m_vServer.stereoMode<<std::endl;
            }
            
            else if (arguments[i]=="--anaglyphmask")
            {
                std::string ck1Input=arguments[i+1];
                std::string ck2Input=arguments[i+2];
                
                if(ck1Input.find_first_of('-')==0 || ck2Input.find_first_of('-')==0)
                {
                    std::cerr<<"Error on --anaglyphmask argument: "<<ck1Input<<" "<<ck2Input<<std::endl;
                    return -1;
                }
                std::stringstream total;
                total<<arguments[i+1]<<" "<<arguments[i+2];
                m_vServer.anaglyphmask=total.str();
                i=i+2;
                viewParameter["anaglyphmask"]=m_vServer.anaglyphmask;
                
            }
            
            else if (arguments[i]=="--anaglyphsat")
            {
                std::string ckInput=arguments[i+1];
                if(ckInput.find_first_of('-')==0)
                {
                    std::cerr<<"Error on --anaglyphsat argument: "<<ckInput<<std::endl;
                    return -1;
                }
                std::stringstream ss1;
                
                ss1<<arguments[++i];
                ss1>>m_vServer.anaglyphsat;
                
                viewParameter["anaglyphsat"]=m_vServer.anaglyphsat;
                
            }
            
            
            else if(arguments[i]=="--isosmooth")
            {
                std::string ckInput=arguments[i+1];
                if(ckInput.find_first_of('-')==0)
                {
                    std::cerr<<"Error on --isosmooth argument: "<<ckInput<<" Ignored."<<std::endl;
                    continue;
                }
                m_vServer.isoSmooth=arguments[++i];
                viewParameter["isosmooth"]=m_vServer.isoSmooth;
                
            }
            else if(arguments[i]=="--isosurfacevalue")
            {
                std::string ckInput=arguments[i+1];
                if(ckInput.find_first_of('-')==0)
                {
                    std::cerr<<"Error on --isosurfacevalue argument: "<<ckInput<<std::endl;
                    return -1;
                }
                
                std::stringstream ss1;
                
                ss1<<arguments[++i];
                ss1>>m_vServer.isosurfaceValue;
                viewParameter["isosurfacevalue"]=m_vServer.isosurfaceValue;
                
            }
            
            else if(arguments[i]=="--color")
            {
                m_vServer.color="yes";
                viewParameter["color"]="yes";
                
            }
            
            else if(arguments[i]=="--colorrangeto")
            {
                std::string ckInput=arguments[i+1];
                if(ckInput.substr(0,2).compare("--")==0)
                {
                    std::cerr<<"Error on --colorrangeto argument: "<<ckInput<<" Ignored."<<std::endl;
                    continue;
                }
                
                std::stringstream ss1;
                m_vServer.isColorRangeTo=true;
                
                ss1<<arguments[++i];
                ss1>>m_vServer.colorRangeTo;
                viewParameter["colorrangeto"]=m_vServer.colorRangeTo;
                
            }
            else if(arguments[i]=="--colorrangefrom")
            {
                std::string ckInput=arguments[i+1];
                
                if(ckInput.substr(0,2).compare("--")==0)
                {
                    std::cerr<<"Error on --colorrangefrom argument: "<<ckInput<<" Ignored."<<std::endl;
                    continue;
                }
                
                std::stringstream ss1;
                m_vServer.isColorRangeFrom=true;
                ss1<<arguments[++i];
                ss1>>m_vServer.colorRangeFrom;
                viewParameter["colorrangefrom"]=m_vServer.colorRangeFrom;
                
            }
            
            else if (arguments[i]=="--out")
            {
                std::string ckInput=arguments[i+1];
                if(ckInput.find_first_of('-')==0)
                {
                    std::cerr<<"Error on --out argument: "<<ckInput<<" Ignored."<<std::endl;
                    continue;
                }
                m_vServer.imageName=arguments[++i];
                viewParameter["out"]=m_vServer.imageName;
                
            }
            
            else if (arguments[i]=="--lfnout")
            {
                std::string ckInput=arguments[i+1];
                if(ckInput.find_first_of('-')==0)
                {
                    std::cerr<<"Error on --lfnout argument: "<<ckInput<<" Ignored."<<std::endl;
                    continue;
                }
                m_vServer.lfnout=arguments[++i];
                viewParameter["lfnout"]=m_vServer.lfnout;
                
            }
            
            else if (arguments[i]=="--VO")
            {
                std::string ckInput=arguments[i+1];
                if(ckInput.find_first_of('-')==0)
                {
                    std::cerr<<"Error on --VO argument: "<<ckInput<<" Ignored."<<std::endl;
                    continue;
                }
                m_vServer.VO=arguments[++i];
                viewParameter["VO"]=m_vServer.VO;
                
            }
            else if (arguments[i]=="--se")
            {
                std::string ckInput=arguments[i+1];
                if(ckInput.find_first_of('-')==0)
                {
                    std::cerr<<"Error on --se argument: "<<ckInput<<" Ignored."<<std::endl;
                    continue;
                }
                m_vServer.se=arguments[++i];
                viewParameter["se"]=m_vServer.se;
                
            }
            
            
            else if (arguments[i]=="--cycle")
            {
                std::string ckInput=arguments[i+1];
                if(ckInput.find_first_of('-')==0)
                {
                    std::cerr<<"Error on --cycle argument: "<<ckInput<<std::endl;
                    return -1;
                }
                m_vServer.cycle=true;
                m_vServer.cycleFile=arguments[++i];
                viewParameter["cycle"]=m_vServer.cycleFile;
                
            }
            else if (arguments[i]=="--cycleoffset")
            {
                std::string ckInput=arguments[i+1];
                if(ckInput.find_first_of('-')==0)
                {
                    std::cerr<<"Error on --cycleoffset argument: "<<ckInput<<std::endl;
                    return -1;
                }
                
                std::stringstream ss1;
                
                ss1<<arguments[++i];
                ss1>>m_vServer.cycleOffset;
                viewParameter["cycleoffset"]=m_vServer.cycleOffset;
                
            }
            else if (arguments[i]=="--cycle_skip_from")
            {
                std::string ckInput=arguments[i+1];
                if(ckInput.find_first_of('-')==0)
                {
                    std::cerr<<"Error on --cycle_skip_from argument: "<<ckInput<<std::endl;
                    return -1;
                }
                std::stringstream ss1;
                
                ss1<<arguments[++i];
                ss1>>m_vServer.cycleSkipFrom;
                viewParameter["cycle_skip_from"]=m_vServer.cycleSkipFrom;
                
            }
            else if (arguments[i]=="--cycle_skip_to")
            {
                std::string ckInput=arguments[i+1];
                if(ckInput.find_first_of('-')==0)
                {
                    std::cerr<<"Error on --cycle_skip_to argument: "<<ckInput<<std::endl;
                    return -1;
                }
                
                std::stringstream ss1;
                
                ss1<<arguments[++i];
                ss1>>m_vServer.cycleSkipTo;
                viewParameter["cycle_skip_to"]=m_vServer.cycleSkipTo;
                
                
            }
            
            else if (arguments[i]=="--camfov")
            {
                std::string ckInput=arguments[i+1];
                
                std::stringstream ss1;
                
                ss1<<arguments[++i];
                ss1>>m_vServer.fov;
                m_vServer.fovIsGiven=true;
                viewParameter["camfov"]=m_vServer.fov;
                
            }
            
            
            else if (arguments[i]=="--camazim")
            {
                std::string ckInput=arguments[i+1];
                if(ckInput.substr(0,2).compare("--")==0)
                {
                    std::clog<<"Error on --camazim argument: "<<ckInput<<std::endl;
                    return -1;
                }
                
                std::stringstream ss1;
                
                ss1<<arguments[++i];
                ss1>>m_vServer.azimuth;
                viewParameter["camazim"]=m_vServer.azimuth;
                
            }
            
            else if (arguments[i]=="--camelev")
            {
                std::string ckInput=arguments[i+1];
                if(ckInput.substr(0,2).compare("--")==0)
                {
                    std::clog<<"Error on --camelev argument: "<<ckInput<<std::endl;
                    return -1;
                }
                
                std::stringstream ss2;
                
                ss2<<arguments[++i];
                ss2>>m_vServer.elevation;
                viewParameter["camelev"]=m_vServer.elevation;
                
                
            }
            
            else if (arguments[i]=="--zoom")
            {
                std::string ckInput=arguments[i+1];
                if(ckInput.find_first_of('-')==0)
                {
                    std::cerr<<"Error on --zoom argument: "<<ckInput<<std::endl;
                    return -1;
                }
                std::stringstream ss3;
                
                ss3<<arguments[++i];
                ss3>>m_vServer.zoom;
                viewParameter["zoom"]=m_vServer.zoom;
                
                
            }
            
            else if (arguments[i]=="--campos")
            {
                std::string ck1Input=arguments[i+1];
                std::string ck2Input=arguments[i+2];
                std::string ck3Input=arguments[i+3];
                
                if(ck1Input.substr(0,2).compare("--")==0 || ck2Input.substr(0,2).compare("--")==0
                   ||ck3Input.substr(0,2).compare("--")==0)
                {
                    std::cerr<<"Error on --campos argument: "<<ck1Input<<" "<<ck2Input <<" "<< ck3Input<<std::endl;
                    return -1;
                }
                
                std::string camposStr;
                
                std::stringstream total;
                total<<ck1Input;
                camposStr.append(ck1Input);
                camposStr.append(" ");

                total>> m_vServer.cameraPos[0];

             //   total>>m_vServer.cameraPos[0];
               // total.str("");

                //total.clear();
                total.str("");
                
                total<<ck2Input;
                
                camposStr.append(ck2Input);
                camposStr.append(" ");

                total>> m_vServer.cameraPos[1];

                //total>>m_vServer.cameraPos[1];
                //total.clear();
                total.str("");

                
                total<<ck3Input;
                camposStr.append(ck3Input);
                total>>m_vServer.cameraPos[2];

                //total<<arguments[i+3];
               // camposStr.append(total.str());
                //total>>m_vServer.cameraPos[2];
                
                i=i+3;
                m_vServer.setCameraPos=true;
                
                viewParameter["campos"]=camposStr;
                
            }
            else if (arguments[i]=="--camfp")
            {
                std::string ck1Input=arguments[i+1];
                std::string ck2Input=arguments[i+2];
                std::string ck3Input=arguments[i+3];
                if(ck1Input.substr(0,2).compare("--")==0 || ck2Input.substr(0,2).compare("--")==0 ||ck3Input.substr(0,2).compare("--")==0)
                {
                    std::cerr<<"Error on --camfp argument: "<<ck1Input<<" "<<ck2Input <<" "<< ck3Input<<std::endl;
                    return -1;
                }
                std::string camfpStr;
                
                std::stringstream total;
                total<<arguments[i+1];
                camfpStr.append(total.str());
                camfpStr.append(" ");
                
                total>>m_vServer.cameraFocalPoint[0];
                total.clear();
                total<<arguments[i+2];
                camfpStr.append(total.str());
                camfpStr.append(" ");
                
                total>>m_vServer.cameraFocalPoint[1];
                total.clear();
                total<<arguments[i+3];
                camfpStr.append(total.str());
                total>>m_vServer.cameraFocalPoint[2];
                i=i+3;
                m_vServer.setCameraFocalPoint=true;
                viewParameter["camfp"]=camfpStr;
                
                
            }
            else if (arguments[i]=="--camroll")
            {
                std::string ckInput=arguments[i+1];
                if(ckInput.find_first_of('-')==0)
                {
                    std::cerr<<"Error on --camroll argument: "<<ckInput<<std::endl;
                    continue;
                }
                
                std::stringstream ss4;
                
                ss4<<arguments[++i];
                ss4>>m_vServer.cameraRoll;
                m_vServer.setCameraRoll=true;
                viewParameter["camroll"]=m_vServer.cameraRoll;
                
            }
            
            
            else if (arguments[i]=="--opacity")
            {
                std::string ckInput=arguments[i+1];
                if(ckInput.find_first_of('-')==0)
                {
                    std::cerr<<"Error on --opacity argument: "<<ckInput<<std::endl;
                    continue;
                }
                
                std::stringstream ss4;
                
                ss4<<arguments[++i];
                ss4>>m_vServer.opacity;
                viewParameter["opacity"]=m_vServer.opacity;
                
                
            }
            else if(arguments[i]=="--colorscalar")
            {
                std::string ckInput=arguments[i+1];
                if(ckInput.find_first_of('-')==0)
                {
                    std::cerr<<"Error on --colorscalar argument: "<<ckInput<<" Ignored."<<std::endl;
                    continue;
                }
                
                m_vServer.colorScalar=arguments[++i];
                viewParameter["colorscalar"]=m_vServer.colorScalar;
                
            }
            
            else if (arguments[i]=="--colortable")
            {
                std::string ckInput=arguments[i+1];
                if(ckInput.find_first_of('-')==0)
                {
                    std::cerr<<"Error on --colortable argument: "<<ckInput<<" Ignored."<<std::endl;
                    continue;
                }
                
                m_vServer.colorTable=arguments[++i];
                viewParameter["colortable"]=m_vServer.colorTable;
                
            }
            
            else if(arguments[i]=="--logscale")
            {
                m_vServer.uselogscale="yes";
                viewParameter["logscale"]="yes";
                
            }

            else if(arguments[i]=="--autorange")
            {
                m_vServer.autoRange="yes";
                viewParameter["autoRange"]="yes";
                
            }
            
            else if(arguments[i]=="--usememory")
            {
                m_vServer.useMemory=true;
                useInMemory = true;
            }

            else if(arguments[i]=="--autorangemin")
            {
                std::string ckInput=arguments[i+1];
                if(ckInput.find_first_of('-')==0)
                {
                    std::cerr<<"Error on --autorange argument: "<<ckInput<<std::endl;
                    return -1;
                }

                std::stringstream ss1;
                
                ss1<<arguments[++i];
                ss1>>m_vServer.autoRangeMin;
                m_vServer.setAutoRangeMin = true;
                viewParameter["autoRangeMin"]=m_vServer.autoRangeMin;
                
            }
            else if(arguments[i]=="--autorangemax")
            {
                std::string ckInput=arguments[i+1];
                if(ckInput.find_first_of('-')==0)
                {
                    std::cerr<<"Error on --autorange argument: "<<ckInput<<std::endl;
                    return -1;
                }
                std::stringstream ss1;

                ss1<<arguments[++i];
                ss1>>m_vServer.autoRangeMax;
                m_vServer.setAutoRangeMax = true;
                viewParameter["autoRangeMax"]=m_vServer.autoRangeMax;
                
            }
            
            else if (arguments[i]=="--glyphs")
            {
                std::string ckInput=arguments[i+1];
                if(ckInput.find_first_of('-')==0)
                {
                    std::cerr<<"Error on --glyphs argument: "<<ckInput<<std::endl;
                    return -1;
                }
                
                m_vServer.glyphs=arguments[++i];
                viewParameter["glyphs"]=m_vServer.glyphs;
                
            }
            
            else if (arguments[i]=="--scaleglyphs")
            {
                m_vServer.scaleGlyphs="yes";
                viewParameter["scaleglyphs"]="yes";
                
            }
            
            else if(arguments[i]=="--radiusscalar")
            {
                std::string ckInput=arguments[i+1];
                if(ckInput.find_first_of('-')==0)
                {
                    std::cerr<<"Error on --radiusscalar argument: "<<ckInput<<std::endl;
                    return -1;
                }
                
                m_vServer.radiusscalar=arguments[++i];
                viewParameter["radiusscalar"]=m_vServer.radiusscalar;
                
            }
            else if (arguments[i]=="--radius")
            {
                std::string ckInput=arguments[i+1];
                if(ckInput.find_first_of('-')==0)
                {
                    std::cerr<<"Error on --radius argument: "<<ckInput<<std::endl;
                    return -1;
                }
                
                std::stringstream ss1;
                
                ss1<<arguments[++i];
                ss1>>m_vServer.radius;
                viewParameter["radius"]=m_vServer.radius;
                
                // std::clog<<"radius="<<m_vServer.radius<<std::endl;
            }
            
            else if (arguments[i]=="--height")
            {
                std::string ckInput=arguments[i+1];
                if(ckInput.find_first_of('-')==0)
                {
                    std::cerr<<"Error on --height argument: "<<ckInput<<std::endl;
                    return -1;
                }
                
                std::stringstream ss1;
                
                ss1<<arguments[++i];
                ss1>>m_vServer.height;
                viewParameter["height"]=m_vServer.height;
                
            }
            
            else if(arguments[i]=="--heightscalar")
            {
                std::string ckInput=arguments[i+1];
                if(ckInput.find_first_of('-')==0)
                {
                    std::cerr<<"Error on --heightscalar argument: "<<ckInput<<std::endl;
                    return -1;
                }
                m_vServer.heightscalar=arguments[++i];
                viewParameter["heightscalar"]=m_vServer.heightscalar;
                
            }
            else if (arguments[i]=="--vectorscalefactor")
            {
                std::string ckInput=arguments[i+1];
                if(ckInput.rfind('-')==0)
                {
                    std::cerr<<"Error on --vectorscalefactor argument: "<<ckInput<<std::endl;
                    return -1;
                }
                
                std::stringstream ss1;
                
                ss1<<arguments[++i];
                ss1>>m_vServer.vectorScalingFactor;
                viewParameter["vectorscalefactor"]=m_vServer.vectorScalingFactor;
                
            }
            else if (arguments[i]=="--vectorscale")
            {
                std::string ckInput=arguments[i+1];
                if(ckInput.find_first_of('-')==0)
                {
                    std::cerr<<"Error on --vectorscale argument: "<<ckInput<<std::endl;
                    return -1;
                }
                
                std::stringstream ss1;
                
                ss1<<arguments[++i];
                ss1>>m_vServer.vectorScale;
                if(m_vServer.vectorScale<0 ||m_vServer.vectorScale>1)  m_vServer.vectorScale=-1;
                
                viewParameter["vectorscale"]=m_vServer.vectorScale;
                
            }
            
            else if (arguments[i]=="--nodefault")
            {
                m_vServer.noDefault="yes";
                viewParameter["nodefault"]="yes";
                
            }
            
            else if (arguments[i]=="--largeimage")
            {
                m_vServer.dataRead=false;
                viewParameter["largeimage"]="yes";
                
            }
            
            else if (arguments[i]=="--savepar")
            {
                std::cout<<"--savepar obsolete option. No more active."<<std::endl;
            }
            
            else if (arguments[i]=="--loadpar")
            {
                std::cout<<"--loadpar obsolete option. No more active."<<std::endl;
                
            }
            else if (arguments[i]=="--mode")
            {
                std::string ckInput=arguments[i+1];
                if(ckInput.find_first_of('-')==0)
                {
                    std::cerr<<"Error on --mode argument: "<<ckInput<<std::endl;
                    return -1;
                }
                m_vServer.mode=ckInput;
                viewParameter["mode"]=m_vServer.mode;
                
            }
            else if (arguments[i]=="--scenario")
            {
                std::string ckInput=arguments[i+1];
                if(ckInput.find_first_of('-')==0)
                {
                    std::cerr<<"Error on --scenario argument: "<<ckInput<<std::endl;
                    return -1;
                }
                m_vServer.scenario=ckInput;
                viewParameter["scenario"]=m_vServer.scenario;
                
            }
            else if (arguments[i]=="--opacityTF")
            {
                
                std::string opacStr;
                
                std::string ckInput0=arguments[i+1];
                if(ckInput0.find_first_of('-')==0) continue;
                std::stringstream ssOpacityTF0(arguments[++i]);
                opacStr.append(ssOpacityTF0.str());
                opacStr.append(" ");
                ssOpacityTF0>>m_vServer.opacityTF[0];
                std::string ckInput1=arguments[i+1];
                if(ckInput1.find_first_of('-')==0) continue;
                std::stringstream ssOpacityTF1(arguments[++i]);
                opacStr.append(ssOpacityTF1.str());
                opacStr.append(" ");
                ssOpacityTF1>>m_vServer.opacityTF[1];
                std::string ckInput2=arguments[i+1];
                if(ckInput2.find_first_of('-')==0) continue;
                std::stringstream ssOpacityTF2(arguments[++i]);
                opacStr.append(ssOpacityTF2.str());
                ssOpacityTF2>>m_vServer.opacityTF[2];
                
                
                
                viewParameter["opacityTF"]=opacStr;
                
            }
            else if (arguments[i]=="--cliplarge")
            {
                m_vServer.cliprange[0]=0.0;
                m_vServer.cliprange[1]=1.0e+13;
                
                viewParameter["cliplarge"]="yes";
                
            }
            else if (arguments[i]=="--cliprange")
            {
                std::string ck1Input=arguments[i+1];
                std::string ck2Input=arguments[i+2];
                
                if(ck1Input.substr(0,2).compare("--")==0 || ck2Input.substr(0,2).compare("--")==0)
                {
                    std::cerr<<"Error on --cliprange argument: "<<ck1Input<<" "<<ck2Input <<std::endl;
                    return -1;
                }
                std::stringstream total;
                std::string cliprangeStr;
                
                
                total<<arguments[i+1];
                total>>m_vServer.cliprange[0];
                cliprangeStr.append(total.str());
                cliprangeStr.append(" ");
                total.clear();
                total<<arguments[i+2];
                total>>m_vServer.cliprange[1];
                cliprangeStr.append(total.str());
                i=i+2;
                
                viewParameter["cliprange"]=cliprangeStr;
                
            }
            
            else  if (arguments[i]=="--history")
            {
                m_vServer.m_historyEnabled=true;
                m_vServer.m_historyFile ="hist.xml";
                
            }
            else if (arguments[i]=="--historyfile")
            {
                
                
                m_vServer.m_historyEnabled=true;
                std::string ckInput=arguments[i+1];
                if(ckInput.find_first_of('-')==0)
                {
                    std::cerr<<"Error on "<<arguments[i]<< " argument: "<<ckInput<<std::endl;
                    return -1;
                }
                m_vServer.m_historyFile=arguments[++i];
                
            }
            else if (arguments[i]=="--file")
            {
                m_vServer.path=arguments[++i];
                fileIsAtTheEnd=false;
 
            }
            
            
            else
            {
                if(i<arguments.size()-1) std::cerr<<"Invalid parameter "<<arguments[i]<<std::endl;
            }
        }
        if (fileIsAtTheEnd)
            m_vServer.path=arguments[(arguments.size()-1)];
        
    } //else if(arguments.size()==1)
    // // std::cerr<<m_vServer.path<<std::endl;
    if(m_vServer.path.find(".bin") ==std::string::npos)
        m_vServer.path +=".bin";
#ifdef GLITE
    
    if(m_vServer.path.substr(0,6) =="lfn://")
    {
        if(m_vServer.VO=="none")
        {
            std::cerr<<"VO option must be given for a logical filename file."<<std::endl;
            return 1;
        }
        std::map<std::string, std::string > downLoad;
        downLoad.insert(make_pair("file",m_vServer.path));
        downLoad.insert(make_pair("VO",m_vServer.VO));
        downLoad.insert(make_pair("out",m_vServer.imageName));
        gLiteGw gLInterface;
        m_vServer.path=gLInterface.readVF(downLoad,-1);
        m_inputLfnGiven=true;
    }
#endif
    
    if(m_vServer.splotch)  //splotch disactive the scaling
        m_vServer.scale="none";
    if(m_vServer.volume=="yes" && m_vServer.splotch)
    {
        m_vServer.splotch=false;
        std::cout<<"--splotch option cannot be used for volumes."<<std::endl;
    }
    
    if(m_vServer.vRendering=="yes")
    {
        m_vServer.slice="none";
        m_vServer.isosurface="none";
    }
    if(m_vServer.isosurface=="yes")
    {
        m_vServer.slice="none";
    }
    
    if(m_vServer.slice=="yes")
    {
        m_vServer.stereo=false;  //stereo not enable for plane view
        m_vServer.noDefault="yes";// No default image can be produced with slice
    }
    
    ///////////////////////
    ///// TEST FOR INPUT PARAMETERS
    //////////////////////
    
    if(m_vServer.slicePlane!="x" && m_vServer.slicePlane!="y"&& m_vServer.slicePlane!="z"&& m_vServer.slicePlane!="none")
    {
        std::cerr<<"Invalid value for sliceplane option "<<m_vServer.slicePlane<<" Option ignored."<<std::endl;
        m_vServer.slicePlane="none";
    } else{
        if(m_vServer.slicePlane!="none")
        {
            m_vServer.sliceOrthoNormal=true;
            m_vServer.slicePlaneNormal="none";
            m_vServer.slicePlanePoint="none";
        }
    }
    if(m_vServer.slice=="yes" && m_vServer.slicePlane=="none" && m_vServer.slicePlaneNormal=="none" && m_vServer.slicePlanePoint=="none")
    {
        m_vServer.slice="none";
        std::cerr<<"Invalid value for sliceplane option "<<m_vServer.slicePlane;
        std::cerr<<" , sliceplanenormal "<<m_vServer.slicePlaneNormal;
        std::cerr<<" or sliceplanepoint "<<m_vServer.slicePlanePoint;
        std::cerr<<" Option ignored."<<std::endl;
    }
    if(m_vServer.slicePlane=="x")
        m_vServer.nSlicePlane=0;
    if(m_vServer.slicePlane=="y")
        m_vServer.nSlicePlane=1;
    if(m_vServer.slicePlane=="z")
        m_vServer.nSlicePlane=2;
    
    if(m_vServer.vtkImage!="none")
    {
        std::stringstream vtkimtmp(m_vServer.vtkImSize);
        for(int i=0;i<3;i++)
        {
            if(vtkimtmp.eof() )
            {
                std::cerr<<"Warning: Invalid field vtkimsize "<<m_vServer.vtkImSize<<" "<<std::endl;
            }
            vtkimtmp>> m_vServer.vtkImSizeValue[i];
        }
    }
    if(m_vServer.nSlicePlane<0)
    {
        if(m_vServer.slicePlaneNormal!="none" )
        {
            std::stringstream sspntmp(m_vServer.slicePlaneNormal);
            for(int i=0;i<3;i++)
            {
                if(sspntmp.eof() )
                {
                    std::cerr<<"Invalid field sliceplanenormal "<<m_vServer.slicePlaneNormal<<" ";
                    m_vServer.slicePlaneNormal=="none";
                }
                sspntmp>> m_vServer.slicePlaneNormalValue[i];
            }
        }
        if(m_vServer.slicePlanePoint!="none" )
        {
            std::stringstream sspptmp(m_vServer.slicePlanePoint);
            for(int i=0;i<3;i++)
            {
                if(sspptmp.eof() )
                {
                    std::cerr<<"Invalid field sliceplanepoint "<<m_vServer.slicePlanePoint<<std::endl;
                    m_vServer.slicePlanePoint=="none";
                }
                sspptmp>> m_vServer.slicePlanePointValue[i];
            }
        }
    }
    
    if(m_vServer.color=="yes")
        setColorLut();
    setGlyphs();
    
    if(!m_vServer.internalData && !m_vServer.useMemory)
    {
        std::fstream inFile;
        //  std::clog<<m_vServer.path<<std::endl;
        inFile.open(m_vServer.path.c_str());
        
        if(!inFile)
        {
            std::cerr<<"the path is incorrect, please try again  or --help for help"<<std::endl;
            return -1;
        }
        
        else
            inFile.close();
    }
    if(m_vServer.cycle)
    {
        m_vServer.noDefault="yes";
        readCycleFile();
    }
    
    if(m_vServer.volume!="yes" && (m_vServer.xField=="none"||m_vServer.yField=="none" || m_vServer.zField=="none"))
        std::cerr<<"Warning x y and z fields not assigned. Unpredictable behavior could occur"<<std::endl;
    
    if(m_vServer.volume=="yes" && (m_vServer.vRendering!="yes" && m_vServer.slice!="yes" && m_vServer.isosurface !="yes"))
    {
        std::cerr<<"Error: no volume visualisation tecnique is specified"<<std::endl;
        return -1;
    }
    
    if(m_vServer.volume=="yes" && m_vServer.vRendering=="yes" && m_vServer.vRenderingField=="none")
    {
        std::cerr<<"Error: --vrenderingfield field not assigned"<<std::endl;
        return -1;
    }
    if(m_vServer.volume=="yes" && m_vServer.slice=="yes" && m_vServer.sliceField=="none")
    {
        std::cerr<<"Error: --slicefield field not assigned"<<std::endl;
        return -1;
    }
    if(m_vServer.volume=="yes" && m_vServer.isosurface =="yes" && m_vServer.isosurfaceField=="none")
    {
        std::cerr<<"Error: --isosurfacefield field not assigned"<<std::endl;
        return -1;
    }
    if(m_vServer.savepar=="yes")
        saveParameters();
    if(m_vServer.splotch)
        readSplocthColumn(); //read Splotch par file and set m_vServer sploth Field names
    /*  if(m_vServer.loadpar!="yes")
     {*/
    if(m_vServer.isColorRangeTo && m_vServer.isColorRangeFrom)
        if(m_vServer.colorRangeTo<m_vServer.colorRangeFrom)
        {
		std:cerr<<"Invalid parameters colorrangeto and  colorrangefrom "<<std::endl;
            m_vServer.isColorRangeTo=false;
            m_vServer.isColorRangeFrom=false;
        }
    if(!m_vServer.internalData && !m_vServer.useMemory)
    {
        if(!readHeader())
            return -1;
        if(!setAxisScalarVectorAndVolumesFields()) return -1; // check for all fields in the table
    }
    //   }
    setNumberImages();
    
    
    
    return 0;
}


//---------------------------------------------------------------------
bool OptionsSetter::setAxisScalarVectorAndVolumesFields ()
//---------------------------------------------------------------------
{
    if(m_vServer.volume!="yes")
    {
        if(m_vServer.xField!="none")
        {
            bool colAssigned=false;
            for (int i=0; i<m_vServer.nCols; i++)
            {
                if (iCompare(m_vServer.fieldNames[i],m_vServer.xField)==0)
                {
                    m_vServer.x=i;
                    colAssigned=true;
                    //         m_vServer.xField= m_vServer.fieldNames[i];
                    break;
                }
            }
            if(!colAssigned)
            {
                std::cerr<<"Invalid column name in --x option: "<<m_vServer.xField<<std::endl;
                return false;
            }
        }  else
        {
            int nOfField=m_vServer.fieldNames.size();
            if(nOfField>0)
                m_vServer.xField= m_vServer.fieldNames[0];
            std::string::size_type caseA;
            std::string::size_type caseB;
            std::string::size_type caseC;
            std::string::size_type caseD;
            std::string::size_type caseE;
            std::string names;
            for (int i=0; i<m_vServer.nCols; i++)
            {
                names=m_vServer.fieldNames[i].c_str();
                caseD =names.find("X",0);
                caseE =names.find("x",0);
                caseA =names.find("DE",0);
                caseB =names.find("De",0);
                caseC =names.find("de",0);
                if (caseD!=std::string::npos||caseE!=std::string::npos||caseA!=std::string::npos||caseB!=std::string::npos||caseC!=std::string::npos)
                {
                    m_vServer.x=i;
                    m_vServer.xField= m_vServer.fieldNames[i];
                    break;
                }
            }
        } //else if m_vServer.xField!="none"
        if(m_vServer.yField!="none")
        {
            bool colAssigned=false;
            for (int i=0; i<m_vServer.nCols; i++)
            {
                if (iCompare(m_vServer.fieldNames[i],m_vServer.yField)==0)
                {
                    m_vServer.y=i;
                    colAssigned=true;
                    //          m_vServer.yField= m_vServer.fieldNames[i];
                    break;
                }
            }
            if(!colAssigned)
            {
                std::cerr<<"Invalid column name in --y option: "<<m_vServer.yField<<std::endl;
                return false;
            }
        }  else
        {
            int nOfField=m_vServer.fieldNames.size();
            if(nOfField>1)
                m_vServer.yField= m_vServer.fieldNames[1];
            std::string::size_type caseA;
            std::string::size_type caseB;
            std::string::size_type caseC;
            std::string::size_type caseD;
            std::string::size_type caseE;
            std::string names;
            for (int i=0; i<m_vServer.nCols; i++)
            {
                names=m_vServer.fieldNames[i].c_str();
                caseD =names.find("Y",0);
                caseE =names.find("y",0);
                caseA =names.find("RA",0);
                caseB =names.find("Ra",0);
                caseC =names.find("ra",0);
                if (caseD!=std::string::npos||caseE!=std::string::npos||caseA!=std::string::npos||caseB!=std::string::npos||caseC!=std::string::npos)
                {
                    m_vServer.y=i;
                    m_vServer.yField= m_vServer.fieldNames[i];
                    break;
                }
            }
        } //else if m_vServer.yField!="none"
        if(m_vServer.zField!="none")
        {
            bool colAssigned=false;
            for (int i=0; i<m_vServer.nCols; i++)
            {
                if (iCompare(m_vServer.fieldNames[i],m_vServer.zField)==0)
                {
                    m_vServer.z=i;
                    colAssigned=true;
                    //          m_vServer.zField= m_vServer.fieldNames[i];
                    break;
                }
            }
            if(!colAssigned)
            {
                std::cerr<<"Invalid column name in --z option: "<<m_vServer.zField<<std::endl;
                return false;
            }
        }  else
        {
            int nOfField=m_vServer.fieldNames.size();
            if(nOfField>2)
                m_vServer.zField= m_vServer.fieldNames[2];
            std::string::size_type caseA;
            std::string::size_type caseB;
            std::string::size_type caseC;
            std::string::size_type caseD;
            std::string::size_type caseE;
            std::string names;
            for (int i=0; i<m_vServer.nCols; i++)
            {
                names=m_vServer.fieldNames[i].c_str();
                caseD =names.find("Z",0);
                caseE =names.find("z",0);
                caseA =names.find("MAG",0);
                caseB =names.find("Mag",0);
                caseC =names.find("mag",0);
                if (caseD!=std::string::npos||caseE!=std::string::npos||caseA!=std::string::npos||caseB!=std::string::npos||caseC!=std::string::npos)
                {
                    m_vServer.z=i;
                    m_vServer.zField= m_vServer.fieldNames[i];
                    break;
                }
            }
        } //else if m_vServer.zField!="none"
        //points
    }
    
    //vectors
    if(m_vServer.vector=="yes" &&( m_vServer.xVectorField!="none"&& m_vServer.yVectorField!="none" && m_vServer.zVectorField!="none"))
    {
        int vectorAssigned=0;
        for (int i=0; i<m_vServer.nCols; i++)
        {
            
            if (iCompare(m_vServer.fieldNames[i],m_vServer.xVectorField)==0)
            {
                m_vServer.vx=i;
                vectorAssigned++;
                //std::clog << m_vServer.vx << std::endl;
            }
            
            else if (iCompare(m_vServer.fieldNames[i],m_vServer.yVectorField)==0)
            {
                m_vServer.vy=i;
                vectorAssigned++;
                //std::clog <<m_vServer.vy << std::endl;
            }
            else if (iCompare(m_vServer.fieldNames[i],m_vServer.zVectorField)==0)
            {
                m_vServer.vz=i;
                vectorAssigned++;
                //std::clog << m_vServer.vz<< std::endl;
            }
        }
        if(vectorAssigned!=3)
            std::cerr<<"Severe Warning: One or more invalid vector fields vx="<<m_vServer.xVectorField<<" vy="<<m_vServer.yVectorField<<" vz="<<m_vServer.zVectorField<<endl;
    }
    
    //radiusscalar
    if(m_vServer.radiusscalar!="none"&& m_vServer.scaleGlyphs=="yes" && m_vServer.nGlyphs!=0)
    {
        bool colAssigned=false;
        for (int i=0; i<m_vServer.nCols; i++)
        {
            if (iCompare(m_vServer.fieldNames[i],m_vServer.radiusscalar)==0)
            {
                m_vServer.nRadius=i;
                colAssigned=true;
                break;
                //        m_vServer.radiusscalar= m_vServer.fieldNames[i];
                //std::clog << m_vServer.nRadius << std::endl;
            }
        }
        if(!colAssigned)
        {
            std::cerr<<"Warning. Inavild radius scalar name: "<<m_vServer.radiusscalar<<std::endl;
            m_vServer.radiusscalar="none";
        }
    }
    //heightscalar
    if(m_vServer.heightscalar!="none"&& m_vServer.scaleGlyphs!="none"&& m_vServer.nGlyphs!=0 && m_vServer.nGlyphs!=1 )
    {
        bool colAssigned=false;
        for (int i=0; i<m_vServer.nCols; i++)
        {
            if (iCompare(m_vServer.fieldNames[i],m_vServer.heightscalar)==0)
            {
                m_vServer.nHeight=i;
                colAssigned=true;
                break;
                //        m_vServer.heightscalar= m_vServer.fieldNames[i];
                //std::clog << m_vServer.nWeigth << std::endl;
            }
        }
        if(!colAssigned)
        {
            std::cerr<<"Warning. Invalid height scalar name: "<<m_vServer.heightscalar<<std::endl;
            m_vServer.heightscalar="none";
        }
    }
    
    //volumerendering
    if(m_vServer.volume=="yes" &&( m_vServer.vRendering=="yes"&& m_vServer.vRenderingField!="none" ))
    {
        bool colAssigned=false;
        for (int i=0; i<m_vServer.nCols; i++)
        {
            if (iCompare(m_vServer.fieldNames[i],m_vServer.vRenderingField)==0)
            {
                m_vServer.nVRenderingField=i;
                colAssigned=true;
                break;
                //        m_vServer.vRenderingField= m_vServer.fieldNames[i];
                //std::clog << m_vServer.volumeField << std::endl;
            }
        }
        if(!colAssigned)
        {
            std::cerr<<"Error. Invalid volume rendering field name: "<<m_vServer.vRenderingField<<std::endl;
            return false;
        }
    }
    
    //slice
    if(m_vServer.volume=="yes"&&( m_vServer.slice=="yes"&& m_vServer.sliceField!="none"))
    {
        bool colAssigned=false;
        for (int i=0; i<m_vServer.nCols; i++)
        {
            if (iCompare(m_vServer.fieldNames[i],m_vServer.sliceField)==0)
            {
                m_vServer.nSliceField=i;
                colAssigned=true;
                break;
                //        m_vServer.sliceField= m_vServer.fieldNames[i];
                //std::clog << m_vServer.volumeField << std::endl;
            }
        }
        if(!colAssigned)
        {
            std::cerr<<"Error. Invalid  slice field name: "<<m_vServer.sliceField<<std::endl;
            return false;
        }
        
    }
    
    //isosurface
    if(m_vServer.volume=="yes"&&( m_vServer.isosurface=="yes"&& m_vServer.isosurfaceField!="none" ))
    {
        bool colAssigned=false;
        for (int i=0; i<m_vServer.nCols; i++)
        {
            if (iCompare(m_vServer.fieldNames[i],m_vServer.isosurfaceField)==0)
            {
                m_vServer.nIsosurfaceField=i;
                colAssigned=true;
                break;
                //        m_vServer.isosurfaceField= m_vServer.fieldNames[i];
                //std::clog << m_vServer.volumeField << std::endl;
            }
        }
        if(!colAssigned)
        {
            std::cerr<<"Error. Invalid  isosurface field name: "<<m_vServer.isosurfaceField<<std::endl;
            return false;
        }
    }
    //lutscalar
    if(m_vServer.color=="yes" && m_vServer.colorScalar!="none")
    {
        bool colAssigned=false;
        for (int i=0; i<m_vServer.nCols; i++)
        {
            
            if (iCompare(m_vServer.fieldNames[i],m_vServer.colorScalar)==0)
            {
                m_vServer.nColorScalar=i;
                colAssigned=true;
                break;
                //        m_vServer.colorScalar= m_vServer.fieldNames[i];
                //std::clog << m_vServer.nColorScalar << std::endl;
            }
        }
        if(!colAssigned)
            std::cerr<<"Warning. Invalid color field name: "<<m_vServer.isosurfaceField<<std::endl;
        
    }
    if(m_vServer.splotch)
    {
        
        if(m_vServer.hsmlField!="none")
        {
            bool colAssigned=false;
            for (int i=0; i<m_vServer.nCols; i++)
            {
                if (iCompare(m_vServer.fieldNames[i],m_vServer.hsmlField)==0)
                {
                    m_vServer.hsml=i;
                    colAssigned=true;
                    break;
                    //        m_vServer.colorScalar= m_vServer.fieldNames[i];
                    //std::clog << m_vServer.nColorScalar << std::endl;
                }
            }
            if(!colAssigned)
            {
                std::cerr<<"Warning. Invalid Splotch  HSML field name: "<<m_vServer.hsmlField<<std::endl;
                m_vServer.hsmlField="none";
            }
        }
        if(m_vServer.spIntensityField!="none")
        {
            bool colAssigned=false;
            for (int i=0; i<m_vServer.nCols; i++)
            {
                if (iCompare(m_vServer.fieldNames[i],m_vServer.spIntensityField)==0)
                {
                    m_vServer.spIntensity=i;
                    colAssigned=true;
                    break;
                    //        m_vServer.colorScalar= m_vServer.fieldNames[i];
                    //std::clog << m_vServer.nColorScalar << std::endl;
                }
            }
            if(!colAssigned)
            {
                std::cerr<<"Warning. Invalid Splotch Intensity field name: "<<m_vServer.spIntensityField<<std::endl;
                m_vServer.spIntensityField="none";
            }
        }
        if(m_vServer.sprField!="none")
        {
            bool colAssigned=false;
            for (int i=0; i<m_vServer.nCols; i++)
            {
                if (iCompare(m_vServer.fieldNames[i],m_vServer.sprField)==0)
                {
                    m_vServer.spr=i;
                    colAssigned=true;
                    break;
                    //        m_vServer.colorScalar= m_vServer.fieldNames[i];
                    //std::clog << m_vServer.nColorScalar << std::endl;
                }
            }
            if(!colAssigned)
            {
                std::cerr<<"Warning. Invalid Splotch smoothing field name: "<<m_vServer.sprField<<std::endl;
                m_vServer.sprField="none";
            }
        }
        
        if(m_vServer.spIField!="none")
        {
            bool colAssigned=false;
            for (int i=0; i<m_vServer.nCols; i++)
            {
                if (iCompare(m_vServer.fieldNames[i],m_vServer.spIField)==0)
                {
                    m_vServer.spI=i;
                    colAssigned=true;
                    break;
                    //        m_vServer.colorScalar= m_vServer.fieldNames[i];
                    //std::clog << m_vServer.nColorScalar << std::endl;
                }
            }
            if(!colAssigned)
            {
                std::cerr<<"Warning. Invalid Splotch Intensity field name: "<<m_vServer.spIField<<std::endl;
                m_vServer.spIField="none";
            }
        }
        
        if(m_vServer.spC1Field!="none")
        {
            bool colAssigned=false;
            for (int i=0; i<m_vServer.nCols; i++)
            {
                if (iCompare(m_vServer.fieldNames[i],m_vServer.spC1Field)==0)
                {
                    m_vServer.spC1=i;
                    colAssigned=true;
                    break;
                    //        m_vServer.colorScalar= m_vServer.fieldNames[i];
                    //std::clog << m_vServer.nColorScalar << std::endl;
                }
            }
            if(!colAssigned)
            {
                std::cerr<<"Warning. Invalid Splotch Color 1 field name: "<<m_vServer.spC1Field<<std::endl;
                m_vServer.spC1Field="none";
            }
        }
        
        if(m_vServer.spC2Field!="none")
        {
            bool colAssigned=false;
            for (int i=0; i<m_vServer.nCols; i++)
            {
                if (iCompare(m_vServer.fieldNames[i],m_vServer.spC2Field)==0)
                {
                    m_vServer.spC2=i;
                    colAssigned=true;
                    break;
                    //        m_vServer.colorScalar= m_vServer.fieldNames[i];
                    //std::clog << m_vServer.nColorScalar << std::endl;
                }
            }
            if(!colAssigned)
            {
                std::cerr<<"Warning. Invalid Splotch Color 2 field name: "<<m_vServer.spC2Field<<std::endl;
                m_vServer.spC2Field="none";
            }
        }
        
        if(m_vServer.spC3Field!="none")
        {
            bool colAssigned=false;
            for (int i=0; i<m_vServer.nCols; i++)
            {
                if (iCompare(m_vServer.fieldNames[i],m_vServer.spC3Field)==0)
                {
                    m_vServer.spC3=i;
                    colAssigned=true;
                    break;
                    //        m_vServer.colorScalar= m_vServer.fieldNames[i];
                    //std::clog << m_vServer.nColorScalar << std::endl;
                }
            }
            if(!colAssigned)
            {
                std::cerr<<"Warning. Invalid Splotch Color 3 field name: "<<m_vServer.spC3Field<<std::endl;
                m_vServer.spC3Field="none";
            }
        }
        
        if(m_vServer.spColorField!="none")
        {
            bool colAssigned=false;
            for (int i=0; i<m_vServer.nCols; i++)
            {
                if (iCompare(m_vServer.fieldNames[i],m_vServer.spColorField)==0)
                {
                    m_vServer.spColor=i;
                    colAssigned=true;
                    break;
                    //        m_vServer.colorScalar= m_vServer.fieldNames[i];
                    //std::clog << m_vServer.nColorScalar << std::endl;
                }
            }
            if(!colAssigned)
            {
                std::cerr<<"Warning. Invalid Splotch Color field name: "<<m_vServer.spColorField<<std::endl;
                m_vServer.spColorField="none";
            }
        }
        
    }
    return true;
}


//---------------------------------------------------------------------
void OptionsSetter::showHelp()
//---------------------------------------------------------------------
{
    
    std::clog<<std::endl;
    std::clog<<"VisIVOViewer Version 3.5.0 March 19th 2025 "<<std::endl<<std::endl;
    
    std::clog<<"     [pathfile] (mandatory) Absolute path file. Path must be the last command(e.g. /home/user/VisivoBinaryTable.bin)"<<std::endl<<std::endl;
    
    std::clog<<" --out  [filename]   (optional) Change default file name and or directory "<<std::endl<<std::endl;
    
    std::clog<<" --cycle [filename]   (optional) Will produce a sequence of images reading data of azimuth elevation and zooming from the file given with this parameter. This option will prevent the production of default images and the camera position and zoom factor of the line command will be ignored. The file for cycle must contain three ascii values for each row. In the followin the main parameters are listed. Please refer to the user guide for a complete description."<<std::endl<<std::endl;
    
    std::clog<<" --cycleoffset [value] (optional). Is the value for the progressive number of files produced with the cycle option. Default value is 0."<<std::endl<<std::endl;
    
    
    std::clog<<" -x    [field] (optional) Select x field to load  (-x x)"<<std::endl<<std::endl;
    
    std::clog<<" -y    [field]   (optional) Select y field to load  (-y y)"<<std::endl<<std::endl;
    
    std::clog<<" -z    [field]  (optional) Select z field to load  (-z z)"<<std::endl<<std::endl;
    
    std::clog<<" --vector    (optional) Mandatory if you want create a vector (--vector)"<<std::endl<<std::endl;
    
    std::clog<<" --vx    [field] (optional) Select x field to load for create a vector "<<std::endl<<std::endl;
    
    
    std::clog<<" --vy    [field]   (optional) Select y field to load  for create a vector "<<std::endl<<std::endl;
    
    
    std::clog<<" --vz    [field]  (optional) Select z field to load for create a vector "<<std::endl<<std::endl;
    
    
    std::clog<<" --nodefault     (optional)  If the user want only one default image or an image where the user select the camera position and/or the zoom"<<std::endl<<std::endl;
    
    std::clog<<" --stereo (optional) May assume RedBlue or CrystalEyes or Anaglyph. Produce stereoscopic images"<<std::endl<<std::endl;
    
    std::clog<<" --anaglyphsat (optional).  Anaglyph color saturation from 0.0 to 1.0:";
    std::clog<<" 0.0 means that no color from the original object is maintained, ";
    std::clog<<"1.0 means all of the color is maintained. The default value is 0.65.";
    std::clog<<" Too much saturation can produce uncomfortable 3D viewing"<<std::endl<<std::endl;
    
    std::clog<<" --anaglyphmask (optional). These two numbers are bits mask that control";
    std::clog<<" which color channels of the original stereo images are used to produce";
    std::clog<<" the final anaglyph image. The first value is the color mask for the left";
    std::clog<<" view, the second the mask for the right view. If a bit in the mask is";
    std::clog<<" on for a particular color for a view, that color is passed on to the";
    std::clog<<" final view; if it is not set, that channel for that view is ignored.";
    std::clog<<" The bits are arranged as r, g, and b, so r = 4, g = 2, and b = 1.";
    std::clog<<" By default, the first value (the left view) is set to 4, and the second";
    std::clog<<" value is set to 3"<<std::endl<<std::endl;
    
    
    std::clog<<" --color (optional) Mandatory if you want use color table "<<std::endl<<std::endl;
    
    std::clog<<" --colorscalar  [field]  (optional) Select lut field"<<std::endl<<std::endl;
    
    std::clog<< " --colortable  [name]  (optional) Select the table that want for color the glyphs"<<std::endl<<std::endl;
    
    std::clog<< " --logscale  (optional) use this command if you want logaritmic scale. If the field selected for looktable have value that are <=0 this option will be ignored"<<std::endl<<std::endl;
    
    std::clog<< " --camazim   [double]  (optional)Rotate the camera about the view up vector centered at the focal point. Note that the view up vector is not necessarily perpendicular to the direction of projection."<<std::endl<<std::endl;
    
    std::clog<<" --camelev   [double]   (optional) Rotate the camera about the cross product of the direction of projection and the view up vector centered on the focal point."<<std::endl<<std::endl;
    
    
    std::clog<<" --zoom   [value] (optional) In perspective mode, decrease the view angle by the specified factor. In parallel mode, decrease the parallel scale by the specified factor. A value greater than 1 is a zoom-in, a value less than 1 is a zoom-out"<<std::endl<<std::endl;
    
    
    std::clog<<" --glyphs   [name]  (optional) Visualize pixel,sphere,cone,cylender or cube if the numer of point are less than 1000 (--glyphs pixel)"<<std::endl<<std::endl;
    
    
    std::clog<<" --scaleglyphs     (optional)  Enables the geometrical form to be scaled with a scalar field."<<std::endl<<std::endl;
    
    std::clog<<" --radius  [value]  (optional) Select radius size"<<std::endl<<std::endl;
    
    std::clog<<" --height [value]  (optional) Select height size"<<std::endl<<std::endl;
    
    std::clog<<" --radiusscalar  [field]  (optional) It is used with scaleglyps and don't work with pixel. Is used to scale the radiusscalar of the glyphs"<<std::endl<<std::endl;
    
    std::clog<<" --heightscalar  [field]  (optional) It is used with scaleglyps and don't work with pixel and sphere. Is used to scale the heightscalar of the glyphs"<<std::endl<<std::endl;
    
    
    std::clog<<" --opacity   [value]   (optional) This option changes the opacity factor. The default value is 0.66"<<std::endl<<std::endl;
    
    std::clog<<" --scale     (optional) Select normalization"<<std::endl<<std::endl;
    
    std::clog<<" --volume  (optional) To visualize volume and have header file with cell size and computational size "<<std::endl<<std::endl;
    
    std::clog<<" --vrendering   (optional) To visualize volume rendering you must select this command with --volume and --vrenderingfield "<<std::endl<<std::endl;
    
    std::clog<<" --vrenderingfield [field]   (optional) To visualize volume rendering you must select this command with --volume and --vrendering (--volume --vrendering --vrenderingfield scalar0)"<<std::endl<<std::endl;
    
    std::clog<<"--shadow (optional) Enables shadow view in the rendering view. "<<std::endl<<std::endl;

    std::clog<<"--autorange (Optional) Enables auto range feature, which allows to get an automatically computed view range for the image. "<<std::endl<<std::endl;
    
    std::clog<<"-autorangemin [value] (Optional) Sets the minimum value of view range for the image, it is ignored if --autorange option is not enabled." << std::endl<<std::endl;
    
    std::clog<<"-autorangemax [value] (Optional) Sets the maximum value of view range for the image, it is ignored if --autorange option is not enabled." << std::endl<<std::endl;
    
    std::clog<<" --slice   (optional) If you want visualize orthoslice you must select this command with --volume --sliceplane and --vslicefield "<<std::endl<<std::endl;
    
    std::clog<<" --slicefield [field]   (optional) If you want visualize orthoslice you must select this command with --volume --sliceplane and --slice  (--volume --slice --sllice plane x --slicefield scalar0)"<<std::endl<<std::endl;
    
    std::clog<<" --sliceplane [plane]   (optional) If you want visualize orthoslice you must select this command with --volume --slicefield and --slice.You shold select x y or z  (--volume --slice --sllice plane y --slicefield scalar0)"<<std::endl<<std::endl;
    
    std::clog<<" --sliceposition [position]   (optional) If you want visualize orthoslice you must select this command with --volume --slicefield --sliceplane and --slice.If you don't select slice position defoult value is zero.This number can have value between zero and computational size value.  (--volume --slice --sliceplane z --slicefield scalar0 --sliceposition 5 )"<<std::endl<<std::endl;
    
    std::clog<<" --isosurface   (optional) If you want visualize isosurface you must select this command with --volume --isosurface and --isosurfacefield "<<std::endl<<std::endl;
    
    std::clog<<" --isosurfacefield [field]   (optional) If you want visualize isosurface you must select this command with --volume  --isosurface and --isosurfacefield  (--volume --isosurface --isosurfacefield scalar0 --isosurfacevalue 0.8 )"<<std::endl<<std::endl;
    
    std::clog<<" --isosurfacevalue [field]   (optional) If you want visualize isosurface you must select this command with --volume  --isosurface --isosurfacevalue and --isosurfacefield  "<<std::endl<<std::endl;
    
    std::cout<<"--history (optional) create an XML file which contains the history of operations performed (default create hist.xml file)"<<std::endl<<std::endl;
    
    std::cout<<"--historyfile [filename]   (optional) Change default history file name  and or directory "<<std::endl<<std::endl;
    
    /*  std::clog<<" --savepar     (optional) Is possible save all the parameter used in a txt file. (--savepar)"<<std::endl<<std::endl;
     
     
     std::clog<<" --loadpar   [path txt file]  (optional)  Is possible load a txt file where the are all the parameter that you want for the visualization.  (--loadpar /home/pippo/ciccio.bin.txt)"<<std::endl<<std::endl;*/
    
    
    /*  std::clog<<" --help "<<std::endl<<std::endl;*/
    
    
    return;
}

//---------------------------------------------------------------------
void OptionsSetter::saveParameters()
//---------------------------------------------------------------------

{
    
    std::string txtParamiters = m_vServer.path + ".txt";
    //std::clog <<"paramiters="<<txtParamiters<<std::endl;
    
    std::ofstream outParam(txtParamiters.c_str());
    if(!outParam)  return;
    if (m_vServer.volume!="none")
    {
        outParam << "--glyphs " << m_vServer.glyphs<< std::endl;
        
        if (m_vServer.nGlyphs!=0 )
            outParam << "--radius " << m_vServer.radius<< std::endl;
        
        if (m_vServer.nGlyphs!=0 && m_vServer.nGlyphs!=1)
            outParam << "--height " << m_vServer.height<< std::endl;
    }
    
    if(m_vServer.color=="yes")
    {
        outParam << "--color" <<  std::endl;
        outParam << "--colortable " << m_vServer.colorTable<< std::endl;
        
        if (m_vServer.volume!="none")
            outParam << "--colorscalar " << m_vServer.colorScalar<< std::endl;
        
        if(m_vServer.uselogscale=="yes")
            outParam << "--logscale" <<  std::endl;
    }
    
    if(m_vServer.elevation!=0 & m_vServer.azimuth!=0 & m_vServer.zoom!=1)
    {
        outParam << "--camelev " << m_vServer.elevation<< std::endl;
        outParam << "--camazim " << m_vServer.azimuth<< std::endl;
        outParam << "--zoom " << m_vServer.zoom<< std::endl;
    }
    if (m_vServer.volume!="none")
    {
        outParam << "--opacity " << m_vServer.opacity<< std::endl;
        
        outParam << "-x " << m_vServer.xField<< std::endl;
        outParam << "-y " << m_vServer.yField<< std::endl;
        outParam << "-z " << m_vServer.zField<< std::endl;
    }
    outParam << "--out " << m_vServer.imageName<< std::endl;
    
    if(m_vServer.vector=="yes" && m_vServer.volume!="none")
    {
        outParam << "--vector" << std::endl;
        outParam << "--vx " << m_vServer.xVectorField<< std::endl;
        outParam << "--vy " << m_vServer.yVectorField<< std::endl;
        outParam << "--vz " << m_vServer.zVectorField<< std::endl;
    }
    
    if (m_vServer.scale=="yes")
        outParam << "--scale" << std::endl;
    
    if (m_vServer.scaleGlyphs=="yes"&& m_vServer.volume=="yes")
        
    {
        outParam << "--scaleglyphs" << std::endl;
        if (m_vServer.nGlyphs!=0 )
            outParam << "--radiusscalar " << m_vServer.radiusscalar<< std::endl;
        
        if (m_vServer.nGlyphs!=0 && m_vServer.nGlyphs!=1)
            outParam << "--heightscalar " << m_vServer.heightscalar<< std::endl;
    }
    if (m_vServer.noDefault=="yes" && m_vServer.slice!="yes")
        outParam << "--nodefault" << std::endl;
    
    if (m_vServer.volume=="yes")
    {
        outParam << "--volume" << std::endl;
        if(m_vServer.vRendering=="yes" && m_vServer.vRenderingField !="none")
        {
            outParam << "--rendering "<< std::endl;
            outParam << "--renderingfield "<<m_vServer.vRenderingField << std::endl;
        }
        else if(m_vServer.slice=="yes" && m_vServer.sliceField !="none")
        {
            outParam << "--slice "<< std::endl;
            outParam << "--slicefield "<<m_vServer.sliceField << std::endl;
            outParam << "--sliceplane "<<m_vServer.slicePlane << std::endl;
            outParam << "--sliceposition "<<m_vServer.slicePosition << std::endl;
        }
        
        else if(m_vServer.isosurface=="yes" && m_vServer.isosurfaceField !="none")
        {
            outParam << "--isosurface "<< std::endl;
            outParam << "--isosurfacefield "<<m_vServer.isosurfaceField << std::endl;
            outParam << "--isosurfacevalue "<<m_vServer.isosurfaceValue << std::endl;
            
        }
    }
    
    outParam <<m_vServer.path<< std::endl;
    
    
    outParam.close();
    
    return;
    
}

//---------------------------------------------------------------------
void OptionsSetter::loadParameters(std::string path)
//---------------------------------------------------------------------
{
    
    
    std::fstream inFile;
    inFile.open(path.c_str());
    
    if(!inFile)
    {
        std::clog<<"the txt path is incorrect, please try again  or --help for help "<<std::endl;
        return ;
    }
    
    std::string ext=getExt(path.c_str());
    if (ext!="txt")
    {
        std::clog<<"please give me a txt file "<<std::endl;
        return ;
    }
    
    std::string::size_type index = std::string::npos;
    std::string fields;
    std::string tmp = "";
    std::vector<std::string> lineData;
    std::vector<std::string> optToLoad;
    int i ,j;
    int indexp = tmp.find('#');
    while(!inFile.eof())
    {
        tmp = "";
        
        getline(inFile, tmp);
        
        findAndReplace(tmp, '\t', ' ');
        tmp = trim(tmp);
        
        indexp = tmp.find('#');
        
        if(indexp != std::string::npos)
            tmp.erase(index);
        
        if(tmp.compare(""))
            lineData.push_back(tmp);
        // std::clog<<tmp<<std::endl;
    }
    
    inFile.close();
    
    for(i = 0; i <lineData.size()-1; i++)
    {
        std::stringstream ss;
        ss << lineData[i];
        
        
        for(j = 0; j < 2; j++)
        {
            tmp="";
            ss >> tmp;
            optToLoad.push_back(tmp);
            // // std::clog<<tmp<<std::endl;
            
        }
        
    }
    std::stringstream ss;
    ss << lineData[lineData.size()-1];
    tmp="";
    ss >> tmp;
    optToLoad.push_back(tmp);
    //std::clog<<tmp<<std::endl;
    
    
    parseOption ( optToLoad);
    
    
}

//---------------------------------------------------------------------
bool OptionsSetter::readHeader()
//---------------------------------------------------------------------
{
    int i = 0;
    
    std::string headerFileName = m_vServer.path + ".head";
    std::ifstream inHeader;
    inHeader.open(headerFileName.c_str());
    
    if(!inHeader)
        return false;
    
    const unsigned long int BUFFER_SIZE = 73;
    char buffer[BUFFER_SIZE];
    
    const int BUFFER_FULL = BUFFER_SIZE - 1;
    std::streamsize sSize = BUFFER_FULL;
    
    
    inHeader >>m_vServer.dataType;
    
    if(m_vServer.dataType=="double" ||m_vServer.dataType=="d")
    {
        std::cerr<<"VisIVOViewer accepts only float. Use VisIVOImporter to convert in float data format."<<std::endl;
        return false;
    }
    
    inHeader >> m_vServer.nCols;
    
    std::string tmp = "";
    getline(inHeader, tmp); //!to remove the carrige return character  (\r) from the last line
    getline(inHeader, tmp);
    std::stringstream sstmp(tmp);
    sstmp >> m_vServer.nRows;
    
    if(!(sstmp.eof()))
    {
        //    m_isVolume=true;
        
        sstmp >>  m_vServer.comp[0];
        sstmp >>  m_vServer.comp[1];
        sstmp >>  m_vServer.comp[2];
        
        sstmp >> m_vServer.size[0];
        sstmp >>  m_vServer.size[1];
        sstmp >>  m_vServer.size[2];
        
        tryToSetDimension();
    }
    if(m_vServer.volume=="yes")
        if(m_vServer.comp[0]*m_vServer.comp[1]*m_vServer.comp[2] != m_vServer.nRows)
        {
            std::cerr<<"Invalid header file for volume: nRows and cells don't match"<<std::endl;
            std::cerr<<"Operation aborted"<<std::endl;
            return false;
        }
    
    inHeader >> m_vServer.endian;
    if(m_vServer.endian=="b") m_vServer.endian="big";
    if(m_vServer.endian=="l") m_vServer.endian="little";
    
#ifdef VSBIGENDIAN
    m_vServer.systemEndianism="big";
#endif
    
    if(m_vServer.endian=="big" && m_vServer.systemEndianism=="little")
        m_vServer.needSwap=true;
    if(m_vServer.endian=="little" && m_vServer.systemEndianism=="big")
        m_vServer.needSwap=true;
    
    
    //std::clog <<"rows="<<m_vServer.nRows<<std::endl;
    
    //to throw away the cr after the number if rows
    inHeader.getline(buffer, BUFFER_SIZE);
    
    //reads the field names
    m_vServer.fieldNames.clear();
    
    std::string name = "";
    
    for(i = 0; i < m_vServer.nCols && !(inHeader.eof()); i++)
    {
        name = "";
        
        sSize = BUFFER_FULL;
        
        while(sSize == BUFFER_FULL)
        {
            inHeader.clear();
            inHeader.getline(buffer, BUFFER_SIZE);
            sSize = inHeader.gcount();
            
            if(sSize < 1)
                break;
            
            
            name.append(buffer);
            size_t position=name.find( '\r');
            //std::clog <<name<<std::endl;
            if(position!=-1)
            {
                name.replace(position,2,"");
                //std::clog << name<<std::endl;
            }
            
            if(sSize == BUFFER_FULL && buffer[BUFFER_SIZE - 2] == 0)
                sSize = 0;
        }
        m_vServer.fieldNames.push_back(name);
    }
    if(m_vServer.fieldNames.size()<3 && m_vServer.volume !="yes")
    {
        std::cerr<<"Warning: No 3D table is given. Images cannot be created"<<std::endl;
        return false;
    }
    
    if(i != m_vServer.nCols)
    {
        m_vServer.nCols = 0;
        m_vServer.nRows = 0;
        m_vServer.fieldNames.clear();
    }
    
    inHeader.close();
    
    return true;
}

//---------------------------------------------------------------------
void OptionsSetter::setGlyphs()
//---------------------------------------------------------------------
{
    
    if(m_vServer.glyphs=="pixel")
        m_vServer.nGlyphs=0;
    
    else if(m_vServer.glyphs=="sphere")
        m_vServer.nGlyphs=1;
    
    else if(m_vServer.glyphs=="cone")
        m_vServer.nGlyphs=2;
    
    else if(m_vServer.glyphs=="cylinder")
        m_vServer.nGlyphs=3;
    
    else if(m_vServer.glyphs=="cube")
        m_vServer.nGlyphs=4;
    
    return;
}

//---------------------------------------------------------------------
void OptionsSetter::setColorLut()
//---------------------------------------------------------------------
{
    
    if(m_vServer.colorTable=="default")
        m_vServer.nColorTable=0;
    
    else if(m_vServer.colorTable=="default_step")
        m_vServer.nColorTable=1;
    
    else if(m_vServer.colorTable=="efield")
        m_vServer.nColorTable=2;
    
    else if(m_vServer.colorTable=="glow")
        m_vServer.nColorTable=3;
    
    else if(m_vServer.colorTable=="gray")
        m_vServer.nColorTable=4;
    
    else if(m_vServer.colorTable=="min_max")
        m_vServer.nColorTable=5;
    
    else if(m_vServer.colorTable=="physics_contour")
        m_vServer.nColorTable=6;
    
    else if(m_vServer.colorTable=="pure_red")
        m_vServer.nColorTable=7;
    
    else if(m_vServer.colorTable=="pure_green")
        m_vServer.nColorTable=8;
    
    else if(m_vServer.colorTable=="pure_blue")
        m_vServer.nColorTable=9;
    
    else if(m_vServer.colorTable=="run1")
        m_vServer.nColorTable=10;
    
    else if(m_vServer.colorTable=="run2")
        m_vServer.nColorTable=11;
    
    else if(m_vServer.colorTable=="sar")
        m_vServer.nColorTable=12;
    
    else if(m_vServer.colorTable=="temperature")
        m_vServer.nColorTable=13;
    
    else if(m_vServer.colorTable=="tensteps")
        m_vServer.nColorTable=14;
    
    else if(m_vServer.colorTable=="volren_glow")
        m_vServer.nColorTable=15;
    
    else if(m_vServer.colorTable=="volren_green")
        m_vServer.nColorTable=16;
    
    else if(m_vServer.colorTable=="volren_rgb")
        m_vServer.nColorTable=17;
    
    else if(m_vServer.colorTable=="volren_twolevel")
        m_vServer.nColorTable=18;
    else if(m_vServer.colorTable=="all_yellow")
        m_vServer.nColorTable=19;
    else if(m_vServer.colorTable=="all_cyane")
        m_vServer.nColorTable=20;
    else if(m_vServer.colorTable=="all_violet")
        m_vServer.nColorTable=21;
    else if(m_vServer.colorTable=="all_white")
        m_vServer.nColorTable=22;
    else if(m_vServer.colorTable=="all_black")
        m_vServer.nColorTable=23;
    else if(m_vServer.colorTable=="all_red")
        m_vServer.nColorTable=24;
    else if(m_vServer.colorTable=="all_green")
        m_vServer.nColorTable=25;
    else if(m_vServer.colorTable=="all_blu")
        m_vServer.nColorTable=26;
    
    else if(m_vServer.colorTable!="none")
    {
        m_vServer.nColorTable=-1; //external file table
        std::ifstream extPal(m_vServer.colorTable.c_str());
        if(!extPal)
        {
            std::cerr<<"Warning: Invalid filename for external palette "<<m_vServer.colorTable<<std::endl;
            std::cerr<<"Defaul palette will be used"<<std::endl;
            m_vServer.nColorTable=0;
            return;
        }
        
        int idOp=0;
        int  palId=-1;
        bool hsv=false;
        while(!extPal.eof())
        {
            std::string tmp = "";
            getline(extPal, tmp);
            std::stringstream sstmp(tmp);
            if(sstmp.str()=="") continue;
            std::string word;
            std::vector<std::string> countLines;
            double  palR=0., palG=0.,palB=0.,palA=1.0;
            while(sstmp >> word)
                countLines.push_back(word);
            if(palId==-1 && trim(word)=="RGB") continue;
            if(palId==-1 && trim(word)=="HSV")
            {
                hsv=true;
                continue;
            }
            if(trim(countLines[0]).substr(0,1) =="#") continue;
            if(idOp==0 && countLines.size()<=5 && countLines.size()>=3)
                idOp=countLines.size(); //Id+R+G+B+A, R+G+B+A; R+G+B (or ID+HSV+A) range 0.0 - 1.0
            if(idOp!=countLines.size())
            {
                std::cerr<<"Warning: Invalid  external palette "<<m_vServer.colorTable<<std::endl;
                std::cerr<<"Palette file must contain the same number of entry (5,4 or 3) for each line"<<std::endl;
                std::cerr<<"5 Values: IdRGBA; 4 Values RGBA; 3 Values RGB only"<<std::endl;
                std::cerr<<"5 Values: IdHSVA; 4 Values HSVA; 3 Values HSV only"<<std::endl;
                std::cerr<<"Default palette will be used"<<std::endl;
                m_vServer.nColorTable=0;
                return;
                
            }
            int indexWord=0;
            if(idOp<5) palId++;
            switch(idOp)
            {
                case 5:
                {
                    std::stringstream sstmpId(countLines[indexWord]);
                    int palIdOld=palId;
                    sstmpId >> palId;
                    if(palId<=palIdOld)
                    {
                        std::cerr<<"Warning invalid palette Id sequence is given: "<< palId;
                        std::cerr<<palId <<" is lower or equal to previous value "<<palIdOld<<std::endl;
                        std::cerr<<"Default palette will be used"<<std::endl;
                        m_vServer.nColorTable=0;
                        return;
                    }
                    if(palId<0)
                    {
                        std::cerr<<"Warning invalid palette Id is given: "<< palId;
                        std::cerr<<"Default palette will be used"<<std::endl;
                        m_vServer.nColorTable=0;
                        return;
                    }
                    indexWord++;
                }
                case 4:
                {
                    std::stringstream sstmpA(countLines[indexWord+3]);
                    sstmpA >> palA;
                    if(palA<0.) palA=0.;
                    if(palA>1.) palA=1.;
                }
                case 3:
                {
                    std::stringstream sstmpR(countLines[indexWord]);
                    sstmpR >> palR;
                    if(palR<0.) palR=0.;
                    if(palR>1.) palR=1.;
                    std::stringstream sstmpG(countLines[indexWord+1]);
                    sstmpG >> palG;
                    if(palG<0.) palG=0.;
                    if(palG>1.) palG=1.;
                    std::stringstream sstmpB(countLines[indexWord+2]);
                    sstmpB >> palB;
                    if(palB<0.) palB=0.;
                    if(palB>1.) palB=1.;
                    if(idOp==3)
                    {
                        palA=1.;
                        break;
                    }
                }
            } //end swwitch
            if(hsv)
            {
                double h=palR;
                double s=palG;
                double v=palB;
                HSVtoRGB(&palR,&palG,&palB,h,s,v);
            }
            m_vServer.extPalId.push_back(palId);
            m_vServer.extPalR.push_back(palR);
            m_vServer.extPalG.push_back(palG);
            m_vServer.extPalB.push_back(palB);
            m_vServer.extPalA.push_back(palA);
            
        } //while
        if(m_vServer.extPalId.size()<=1)
        {
            std::cerr<<"Warning: Invalid filename for external palette "<<m_vServer.colorTable<<std::endl;
            std::cerr<<"Defaul palette will be used"<<std::endl;
            m_vServer.nColorTable=0;
        }
    } //else if(m_vServer.colorTable!="none")
    return;
}
//---------------------------------------------------------------------
void OptionsSetter::setNumberImages()
//---------------------------------------------------------------------
{
    if (m_vServer.noDefault!="yes")
    {
        int numImage=5;
        
        if(m_vServer.azimuth==0 && m_vServer.elevation==0 && m_vServer.zoom==1)
            numImage=4;
        
        m_vServer.numImageToLoad=numImage-1;
        
    }
    
    else
    {
        
        if (m_vServer.azimuth!=0 || m_vServer.elevation!=0 || m_vServer.zoom!=1)
            m_vServer.numImageToLoad=0;
        
        else
            m_vServer.numImageToLoad=0;
        
    }
    if(m_vServer.cycle)
        if(m_vServer.cycleSkipFrom>=0 && m_vServer.cycleSkipTo >m_vServer.cycleSkipFrom)
            m_vServer.numImageToLoad=m_vServer.cycleSkipTo-m_vServer.cycleSkipFrom;
        else{
            m_vServer.numImageToLoad=m_cycleAzimuth.size();
            if(m_vServer.slice=="yes")
            {
                if(m_vServer.sliceOrthoNormal)
                    m_vServer.numImageToLoad=m_cyclePlane.size();
                else
                    m_vServer.numImageToLoad=m_cyclePlanePointNormal.size();
                
            }
        }
    return;
}

//---------------------------------------------------------------------
int OptionsSetter::images()
//---------------------------------------------------------------------
{
    // check for output directory
    std::ofstream testout;
    testout.open(m_vServer.imageName.c_str());
    if(!testout.good())
    {
        std::cerr<<"Invalid output filename or directory. Viewer aborted"<<std::endl;
        return -1;
    } else
        remove(m_vServer.imageName.c_str());
    //
    
    Pipe *pp;
    int errPipe=0;
    int countCycle=0;
    std::string VSCycleFileList;
    std::string cycleRootName,lfnCycleRootName;
    cycleRootName=m_vServer.imageName;

    
    size_t ext;
    ext=cycleRootName.find(".png");
    if(ext!=std::string::npos) cycleRootName=cycleRootName.substr(0,ext);
    ext=cycleRootName.find(".tga");
    if(ext!=std::string::npos) cycleRootName=cycleRootName.substr(0,ext);
    if(m_vServer.lfnout!="none") lfnCycleRootName=m_vServer.lfnout;
    
    std::ofstream ouFileList;
    bool onlyOneImage=false;
    //  std::clog<<"Start IMAGES!"<<std::endl;
    if(m_vServer.cycle) 
    {	std::stringstream tmp;
        tmp<<m_vServer.cycleOffset;
        VSCycleFileList=getDir(cycleRootName)+"VSCycleImage"+tmp.str()+".txt";
        ouFileList.open(VSCycleFileList.c_str());
    }
    
    if (m_vServer.splotch)
    {
#ifndef SPLVISIVO
        std::cerr<<"VisIVO must be compiled with SPLVISIVO flag to enable Splotch"<<std::endl;
        return 1;
#else
        if(!m_vServer.cycle)
            m_vServer.numImageToLoad=1;
        if(m_vServer.cycle && (m_vServer.cycleSkipFrom>=0 && m_vServer.cycleSkipTo >m_vServer.cycleSkipFrom))
			countCycle=m_vServer.cycleSkipFrom;
        for(int i=0;i<m_vServer.numImageToLoad;i++)
        {
            if(m_vServer.cycle)
            {
                if(countCycle==m_cycleAzimuth.size())
                    break;
                std::stringstream suffix;
                int isuffix;
                m_vServer.azimuth=m_cycleAzimuth[countCycle];
                m_vServer.elevation=m_cycleElevation[countCycle];
                m_vServer.zoom=m_cycleZoom[countCycle];
                if(m_vServer.zoom<0.0)
                    std::cerr<<"Warning: negative zoom value is given, could be ignored."<<std::endl;
                if(m_vServer.setCameraPos)
                {
                    m_vServer.cameraPos[0]=m_cyclecamPos0[countCycle];
                    m_vServer.cameraPos[1]=m_cyclecamPos1[countCycle];
                    m_vServer.cameraPos[2]=m_cyclecamPos2[countCycle];
                }
                if(m_vServer.setCameraFocalPoint)
                {
                    m_vServer.cameraFocalPoint[0]=m_cyclecamFp0[countCycle];
                    m_vServer.cameraFocalPoint[1]=m_cyclecamFp1[countCycle];
                    m_vServer.cameraFocalPoint[2]=m_cyclecamFp2[countCycle];
                }
                if(m_vServer.setCameraRoll) 
                    m_vServer.cameraRoll=m_cyclecamRoll[countCycle];
                
                isuffix=m_vServer.cycleOffset+i;
                suffix<<isuffix;
                std::string zeropad;
                if(isuffix<10) zeropad="000000";
                if(isuffix>=10 && isuffix<100 ) zeropad="00000";
                if(isuffix>=100 && isuffix<1000 ) zeropad="0000";
                if(isuffix>=1000 && isuffix<10000 ) zeropad="000";
                if(isuffix>=10000 && isuffix<100000 ) zeropad="00";
                if(isuffix>=100000 && isuffix<1000000 ) zeropad="0";
                if(isuffix>=1000000) zeropad="";
                m_vServer.imageName=cycleRootName+zeropad+suffix.str()+".tga";
#ifdef GLITE
                if(m_vServer.lfnout!="none")
                    m_vServer.lfnout=lfnCycleRootName+zeropad+suffix.str()+".tga";
#endif
                m_vServer.numImage=4;  
                ouFileList<<m_vServer.imageName<<std::endl;
                countCycle++;
            }else
            {
                if ((m_vServer.azimuth!=0 || m_vServer.elevation!=0 || m_vServer.zoom!=1))
                    m_vServer.numImage=4; 
                if(m_vServer.imageName.find(".tga") == std::string::npos)
                    m_vServer.imageName=m_vServer.imageName+".tga";
                if(m_vServer.lfnout!="none" && m_vServer.lfnout.find(".tga") == std::string::npos)
                    m_vServer.lfnout=m_vServer.lfnout+".tga";
            }
            
            Pipe *ppcamera;
            ppcamera=new SplotchPipeCamera(m_vServer); 
            ppcamera->getCamera(&m_spCamera); 
            
            for (int i=0;i<3;i++)
            {
                m_vServer.spPosition[i]=m_spCamera.position[i];
                m_vServer.spLookat[i]=m_spCamera.lookat[i];
                m_vServer.spSky[i]=m_spCamera.sky[i];
            }
            
            if(!m_vServer.fovIsGiven)
                m_vServer.fov=m_spCamera.fov;
            
            m_vServer.spClip[0]=m_spCamera.clip[0];
            m_vServer.spClip[1]=m_spCamera.clip[1];
            m_vServer.spClip[0]=0;
            m_vServer.spClip[1]=1.0e+13;
            
            delete ppcamera;
            pp=new SplotchPipe(m_vServer);  
            errPipe=pp->createPipe();
#ifdef GLITE	  
            if(m_vServer.lfnout!="none")
            { 
                if(m_vServer.VO=="none")
                {
                    std::cerr<<"Invalid VO option"<<std::endl;
                    return 1;
                }
                //	    std::clog<<m_vServer.imageName<<" "<<m_vServer.lfnout<<std::endl;
                if(!remoteLfn(m_vServer.imageName, m_vServer.se, m_vServer.VO, m_vServer.lfnout, false))
                {
                    std::cerr<<"Invalid copy on gLite Catalogue"<<std::endl;
                    return 1;
                }
            }
#endif
            delete pp;
            if(errPipe<0) return invalidImage;
        }   
        if(m_vServer.cycle) 
            ouFileList.close();
        
        return 0;
#endif
    }
    
    if(m_vServer.cycle && (m_vServer.cycleSkipFrom>=0 && m_vServer.cycleSkipTo >m_vServer.cycleSkipFrom))
        countCycle=m_vServer.cycleSkipFrom;
    for (int i=0;i<=m_vServer.numImageToLoad;i++)
    {
        if(m_vServer.cycle)
        {
            if(i==m_vServer.numImageToLoad)
                break;
            if(m_vServer.slice!="yes" && countCycle==m_cycleAzimuth.size())
                break;
            if(m_vServer.slice=="yes" && m_vServer.sliceOrthoNormal && countCycle==m_cyclePlane.size())
                break;
            if(m_vServer.slice=="yes" && !m_vServer.sliceOrthoNormal && countCycle==m_cyclePlanePointNormal.size())
                break;
            std::stringstream suffix;
            int isuffix;
            if(m_vServer.slice!="yes")
            {
                m_vServer.azimuth=m_cycleAzimuth[countCycle];
                m_vServer.elevation=m_cycleElevation[countCycle];
                m_vServer.zoom=m_cycleZoom[countCycle];
                if(m_vServer.zoom<0.0)
                    std::cerr<<"Warning: negative zoom value is given, could be ignored."<<std::endl;
                
                if(m_vServer.setCameraPos)
                {
                    m_vServer.cameraPos[0]=m_cyclecamPos0[countCycle];
                    m_vServer.cameraPos[1]=m_cyclecamPos1[countCycle];
                    m_vServer.cameraPos[2]=m_cyclecamPos2[countCycle];
                }
                if(m_vServer.setCameraFocalPoint)
                {
                    m_vServer.cameraFocalPoint[0]=m_cyclecamFp0[countCycle];
                    m_vServer.cameraFocalPoint[1]=m_cyclecamFp1[countCycle];
                    m_vServer.cameraFocalPoint[2]=m_cyclecamFp2[countCycle];
                }
                if(m_vServer.setCameraRoll) 
                    m_vServer.cameraRoll=m_cyclecamRoll[countCycle];
                
                
            } else { 
                if(m_vServer.sliceOrthoNormal) 
                    m_vServer.slicePosition=m_cyclePlane[countCycle];
                else{
                    std::stringstream ssppn(m_cyclePlanePointNormal[countCycle]);
                    for(int ipp=0;ipp<3;ipp++)
                        ssppn>>m_vServer.slicePlanePointValue[ipp];
                    for(int ipp=0;ipp<3;ipp++)
                        ssppn>>m_vServer.slicePlaneNormalValue[ipp];
                }
            }   
            isuffix=m_vServer.cycleOffset+i;
            suffix<<isuffix;
            std::string zeropad;
            if(isuffix<10) zeropad="000000";
            if(isuffix>=10 && isuffix<100 ) zeropad="00000";
            if(isuffix>=100 && isuffix<1000 ) zeropad="0000";
            if(isuffix>=1000 && isuffix<10000 ) zeropad="000";
            if(isuffix>=10000 && isuffix<100000 ) zeropad="00";
            if(isuffix>=100000 && isuffix<1000000 ) zeropad="0";
            if(isuffix>=1000000) zeropad="";
            m_vServer.imageName=cycleRootName+zeropad+suffix.str();
            m_vServer.numImage=4;
#ifdef GLITE
            if(m_vServer.lfnout!="none")
                m_vServer.lfnout=lfnCycleRootName+zeropad+suffix.str()+".png";
#endif
            ouFileList<<m_vServer.imageName<<".png"<<std::endl;
            countCycle++;
        }	
        else if ((m_vServer.azimuth!=0 || m_vServer.elevation!=0 || m_vServer.zoom!=1) && m_vServer.noDefault=="yes")
            m_vServer.numImage=4;
        else
        {
            m_vServer.numImage=i;
            if(m_vServer.lfnout!="none" && m_vServer.lfnout.find(".png") == std::string::npos)
                m_vServer.lfnout=m_vServer.lfnout+".png";
        }
        
        int cycleStereo=1;
        std::string rootName=m_vServer.imageName;
        //      std::clog<<rootName<<std::endl;
        
        if(m_vServer.stereo && m_vServer.stereoMode=="CrystalEyes") cycleStereo=2;
        
        for(int k=0;k<cycleStereo;k++)
        {
            if(m_vServer.stereoMode=="CrystalEyes")
            {
                m_vServer.stereoImg=k;
                if(k==0) m_vServer.imageName=rootName+"_r";
                if(k==1) m_vServer.imageName=rootName+"_l";
            }  
            
            if(m_vServer.vtkImage=="yes")
            {
                pp=new VtkImagePipe(m_vServer);
                errPipe=pp->createPipe();
                pp->destroyAll();
                delete pp; 
                if(errPipe<0) return invalidImage;
                break;
            }
            if(m_vServer.volume=="yes" && m_vServer.vRendering=="yes" && m_vServer.vRenderingField !="none")
                pp=new VolumePipe(m_vServer);
            
            else if(m_vServer.volume=="yes" && m_vServer.isosurface=="yes" && m_vServer.isosurfaceField !="none")
                pp=new IsosurfacePipe(m_vServer);
            
            else  if (m_vServer.volume=="yes" && m_vServer.slice=="yes" && m_vServer.sliceField !="none" && (m_vServer.slicePlane !="none" || m_vServer.slicePlanePoint !="none" || m_vServer.slicePlaneNormal !="none"))
            {
                m_vServer.numImage=4; // no default images
                pp=new SlicerPipe(m_vServer);
                onlyOneImage=true;
            }
            
            else  if (m_vServer.vector=="yes")
                pp=new VectorPipe(m_vServer);
            else  if (m_vServer.mode=="smooth")
                pp=new PointsSmoothPipe(m_vServer);
            
            else{
                pp=new PointsPipe(m_vServer);
            }
            if(useInMemory){
                pp->setUseMemory(true);
                pp->setMemTable(table);
            } 
            errPipe=pp->createPipe();
            std::string fName = pp->saveImageAsPng(m_vServer.numImage);
            
            //Filename for history!
            outFilename.insert(outFilename.end(), fName);
            //std::cout<<"************"<<fName<<"*****************"<<std::endl;
            
#ifdef GLITE
            if(m_vServer.lfnout!="none")
            { 
                if(m_vServer.VO=="none")
                {
                    std::cerr<<"Invalid VO option"<<std::endl;
                    return 1;
                }
                std::string tmpLfn=m_vServer.lfnout;
                std::string tmpImgName=m_vServer.imageName;
                std::stringstream ssnumImage;
                std::string numImage="";
                if(m_vServer.noDefault!="yes")
                {
                    ssnumImage<<m_vServer.numImage;
                    ssnumImage>>numImage;
                }
                size_t found=tmpLfn.find(".png");
                
                if(found== std::string::npos)
                    tmpLfn += numImage+".png";
                else
                    tmpLfn=tmpLfn.substr(0,found)+numImage+".png"; 
                
                found=tmpImgName.find(".png");
                if(found== std::string::npos)
                    tmpImgName += numImage+".png";
                else
                    tmpImgName=tmpImgName.substr(0,found)+numImage+".png";
                //	      std::clog<<tmpImgName<<" "<<tmpLfn<<std::endl;
                if(!remoteLfn(tmpImgName, m_vServer.se, m_vServer.VO, tmpLfn, false))
                {
                    std::cerr<<"Invalid copy on gLite Catalogue"<<std::endl;
                    return 1;
                }
            }
#endif
            
            pp->destroyAll();
            delete pp; 
            if(errPipe<0) return invalidImage;
            
            if(onlyOneImage) break;
        } // for( k=0;
    }
    
    if(m_vServer.cycle) 
        ouFileList.close();
#ifdef GLITE
    if(m_inputLfnGiven)
    {
        remove(m_vServer.path.c_str());
        m_vServer.path +=".head";
        remove(m_vServer.path.c_str());
        m_inputLfnGiven=false;
    }
#endif
    
    return 0;
}


//----------------------------
void OptionsSetter::tryToSetDimension()

//---------------------------
{
    if (m_vServer.comp[0]*m_vServer.comp[1]*m_vServer.comp[2]==m_vServer.nRows)
        return;
    
    double size=pow((double)m_vServer.nRows,1.0/3);
    
    double dimensionUp=ceil(size);
    //std::clog<<"dimensionUp="<<dimensionUp<<std::endl;
    
    double dimensionDown= floor(size);
    //std::clog<<"dimensionDown="<<dimensionDown<<std::endl;
    
    if(m_vServer.nRows==pow(dimensionUp,3))
        m_vServer.comp[0]=m_vServer.comp[1]=m_vServer.comp[2]=dimensionUp;
    
    if (m_vServer.nRows==pow(dimensionDown,3))
        m_vServer.comp[0]=m_vServer.comp[1]=m_vServer.comp[2]=dimensionDown;
    
    else 
        m_vServer.comp[0]=m_vServer.comp[1]=m_vServer.comp[2]=0.0 ;
    
}
//----------------------------
void OptionsSetter::readCycleFile()

//---------------------------
{
	std::ifstream inFile;
	double camPathValue[7];
	inFile.open(m_vServer.cycleFile.c_str());
	if(!inFile)
	{
		std::cerr<<"Cannot open cycle file."<<std::endl;
		return ;
	}
	if(m_vServer.slice=="yes")
	{
        if(m_vServer.sliceOrthoNormal)
        {  
            int previousValue=INVALID_READ;
            while(!inFile.eof())
            {
                int plane=INVALID_READ;
                
                std::string tmp = "";
                //  getline(inFile, tmp); //to remove the carrige return character  (\r) from the last line
                getline(inFile, tmp);
                std::stringstream sstmp(tmp);
                sstmp >> plane;
                
                if(plane >=0 && plane!=INVALID_READ && plane!=previousValue )
                {
                    m_cyclePlane.push_back(plane);
                    previousValue=plane;
                }
            }
            if(m_cyclePlane.size() <=0)
            {		
                std::cerr<<"Invalid values in cycle file."<<std::endl;
                return ;
            }
        } else{  // if Orthonormal
            while(!inFile.eof())
            {
                std::string tmp;
                //  getline(inFile, tmp); //to remove the carrige return character  (\r) from the last line
                getline(inFile, tmp);
                std::stringstream sstmp(tmp);
                int tmpcount=0;
                while(!sstmp.eof())
                { 
                    std::string tmp1;
                    sstmp>>tmp1;
                    tmpcount++;
                }  
                if(!tmp.empty() && tmpcount>=6 ) // ignore lines with lower than 6 values
                    m_cyclePlanePointNormal.push_back(tmp);
            }
            if(m_cyclePlanePointNormal.size() <=0)
            {		
                std::cerr<<"Invalid values in cycle file."<<std::endl;
                return ;
            }
            
        }
        
	}else{
        
        int fileType=-1;
        while(!inFile.eof())
        {
            float azimuth=INVALID_READ,elevation=INVALID_READ,zoom=INVALID_READ;
            std::string tmp = "";
            std::string dummy;
            //  getline(inFile, tmp); //to remove the carrige return character  (\r) from the last line
            getline(inFile, tmp);
            std::stringstream sstmp(tmp);
	        std::stringstream sstmpckft(tmp);
	        std::stringstream sstmpck(tmp);
            std::string tmpck;
            if(fileType==-1)
            {
                int wc=0;
                while(!sstmpckft.eof())
                {
                    sstmpckft>>dummy;
                    wc++;
                }  
                if(wc==3) fileType=wc;
                if(wc==4) fileType=wc;
                if(wc==7) fileType=wc;
                if(wc==8) fileType=wc;
                if(wc==10) fileType=wc;
                if(fileType==-1)
                {
                    std::cerr<<"Invalid file type for cycle, ignored"<<std::endl;
                    m_vServer.cycle=false;
                    return;
                }
            }
            if(fileType>3) m_vServer.setCameraRoll=true;
            switch(fileType)
            {
                case 3:
                {
                    std::string tmpck;
                    double camRoll=INVALID_CAM;
                    sstmpck >> tmpck;
                    if(tmpck.compare("NULL")!=0) sstmp >> azimuth;
                    sstmpck >> tmpck;
                    if(tmpck.compare("NULL")!=0) sstmp >> elevation;
                    sstmpck >> tmpck;
                    if(tmpck.compare("NULL")!=0) sstmp >> zoom;
                    if(azimuth==INVALID_READ || elevation==INVALID_READ || zoom==INVALID_READ)
                    {
                        std::cerr<<"Invalid line for cycle file , ignored"<<std::endl;
                        break;
                    }
                    m_cycleAzimuth.push_back(azimuth);
                    m_cycleElevation.push_back(elevation);
                    m_cycleZoom.push_back(zoom);
                    break;
                }
                case 4:
                {
                    double camRoll;
                    sstmpck >> tmpck;
                    if(tmpck.compare("NULL")!=0) sstmp >> azimuth;
                    sstmpck >> tmpck;
                    if(tmpck.compare("NULL")!=0) sstmp >> elevation;
                    sstmpck >> tmpck;
                    if(tmpck.compare("NULL")!=0) sstmp >> zoom;
                    sstmpck >> tmpck;
                    if(tmpck.compare("NULL")!=0) sstmp >> camRoll;
                    if(azimuth==INVALID_READ || elevation==INVALID_READ || zoom==INVALID_READ)
                    {
                        std::cerr<<"Invalid line for cycle file , ignored"<<std::endl;
                        break;
                    }
                    
                    m_cycleAzimuth.push_back(azimuth);
                    m_cycleElevation.push_back(elevation);
                    m_cycleZoom.push_back(zoom);
                    m_cyclecamRoll.push_back(camRoll);
                    break;
                }
                case 7:
                {
                    m_vServer.setCameraFocalPoint=true;
                    double camRoll=INVALID_CAM;
                    double fp[3];
                    fp[0]=INVALID_CAM;
                    fp[1]=INVALID_CAM;
                    fp[2]=INVALID_CAM;
                    sstmpck >> tmpck;
                    if(tmpck.compare("NULL")!=0) sstmp >> azimuth;
                    sstmpck >> tmpck;
                    if(tmpck.compare("NULL")!=0) sstmp >> elevation;
                    sstmpck >> tmpck;
                    if(tmpck.compare("NULL")!=0) sstmp >> zoom;
                    sstmpck >> tmpck;
                    if(tmpck.compare("NULL")!=0) sstmp >> fp[0];
                    sstmpck >> tmpck;
                    if(tmpck.compare("NULL")!=0) sstmp >> fp[1];
                    sstmpck >> tmpck;
                    if(tmpck.compare("NULL")!=0) sstmp >> fp[2];
                    sstmpck >> tmpck;
                    if(tmpck.compare("NULL")!=0) sstmp >> camRoll;
                    if(azimuth==INVALID_READ || elevation==INVALID_READ || zoom==INVALID_READ)
                    {
                        std::cerr<<"Invalid line for cycle file , ignored"<<std::endl;
                        break;
                    }
                    
                    m_cycleAzimuth.push_back(azimuth);
                    m_cycleElevation.push_back(elevation);
                    m_cycleZoom.push_back(zoom);
                    m_cyclecamFp0.push_back(fp[0]);
                    m_cyclecamFp1.push_back(fp[1]);
                    m_cyclecamFp2.push_back(fp[2]);
                    m_cyclecamRoll.push_back(camRoll);
                    break;
                }
                case 8:
                {
                    m_vServer.setCameraFocalPoint=true;
                    m_vServer.setCameraPos=true;
                    double camRoll=INVALID_CAM;
                    double cp[3];
                    cp[0]=INVALID_CAM;
                    cp[1]=INVALID_CAM;
                    cp[2]=INVALID_CAM;
                    double fp[3];
                    fp[0]=INVALID_CAM;
                    fp[1]=INVALID_CAM;
                    fp[2]=INVALID_CAM;
                    sstmpck >> tmpck;
                    //		      std::clog<<tmpck<<std::endl<<tmpck.compare("NULL")<<std::endl;
                    if(tmpck.compare("NULL")!=0) sstmp >> zoom;
                    sstmpck >> tmpck;
                    if(tmpck.compare("NULL")!=0) sstmp >> cp[0];
                    sstmpck >> tmpck;
                    if(tmpck.compare("NULL")!=0) sstmp >> cp[1];
                    sstmpck >> tmpck;
                    if(tmpck.compare("NULL")!=0)  sstmp >> cp[2];
                    sstmpck >> tmpck;
                    if(tmpck.compare("NULL")!=0) sstmp >> fp[0];
                    sstmpck >> tmpck;
                    if(tmpck.compare("NULL")!=0) sstmp >> fp[1];
                    sstmpck >> tmpck;
                    if(tmpck.compare("NULL")!=0) sstmp >> fp[2];
                    sstmpck >> tmpck;
                    if(tmpck.compare("NULL")!=0) sstmp >> camRoll;
                    if(zoom==INVALID_READ)
                    {
                        std::cerr<<"Invalid line for cycle file , ignored"<<std::endl;
                        break;
                    }
                    
                    m_cycleAzimuth.push_back(0);
                    m_cycleElevation.push_back(0);
                    m_cycleZoom.push_back(zoom);
                    m_cyclecamPos0.push_back(cp[0]);
                    m_cyclecamPos1.push_back(cp[1]);
                    m_cyclecamPos2.push_back(cp[2]);
                    m_cyclecamFp0.push_back(fp[0]);
                    m_cyclecamFp1.push_back(fp[1]);
                    m_cyclecamFp2.push_back(fp[2]);
                    m_cyclecamRoll.push_back(camRoll);
                    break;
                }
                case 10:
                {
                    m_vServer.setCameraFocalPoint=true;
                    m_vServer.setCameraPos=true;
                    double camRoll=INVALID_CAM;
                    double cp[3];
                    cp[0]=INVALID_CAM;
                    cp[1]=INVALID_CAM;
                    cp[2]=INVALID_CAM;
                    double fp[3];
                    fp[0]=INVALID_CAM;
                    fp[1]=INVALID_CAM;
                    fp[2]=INVALID_CAM;
                    sstmpck >> tmpck;
                    if(tmpck.compare("NULL")!=0) sstmp >> azimuth;
                    sstmpck >> tmpck;
                    if(tmpck.compare("NULL")!=0) sstmp >> elevation;
                    sstmpck >> tmpck;
                    if(tmpck.compare("NULL")!=0) sstmp >> zoom;
                    sstmpck >> tmpck;
                    if(tmpck.compare("NULL")!=0) sstmp >> cp[0];
                    sstmpck >> tmpck;
                    if(tmpck.compare("NULL")!=0) sstmp >> cp[1];
                    sstmpck >> tmpck;
                    if(tmpck.compare("NULL")!=0)  sstmp >> cp[2];
                    sstmpck >> tmpck;
                    if(tmpck.compare("NULL")!=0) sstmp >> fp[0];
                    sstmpck >> tmpck;
                    if(tmpck.compare("NULL")!=0) sstmp >> fp[1];
                    sstmpck >> tmpck;
                    if(tmpck.compare("NULL")!=0) sstmp >> fp[2];
                    sstmpck >> tmpck;
                    if(tmpck.compare("NULL")!=0) sstmp >> camRoll;
                    if(azimuth==INVALID_READ || elevation==INVALID_READ || zoom==INVALID_READ)
                    {
                        std::cerr<<"Invalid line for cycle file , ignored"<<std::endl;
                        break;
                    }		      
                    m_cycleAzimuth.push_back(azimuth);
                    m_cycleElevation.push_back(elevation);
                    m_cycleZoom.push_back(zoom);
                    m_cyclecamPos0.push_back(cp[0]);
                    m_cyclecamPos1.push_back(cp[1]);
                    m_cyclecamPos2.push_back(cp[2]);
                    m_cyclecamFp0.push_back(fp[0]);
                    m_cyclecamFp1.push_back(fp[1]);
                    m_cyclecamFp2.push_back(fp[2]);
                    m_cyclecamRoll.push_back(camRoll);
                    break;
                }
                    
            }//switch
        }
        
        
        if(m_cycleAzimuth.size() <=0)
        {		
            std::cerr<<"Invalid values in cycle file."<<std::endl;
            return ;
        }
	}
}
//----------------------------
bool OptionsSetter::setInternalData(VisIVOViewer *env)
//---------------------------
{
    m_vServer.path="NoFileApi";
    m_vServer.tableData= new float*[env->nCols];
    m_vServer.goodAllocation=true;
    m_vServer.nRows=env->nRows;
    m_vServer.nCols=env->nCols;
    
    if(env->setatt[VV_SET_VOLUME-VV_PARAM] !=1 &&  env->nCols<3)
    {
        std::cerr<<"Invalid 3D data table. Check for SET_VOLUME OPTION"<<std::endl;
        return false;
    }  
    m_vServer.comp[0]=env->comp[0];
    m_vServer.comp[1]=env->comp[1];
    m_vServer.comp[2]=env->comp[2];
    m_vServer.size[0]=env->size[0];
    m_vServer.size[1]=env->size[1];
    m_vServer.size[2]=env->size[2];
    for(int i=0;i<m_vServer.nCols;i++)
    {  
        m_vServer.tableData[i]=env->internalPointers[i];
        std::string sTmp(env->internalNames[i]);
        m_vServer.columns.insert(make_pair(sTmp,i));
    }
    m_vServer.internalData=true; //QUI
    return true;
}
//----------------------------
bool OptionsSetter::internalData()
//---------------------------
{
    if(m_vServer.internalData) return true; 
    return false; 
}


//---------------------------------------------------------------------
void OptionsSetter::writeHistory()
//---------------------------------------------------------------------
{
    if (m_vServer.m_historyEnabled) {
        viewParameter["file"]=m_vServer.path;
        
        new HistoryXmlWriter(m_vServer.m_historyFile.c_str(),viewParameter,outFilename);
        
    }
    
}
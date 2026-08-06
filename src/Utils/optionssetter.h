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

#ifndef OPTIONSSETTER_H
#define OPTIONSSETTER_H

#include <string>
#include <vector>
#include <map>
#include "visivodef.h"
#include "vstable.h"
#include "vstablemem.h"

struct SplotchCamera
{
	double position[3];
	double lookat[3];
	double sky[3];
	double roll;
	double fov;
	double clip[2];
};

struct VisIVOServerOptions
{
 
  int x,y,z;//!colums number associated to each axes 
  int numImage;//! id number of image
  int  numImageToLoad; //!number of image to generete
  int vx,vy,vz;//!colums number associated to each vectore axes
  int spr,spI,spC1,spC2,spC3;//!colums number associated used in splotch
  int hsml,spColor,spIntensity;//!colums number associated used in splotch
  int nRows, nCols;//!numeber of rows and colums 
  int nRadius, nHeight,nColorScalar,nVRenderingField ,nSlicePlane,
  nIsosurfaceField;//!colums number associated to radius scaling,hiegth scaling, lut, rendering, slice and isosurface 
  int nGlyphs,nColorTable;//!set the number of glyphs and luttable 
  int nSliceField;//!set the number of plane that the user want for slice
  
  double elevation,zoom,azimuth;//! value for camera position default il 0, 0,1
  bool setCameraPos, setCameraFocalPoint,setCameraRoll;
  double cameraPos[3],cameraFocalPoint[3],cameraRoll;//! value for camera position default il 0, 0,1
  double *cameraPosPrev,*cameraFocalPointPrev,*cameraRollPrev;//! value for camera position default il 0, 0,1
  bool cycle;//! boolean value for cycle in azimuth zooming and elevation
  std::string cycleFile; //! set if the system if little or big endian
  int cycleOffset; //! set if the system if little or big endian
  int cycleSkipFrom; //! skip the first # lines from the begin in cycle file
  int cycleSkipTo; //! read up to # line number in cycle file
 double opacity;//!value of opacity dafault is 0.666
  double radius,height ;//! value of radius and hieght for glyphs. the usere can use this with scaling or not. default is 1 for both 
  double size[3];//! cell size of volume
  double comp[3];//! resolution of volume 
  int slicePosition;//! value for the position of slice 
  double isosurfaceValue;//! value for isosurface
  
  bool needSwap;//! set if the colum need swap or not
  bool splotch;//! set if splotch is required
  
  std::string systemEndianism; //! set if the system if little or big endian
  std::string path;//! absolute path with filename of file that you want visualize
  std::string splotchpar;//! absolute path with filename of file used for splotch
  std::string xField,yField,zField;//! name of x, y and z field
  std::string sprField,spIField,spC1Field,spC2Field,spC3Field;//! name of r,I,C1-3 fields of Splotch
  std::string colorScalar;//! name of field used for lut 
  float colorRangeFrom; //!set range FROM in volume to be used for the lut
  float colorRangeTo; //!set range TO in volume to be used for the lut
  bool isColorRangeFrom; //!set range FROM in volume to be used for the lut
  bool isColorRangeTo; //!set range TO in volume to be used for the lut
  bool stereo;
  std::string stereoMode;
  int stereoImg;
  float anaglyphsat;
  std::string anaglyphmask;
  std::string imageName;//! name of output image 
  std::string xVectorField,yVectorField,zVectorField;//! name of x, y, z and field for create a vector
  std::string hsmlField,spColorField,spIntensityField;//!colums names used in splotch
  std::string scale;//! if is yes there is the scaling of axes
  std::string noDefault;//! if is yes you have only one image
  std::string  savepar;//! if is yes save all parameters in a txt file
  std::string loadpar; //! if is yes load a txt file with all parameteres
  std::string scaleGlyphs;//! if is yes the user can select the scaling for radius and/or glyphs 
  std::string   vector;//! if is yes the user can create vector 
  std::string color;//! if is yes you can use lut
  std::string oneColor;//! set the color for points and isosurfaces
  std::string uselogscale;//! if is yes you have the lut with a log scale
  std::string volume;//! if isi yes you can create a volume
  std::string vRendering;//! if is yes can visualize a volume rendering
  bool vShadow;//! if  yes  visualize a shadow volume
  bool showBox; //! if no the box is not shown
  bool showLut; //! if no the color bar is not shown
  bool showAxes; //! if yes the axes are  shown
  bool wireframe; //! if yes the isosurface is   shown using a wireframe
  bool vectorLine;//! if yes vectors are displayed with lines
  double vectorScalingFactor;
  int vectorScale;
  std::string backColor;//! assumes values: white (default), yellow, red, green,  blu, cyane, violet, black
  std::string imageSize;//! assume values: small, medium (default), large
  std::string isoSmooth;//! assume values: none (default), low, medium, high
  std::string slice;//! if is yes can visualize a slice
  std::string autoRange;//! if is yes volume renderer uses autorange filter
  float autoRangeMin; //!force autorange minimum value
  float autoRangeMax; //!force autorange maximum value
  bool setAutoRangeMin, setAutoRangeMax;
  std::string isosurface;//! if is yes can visualize a isosurface 
  std::string radiusscalar,heightscalar;//! name of field the isd used for scale the glyphs by radius and/or heigth
  std::string glyphs;//! glyphs selected
  std::string colorTable;//! luttable selected
  bool  colorIsFile;
  std::string endian;//! if the file is big o little endian
  std::string dataType;//! =float
  std::string vRenderingField,sliceField,isosurfaceField;//! name of field selected for rendering , slice and isosurface 
  std::string slicePlane;//! plane selected for slice 
  bool sliceOrthoNormal;
  std::string slicePlaneNormal;//! normal to a slicer plane
  std::string slicePlanePoint;//! application point for slice  
  double slicePlaneNormalValue[3];//! normal point coordinate for a slicer plane
  double slicePlanePointValue[3];//! application point coordinate for a slicer plane

  std::vector<std::string> fieldNames;//! names of fields in the original VBT  

  bool dataRead;//! data read in tabled data (if goodAllocation)
  bool goodAllocation;//! data already read in tableData
  std::map<std::string, int> columns;//! data read from table in tableData
  float **tableData; //! data read from table
  std::vector<int> extPalId;
  std::vector<double> extPalR,extPalG,extPalB,extPalA;//! palette color
  double spPosition[3];//! Splotch camera
  double spLookat[3];//! Splotch camera
  double roll;//! Splotch camera
  double spSky[3];//! Splotch camera
  double fov;
  bool fovIsGiven;
  double spClip[2];//! Splotch camera

// vtk image creation
   std::string vtkImage; //! Enables vtk image creation
   std::string vtkImSize; //! String Array Image size (grid dim)
   int vtkImSizeValue[3]; //! Array Image size (grid dim)
   float vtkGausExp;  //! Gaussian exp, negative value
   float vtkSpacingFact; //! Radius of influence  suggested [1,2] 
   float vtkScale; //! Scale (amplification) factor
   float vtkEcc; //! Exccentricity
   bool internalData;
   bool useMemory;

// vtk image creation
  std::string mode;
  std::string scenario;
  float opacityTF[3];
  std::string VO;
  std::string lfnout;
  std::string se;

// clipping plane
    bool clipset;
    double cliprange[2];
    
    //history
    bool m_historyEnabled;
    std::string m_historyFile;
    
    //splotch recalc cam
    bool autocalcCam;

};

class OptionsSetter
{
  static const unsigned int MAX_ELEMENT_TO_LOAD;
  static const float INVALID_READ;
  static const double INVALID_DREAD;
  static const double INVALID_CAM;
  public:
   
    OptionsSetter ( );
    ~OptionsSetter ( );
    
    void writeHistory();

    int parseOption (const std::vector<std::string> arguments );
    bool setInternalData (VisIVOViewer *env);
    void showHelp();
    int images();
    int readData();
    bool internalData();
    VisIVOServerOptions returnOptions(){return m_vServer;};
    void setUseMemory(bool use){useInMemory=use;}
    bool getUseMemory(){return useInMemory;}
    void setTable(VSTable* t){table=t;}
  protected:
  
    VisIVOServerOptions m_vServer;
    SplotchCamera m_spCamera;
    
  private:

    bool setAxisScalarVectorAndVolumesFields ();
    void saveParameters();
    void loadParameters(std::string path);
    bool readHeader();
    bool readSplocthColumn();
    void setGlyphs();
    void setColorLut();
    void setNumberImages();
    void tryToSetDimension();
    void readCycleFile();
    std::vector<float> m_cycleAzimuth;
    std::vector<float> m_cycleElevation;
    std::vector<float> m_cycleZoom;
    std::vector<double> m_cyclecamPos0;
    std::vector<double> m_cyclecamPos1;
    std::vector<double> m_cyclecamPos2;
    std::vector<double> m_cyclecamFp0;
    std::vector<double> m_cyclecamFp1;
    std::vector<double> m_cyclecamFp2;
    std::vector<double> m_cyclecamRoll;
    std::vector<int> m_cyclePlane;
    std::vector<std::string> m_cyclePlanePointNormal;
    bool m_inputLfnGiven;
    bool fileIsAtTheEnd;
    bool useInMemory = false;
    VSTable *table;
    
    std::map<std::string,std::string> viewParameter;
    std::vector <std::string> outFilename;
};

#endif
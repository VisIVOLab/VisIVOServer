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
#include <cstdlib>
#include <cstring>

#include <chrono>
#include "pointspipe.h"

#include "visivoutils.h"
#include "luteditor.h"

#include "extendedglyph3d.h"

#include <sstream>
#include <algorithm>

#include "vtkSphereSource.h"
#include "vtkConeSource.h"
#include "vtkCylinderSource.h"
#include "vtkCubeSource.h"

#include "vtkCamera.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkLookupTable.h"

#include "vtkFloatArray.h"
#include "vtkCellArray.h"
#include"vtkGlyph3D.h"
#include "vtkScalarBarActor.h"
#include "vtkOutlineCornerFilter.h"
#include "vtkProperty.h"
#include <vtkImageData.h>

#include "vtkGenericRenderWindowInteractor.h" 
#include "vtkRenderWindow.h"
#include "vtkRenderer.h"
#include "vtkActor.h"
#include "vtkAxesActor.h"
#include "vtkPolyDataWriter.h"
#include "vtkImageHistogramStatistics.h"
#include "vtkImageThreshold.h"

//---------------------------------------------------------------------
PointsPipe::PointsPipe ( VisIVOServerOptions options)
//---------------------------------------------------------------------
{
  m_visOpt=options;
  constructVTK();
  m_glyphFilter   = ExtendedGlyph3D::New();
  m_glyph         = vtkGlyph3D::New();
  m_pConeActor    = vtkActor::New();
  m_polyData      = vtkPolyData::New();
  m_pConeMapper   = vtkPolyDataMapper::New();
}
//---------------------------------
PointsPipe::~PointsPipe()
//---------------------------------
{
  destroyVTK();
  if ( m_glyph!=0)
    m_glyph->Delete() ;
  if ( m_glyphFilter!=0)
    m_glyphFilter->Delete() ;
  if ( m_pConeMapper != 0 )
    m_pConeMapper->Delete();
  if ( m_pConeActor != 0 )
    m_pConeActor->Delete();
  if ( m_polyData!=0)
    m_polyData->Delete() ;


}
//---------------------------------
void PointsPipe::destroyAll()
//---------------------------------
{
  if ( m_glyph!=0)
    m_glyph->Delete() ;
  if ( m_glyphFilter!=0)
    m_glyphFilter->Delete() ;
  if ( m_pConeMapper != 0 )
    m_pConeMapper->Delete();
  if ( m_pConeActor != 0 )
    m_pConeActor->Delete();
  if ( m_polyData!=0)
    m_polyData->Delete() ;


}

vtkSmartPointer<vtkFloatArray> PointsPipe::createAxisArray(const std::string& name, unsigned long long int nRows) {
    vtkSmartPointer<vtkFloatArray> arr = vtkSmartPointer<vtkFloatArray>::New();
    arr->SetName(name.c_str());
    arr->SetNumberOfTuples(nRows);
    return arr;
}

int PointsPipe::getColumnIndex(const std::map<std::string, int>& columns, const std::string& name) {
    auto it = columns.find(name);
    if (it == columns.end()) {
        throw std::runtime_error("Column not found: " + name);
    }
    return it->second;
}


void PointsPipe::setupVertexCells() {
    int nPoints = m_polyData->GetNumberOfPoints();

    vtkSmartPointer<vtkCellArray> newVerts = vtkSmartPointer<vtkCellArray>::New();
    newVerts->EstimateSize(nPoints, 1);

    for (int i = 0; i < nPoints; ++i) {
        newVerts->InsertNextCell(1);
        newVerts->InsertCellPoint(i);
    }

    m_polyData->SetVerts(newVerts);
}

bool PointsPipe::assignColumnToArray(vtkFloatArray* array, const std::string& fieldName, std::ifstream* inFile, int fileColumnIndex)
{
    array->SetNumberOfTuples(m_visOpt.nRows);
    if (m_visOpt.dataRead && m_visOpt.goodAllocation)
    {
        int colIndex = getColumnIndex(m_visOpt.columns, fieldName);
        if (colIndex == -1) return false;

        for (int i = 0; i < m_visOpt.nRows; ++i)
            array->SetValue(i, m_visOpt.tableData[colIndex][i]);

        return true;
    }
    else if (useMemory)
    {
        int colId = memTable->getColId(fieldName);
        
        if (colId == -1) return false;

        std::vector<float> buffer(m_visOpt.nRows);
        memTable->getColumn(colId, buffer.data(), 0, m_visOpt.nRows - 1);
        for (int i = 0; i < m_visOpt.nRows; ++i)
            array->SetValue(i, buffer[i]);

        return true;
    }
    else  
    {
      const size_t chunkSize = 4096;

      std::vector<float> buffer(chunkSize);
      std::streamoff columnOffset = fileColumnIndex * static_cast<std::streamoff>(m_visOpt.nRows) * sizeof(float);
      inFile->seekg(columnOffset);
      size_t remaining = m_visOpt.nRows;
      size_t index = 0;
	    const auto start{std::chrono::steady_clock::now()}; 
      while (remaining > 0) {
          size_t currentReadSize = std::min(chunkSize, remaining);
          inFile->read(reinterpret_cast<char*>(buffer.data()), currentReadSize * sizeof(float));

          if (m_visOpt.needSwap) {
              for (size_t i = 0; i < currentReadSize; ++i) {
                  buffer[i] = floatSwap(reinterpret_cast<char*>(&buffer[i]));
              }
          }

          for (size_t i = 0; i < currentReadSize; ++i) {
              array->SetValue(index++, buffer[i]); 
          }

          remaining -= currentReadSize;
      }

	      const auto end{std::chrono::steady_clock::now()}; 
        const std::chrono::duration<double> elapsed_seconds{end - start}; 
        std::clog << "time for reading: " << elapsed_seconds.count() << std::endl; 
        return true;
    }

    return false;
}

void PointsPipe::applySingleColor(vtkActor* actor, const std::string& colorName)
{
    vtkSmartPointer<vtkProperty> prop = vtkSmartPointer<vtkProperty>::New();
    if (colorName == "yellow")       prop->SetColor(1, 1, 0);
    else if (colorName == "red")     prop->SetColor(1, 0, 0);
    else if (colorName == "green")   prop->SetColor(0, 1, 0);
    else if (colorName == "blu")     prop->SetColor(0, 0, 1);
    else if (colorName == "cyane")   prop->SetColor(0, 1, 1);
    else if (colorName == "violet")  prop->SetColor(1, 0, 1);
    else if (colorName == "black")   prop->SetColor(0, 0, 0);
    else                             prop->SetColor(1, 1, 1); // default white

    actor->SetProperty(prop);
}

void PointsPipe::setBackgroundColor(const std::string& colorName)
{
    double r = 0.0, g = 0.0, b = 0.0; // default black
    if (colorName == "yellow")       r = g = 1.0;
    else if (colorName == "red")     r = 1.0;
    else if (colorName == "green")   g = 1.0;
    else if (colorName == "blue")    b = 1.0;
    else if (colorName == "cyan")    g = b = 1.0;
    else if (colorName == "violet")  r = b = 1.0;
    else if (colorName == "white")   r = g = b = 1.0;

    m_pRenderer->SetBackground(r, g, b);
}

void PointsPipe::setRenderWindow(const std::string& size)
{
    if (size == "small")       m_pRenderWindow->SetSize(512, 365);
    else if (size == "large")  m_pRenderWindow->SetSize(1024, 731);
    else                       m_pRenderWindow->SetSize(792, 566);

    m_pRenderWindow->SetWindowName ("VisIVOServer View");
}

void PointsPipe::configureStereoRender()
{
    m_pRenderWindow->StereoRenderOn();

    if (m_visOpt.stereoMode == "RedBlue") {
        m_pRenderWindow->SetStereoTypeToRedBlue();
    }
    else if (m_visOpt.stereoMode == "Anaglyph") {
        m_pRenderWindow->SetStereoTypeToAnaglyph();
        m_pRenderWindow->SetAnaglyphColorSaturation(m_visOpt.anaglyphsat);

        std::string trimmedMask = trim(m_visOpt.anaglyphmask);
        std::stringstream anatmp(trimmedMask);
        int rightColorMask = 4, leftColorMask = 3;
        anatmp >> rightColorMask;
        anatmp >> leftColorMask;

        m_pRenderWindow->SetAnaglyphColorMask(rightColorMask, leftColorMask);
    }
    else {
        m_pRenderWindow->SetStereoTypeToCrystalEyes();
        if (m_visOpt.stereoImg == 0)
            m_pRenderWindow->SetStereoTypeToRight();
        else
            m_pRenderWindow->SetStereoTypeToLeft();
    }

    m_pRenderWindow->StereoUpdate();
}

int PointsPipe::getColorTableIndex(const std::string& color) {
    static const std::map<std::string, int> colorMap = {
        {"white", 22}, {"yellow", 19}, {"red", 24},
        {"green", 25}, {"blue", 26}, {"cyan", 20},
        {"violet", 21}, {"black", 23}
    };
    auto it = colorMap.find(color);
    return it != colorMap.end() ? it->second : 22;
}

bool PointsPipe::loadHeightScalars(std::ifstream& inFile) {
  auto heightArrays = vtkSmartPointer<vtkFloatArray>::New();
  heightArrays->SetName(m_visOpt.heightscalar.c_str());

  if (!assignColumnToArray(heightArrays, m_visOpt.heightscalar, &inFile, m_visOpt.nHeight))
    return false;

  m_polyData->GetPointData()->SetScalars(heightArrays);
  return true;
}

bool PointsPipe::shouldLoadHeightScalars() const {
  return m_visOpt.heightscalar != "none" &&
         m_visOpt.scaleGlyphs != "none" &&
         m_visOpt.nGlyphs != 0 &&
         m_visOpt.nGlyphs != 1;
}

bool PointsPipe::shouldLoadRadiusScalars() const {
  return m_visOpt.radiusscalar != "none" &&
         m_visOpt.scaleGlyphs != "none" &&
         m_visOpt.nGlyphs != 0;
}

bool PointsPipe::loadRadiusScalars(std::ifstream& inFile, vtkFloatArray* radiusArrays) {
  if (m_visOpt.color == "none") {
    m_visOpt.color = "yes";
    m_visOpt.nColorTable = getColorTableIndex(m_visOpt.oneColor);
    m_visOpt.colorScalar = m_visOpt.radiusscalar;
    m_visOpt.showLut = false;
  }

  radiusArrays->SetName(m_visOpt.radiusscalar.c_str());
  return assignColumnToArray(radiusArrays, m_visOpt.radiusscalar, &inFile, m_visOpt.nRadius);
}

vtkDataArray* PointsPipe::setAutorange(vtkFloatArray* arr){
  auto imageData = vtkSmartPointer<vtkImageData>::New();
  imageData->SetDimensions(m_visOpt.nRows, 1, 1);
  imageData->AllocateScalars(VTK_FLOAT, 1); 
  imageData->GetPointData()->SetScalars(arr);
  double percentileValue[2];
  vtkImageHistogramStatistics *percentile = vtkImageHistogramStatistics::New();
  percentile->SetInputData(imageData);
  percentile->SetAutoRangePercentiles(97, 99);
  percentile->Update();
  percentile->GetAutoRange(percentileValue);

  std::clog << std::fixed << std::setprecision(25)
          << "Autorange " << percentileValue[0] << " " << percentileValue[1] << std::endl;
  vtkImageThreshold *threshold = vtkImageThreshold::New();
  threshold->ThresholdByUpper(percentileValue[0]);
  threshold->SetOutValue(percentileValue[0]);
  threshold->ReplaceOutOn();
  threshold->SetInputData(imageData);
  threshold->Update();

  vtkImageThreshold *threshold2 = vtkImageThreshold::New();
  threshold2->ThresholdByLower(percentileValue[1]);
  threshold2->SetOutValue(percentileValue[1]);
  threshold2->ReplaceOutOn();
  threshold2->SetInputConnection(threshold->GetOutputPort());
  threshold2->Update();
  return threshold2->GetOutput()->GetPointData()->GetScalars();
}

bool PointsPipe::loadColorScalars(std::ifstream& inFile) {
    auto lutArrays = vtkSmartPointer<vtkFloatArray>::New();
    if (!assignColumnToArray(lutArrays, m_visOpt.colorScalar, &inFile, m_visOpt.nColorScalar))
        return false;
    
    //lutArrays = vtkFloatArray::SafeDownCast(setAutorange(lutArrays));
    lutArrays->SetName(m_visOpt.colorScalar.c_str());
    m_polyData->GetPointData()->SetScalars(lutArrays);
    
    double range[2];
    lutArrays->GetRange(range);
    if (range[0] <= 0)
        m_visOpt.uselogscale = "none";

    return true;
}

void PointsPipe::configureGlyphs() {
  if (m_visOpt.nGlyphs != 0)
    setGlyphs();

  if (m_visOpt.scaleGlyphs != "none")
    setScaling();

  if(m_visOpt.scaleGlyphs!="none"&& m_visOpt.nGlyphs!=0 ||( (m_visOpt.heightscalar!="none" ||  (m_visOpt.radiusscalar!="none" && m_visOpt.nGlyphs!=1)) ))
    m_glyph->ScalingOn();
  else
    m_glyph->ScalingOff();
}

void PointsPipe::configureActor() {
  m_pConeMapper->SetInputData(m_polyData);
  m_pConeActor->SetMapper(m_pConeMapper);

  if (m_visOpt.color == "none") {
    applySingleColor(m_pConeActor, m_visOpt.oneColor);
  }

  m_pRenderer->AddActor(m_pConeActor);

  m_visOpt.opacity = std::clamp(m_visOpt.opacity, 0.0, 1.0);
  m_pConeActor->GetProperty()->SetOpacity(m_visOpt.opacity);
}

void PointsPipe::finalizeRender(vtkSmartPointer<vtkGenericRenderWindowInteractor> inter) {
  m_pRenderWindow->SetInteractor(inter);
  m_pRenderWindow->Render();

  setCamera();

  double bounds[6] = {
    m_xRange[0], m_xRange[1],
    m_yRange[0], m_yRange[1],
    m_zRange[0], m_zRange[1]
  };

  if (m_visOpt.showAxes)
    setAxes(m_polyData, bounds);

  inter->Start();
  inter->ExitEvent();
}

//-----------------------------------------------------------------------------------
int PointsPipe::createPipe ()
//------------------------------------------------------------------------------------
{
  int i = 0;
  std::ifstream inFile;
  if(useMemory) m_visOpt.nRows = memTable->getNumberOfRows();
  vtkSmartPointer<vtkFloatArray> xAxis = createAxisArray(m_visOpt.xField, m_visOpt.nRows);
  vtkSmartPointer<vtkFloatArray> yAxis = createAxisArray(m_visOpt.yField, m_visOpt.nRows);
  vtkSmartPointer<vtkFloatArray> zAxis = createAxisArray(m_visOpt.zField, m_visOpt.nRows);
  vtkSmartPointer<vtkFloatArray> radiusArrays = vtkSmartPointer<vtkFloatArray>::New();

  if (!m_visOpt.dataRead && !useMemory) {
      inFile.open(m_visOpt.path.c_str(), std::ios::binary);
      if (!inFile.is_open()) return -1;
  }
  inFile.open(m_visOpt.path.c_str(), std::ios::binary);
  if (!assignColumnToArray(xAxis, m_visOpt.xField, &inFile, m_visOpt.x)) return -1;
  if (!assignColumnToArray(yAxis, m_visOpt.yField, &inFile, m_visOpt.y)) return -1;
  if (!assignColumnToArray(zAxis, m_visOpt.zField, &inFile, m_visOpt.z)) return -1;

  //Retrieving Range
  xAxis->GetRange(m_xRange);  //!minimum and maximum value
  yAxis->GetRange(m_yRange);  //!minimum and maximum value
  zAxis->GetRange(m_zRange);  //!minimum and maximum value

  //assigning points
  SetXYZ(xAxis,yAxis,zAxis);  
  m_polyData->SetPoints(m_points);

    
  // connect m_pRendererderer and m_pRendererder window and configure m_pRendererder window
  m_pRenderWindow->AddRenderer(m_pRenderer);
  
  setupVertexCells();
  
  float tmp[1];
    
  if (shouldLoadRadiusScalars()) {
    if (!loadRadiusScalars(inFile, radiusArrays))
      return -1;
    m_polyData->GetPointData()->SetScalars(radiusArrays);
  }
  //UPDATE
  if (m_visOpt.colorScalar != "none") {
    if (!loadColorScalars(inFile))
      return -1;
  } 
  //UPDATE
  if (shouldLoadHeightScalars())
  {
    if (!loadHeightScalars(inFile))
      return -1;
  }
 
  inFile.close();

  setBoundingBox (m_polyData);

  configureActor();

  if (m_visOpt.colorScalar!="none" && m_visOpt.color!="none")
    setLookupTable ();

  configureGlyphs();
    
  setBackgroundColor(m_visOpt.backColor);
  
  setRenderWindow(m_visOpt.imageSize);
  //Valutazione sull'ordine
  if (m_visOpt.stereo)
    configureStereoRender();
  
  //open view
  vtkSmartPointer<vtkGenericRenderWindowInteractor> inter = vtkSmartPointer<vtkGenericRenderWindowInteractor>::New();
  finalizeRender(inter);

  return 0;
}


//---------------------------------------------------------------------
void PointsPipe::setGlyphs ( )
//---------------------------------------------------------------------
{
  int max=1000;
  
  if ( m_visOpt.nRows<max )
  {    
    /* VTK9 migration
    m_glyph->SetInput (m_polyData );
    replaced
    m_glyph->SetInputData (m_polyData );

    */
    m_glyph->SetInputData (m_polyData );
    
    
    if (m_visOpt.scale=="yes")      
      m_glyph->SetScaleFactor ( 0.04 );
    
    else
      m_glyph->SetScaleFactor ( 2.5 );
    
    m_pConeMapper->SetInputConnection( m_glyph->GetOutputPort() );
               
       
    if (m_visOpt.nGlyphs==1)
    {
      m_sphere   = vtkSphereSource::New();
      setResolution ( );
      setRadius ();
      /* VTK9 migration
      m_glyph->SetSource ( m_sphere->GetOutput() );
      replaced
            m_glyph->SetSourceData ( m_sphere->GetOutput() );
*/
      m_glyph->SetSourceData ( m_sphere->GetOutput() );
      m_sphere->Delete();
    }
    
    else if (m_visOpt.nGlyphs==2)
    {
      m_cone   = vtkConeSource::New();
      setResolution ( ); 
      setRadius ();
            /* VTK9 migration
      m_glyph->SetSource ( m_cone->GetOutput() );
      replaced
      m_glyph->SetSourceData ( m_cone->GetOutput() );
*/
      m_glyph->SetSourceData ( m_cone->GetOutput() );
      m_cone->Delete();
    } 
     
    else if (m_visOpt.nGlyphs==3)
    {  
      m_cylinder   = vtkCylinderSource::New();
      setResolution ( ); 
      setRadius ();
                  /* VTK9 migration
      m_glyph->SetSource ( m_cylinder->GetOutput() ); 
      replaced
      m_glyph->SetSourceData ( m_cylinder->GetOutput() ); 
*/
      m_glyph->SetSourceData ( m_cylinder->GetOutput() ); 
      m_cylinder->Delete(); 
    }
    
    else if (m_visOpt.nGlyphs==4)
    {
      m_cube   = vtkCubeSource::New();
      setRadius (); 
                        /* VTK9 migration
      m_glyph->SetSource ( m_cube->GetOutput() );  
      replaced
      m_glyph->SetSourceData ( m_cube->GetOutput() );  
*/
      m_glyph->SetSourceData ( m_cube->GetOutput() );  
      m_cube->Delete();
    }
    
   
  }
  return ;
}


//---------------------------------------------------------------------
void PointsPipe::setLookupTable ()
//---------------------------------------------------------------------
{ 

  double b[2];
  m_polyData->GetPointData()->SetActiveScalars(m_visOpt.colorScalar.c_str());
   
  m_polyData->GetPointData()->GetScalars(m_visOpt.colorScalar.c_str())->GetRange(b);
  
  m_lut->SetTableRange(m_polyData->GetPointData()->GetScalars()->GetRange());
  m_lut->GetTableRange(b);
  if(m_visOpt.isColorRangeFrom) b[0]=m_visOpt.colorRangeFrom;
  if(m_visOpt.isColorRangeTo) b[1]=m_visOpt.colorRangeTo;
  if(b[1]<=b[0]) b[1]=b[0]+0.0001;
  m_lut->SetTableRange(b[0],b[1]);
  
  if(m_visOpt.uselogscale=="yes")
    m_lut->SetScaleToLog10();
  else
    m_lut->SetScaleToLinear();
  
  m_lut->Build();
  
  SelectLookTable(&m_visOpt, m_lut);

  m_pConeMapper->SetLookupTable(m_lut); 
  m_pConeMapper->SetScalarVisibility(1);
  m_pConeMapper->UseLookupTableScalarRangeOn();

  m_pConeActor->SetMapper(m_pConeMapper);

  if(m_visOpt.showLut)  colorBar();

}


//---------------------------------------------------------------------
    void PointsPipe::setRadius ()
//---------------------------------------------------------------------
{
  if (m_visOpt.nGlyphs==1)
    m_sphere->SetRadius ( m_visOpt.radius);
   
    
  else if (m_visOpt.nGlyphs==2)
  {
    
    m_cone->SetRadius ( m_visOpt.radius );
    m_cone->SetHeight (m_visOpt.height );
  } 
     
  else if (m_visOpt.nGlyphs==3)
  {  
    
    m_cylinder->SetRadius (m_visOpt.radius ); 
    m_cylinder->SetHeight ( m_visOpt.height ); 
  }
  else if (m_visOpt.nGlyphs==4)
  {
    m_cube->SetXLength ( m_visOpt.radius );    
    m_cube->SetYLength ( m_visOpt.height );   
    m_cube->SetZLength ( 1 );  
    
  }
}
    
//---------------------------------------------------------------------
void PointsPipe::setResolution ()
//---------------------------------------------------------------------
{
  if (m_visOpt.nGlyphs==1)
  {
    m_sphere->SetPhiResolution ( 10 );
    m_sphere->SetThetaResolution ( 20 );
  }
    
  else if (m_visOpt.nGlyphs==2)
    m_cone->SetResolution ( 10 );
    
  else if (m_visOpt.nGlyphs==3)
    m_cylinder->SetResolution ( 10); 
   
   
}


//-------------------------------------------------------------------------
bool PointsPipe::SetXYZ(vtkFloatArray *xField, vtkFloatArray *yField, vtkFloatArray *zField  )
//-------------------------------------------------------------------------
{
  double scalingFactors[3];
  scalingFactors[0]=scalingFactors[1]=scalingFactors[2]=0;
         
  m_points=vtkPoints::New();
  m_points->SetNumberOfPoints(m_visOpt.nRows);
    
  
  if(xField->GetNumberOfComponents() != yField->GetNumberOfComponents())
  {
    if(zField && (xField->GetNumberOfComponents() != zField->GetNumberOfComponents() \
       || yField->GetNumberOfComponents() != zField->GetNumberOfComponents()))
    {
      return false;
    }
    return false; // component mismatch, do nothing
  }
  
  
  if(m_visOpt.scale=="yes")
  {
  
    double size = 0;

 
    size = (m_xRange[1] - m_xRange[0] != 0 ? m_xRange[1] - m_xRange[0] : m_xRange[1]);
    scalingFactors[0] = size * 0.1;


   
    size = (m_yRange[1] - m_yRange[0] != 0 ? m_yRange[1] - m_yRange[0] : m_yRange[1]);
    scalingFactors[1] = size * 0.1;

    
    size = (m_zRange[1] - m_zRange[0] != 0 ? m_zRange[1] - m_zRange[0] : m_zRange[1]);
    scalingFactors[2] = size * 0.1;
  }

  double scalingFactorsInv[3];

  int i = 0;
  for(i = 0; i < 3; i++)
    scalingFactorsInv[i] = ((scalingFactors && scalingFactors[i] != 0) ? 1/scalingFactors[i] : 0);

  // Set the points data
  if(m_visOpt.scale=="yes")
  {
    for(i = 0; i < m_visOpt.nRows; i++)
    {
      float inPoint[3];
      float outPoint[3];
      inPoint[0] = outPoint[0] = xField->GetValue(i) * scalingFactorsInv[0];
      inPoint[1] = outPoint[1] = yField->GetValue(i) * scalingFactorsInv[1];
      inPoint[2] = outPoint[2] = zField->GetValue(i) * scalingFactorsInv[2];

      m_points->SetPoint(i,outPoint);
    }
  }
  else
    for(i = 0; i < m_visOpt.nRows; i++)
  {
    float outPoint[3];
    
    outPoint[0] = xField->GetValue(i) ;
    outPoint[1] = yField->GetValue(i) ;
    outPoint[2] = zField->GetValue(i) ;

    m_points->SetPoint(i,outPoint);
  }
  
  return true;
}

//---------------------------------------------------------------------
void PointsPipe::setScaling ()
//---------------------------------------------------------------------
{
  m_glyphFilter->SetUseSecondScalar(true);
  m_glyphFilter->SetUseThirdScalar(true);
  
  m_glyphFilter->SetScaling(1);
  
  if( m_visOpt.heightscalar!="none" && m_visOpt.scaleGlyphs!="none" && m_visOpt.nGlyphs!=0 && m_visOpt.nGlyphs!=1)
    m_glyphFilter->SetInputScalarsSelectionY(m_visOpt.heightscalar.c_str());
      
  if( m_visOpt.radiusscalar!="none" && m_visOpt.scaleGlyphs!="none" && m_visOpt.nGlyphs!=0)  
    m_glyphFilter->SetInputScalarsSelectionXZ(m_visOpt.heightscalar.c_str());
 
  
  if( m_visOpt.nGlyphs!=0)
    m_glyphFilter->SetScaleModeToScaleByScalar();
  else 
    m_glyphFilter->ScalarVisibilityOff();
    

}
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

#ifndef PIPE_H
#define PIPE_H

#include "vtkPolyDataMapper.h"


#include "optionssetter.h"
#include "vstable.h"
#include "vstablemem.h"

class vtkRenderer;
class vtkRenderWindow;
class vtkCamera;
class vtkLookupTable;



class Pipe
{
  static const double INVALID_CAM;
 
  public:

    //void saveImageAsPng(int num );
    std::string saveImageAsPng(int num );

    virtual  int createPipe();
    virtual  void destroyAll(){};
    virtual  bool readData();
    virtual  int getCamera(SplotchCamera *splCamera);
    void setUseMemory ( bool use ){useMemory = use;};
    void setMemTable ( VSTable *table ){memTable = table;};
   
    
  protected:
    
    void setCamera (SplotchCamera *splCamera=NULL);
    void constructVTK();
     void destroyVTK();
     void setBoundingBox ( vtkDataObject *data );
     virtual void colorBar ();
     virtual  void setAxes(vtkDataSet *data,double *bounds);
    
    VisIVOServerOptions m_visOpt;
    
     vtkCamera          *m_camera;
     
    vtkRenderer       *m_pRenderer;
    vtkRenderWindow   *m_pRenderWindow;
    vtkLookupTable      *m_lut;
    VSTable* memTable;
    bool useMemory = false;

};

#endif
#ifndef VTKSOURCE_h
#define VTKSOURCE_h

#include "abstractsource.h"

#include <string>
#include <vector>
#include <fstream>

class VTKSource : public AbstractSource
{

  enum FileType
  {
    structured_points,
    polydata,
    structured_grid,
    unstructured_grid,
    rectilinear_grid
  };

public:
  int readHeader();
  int readData();
  VTKSource();
  ~VTKSource() = default;

private:
  FileType type;
};

#endif
/***************************************************************************
 *   Copyright (C) 2008 by Gabriella Caniglia, Roberto Munzone *
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

#ifndef CHANGASOURCE_H
#define CHANGASOURCE_H

#include "abstractsource.h"

#include <vector>
#include <string>
#include <rpc/xdr.h>
#include <rpc/rpc.h>

struct header
{
  double time;
  int nbodies;
  int ndim;
  int nsph;
  int ndark;
  int nstar;
  int pad;
};

struct gas_particle
{
  float mass;
  float pos[3];
  float vel[3];
  float rho;
  float temp;
  float eps;
  float metals;
  float phi;
};

struct dark_particle
{
  float mass;
  float pos[3];
  float vel[3];
  float eps;
  float phi;
};

struct star_particle
{
  float mass;
  float pos[3];
  float vel[3];
  float metals;
  float tform;
  float eps;
  float phi;
};

class ChangaSource : public AbstractSource

{
public:
  int readHeader();
  int readData();
  ~ChangaSource();
  ChangaSource();

private:
  int xdr_header(struct header *, XDR);
  std::vector<std::string> m_fieldsNames;
  unsigned int npart_total[6];
  XDR xdrread;
  FILE *fpread;
  int nsph;
  int ndark;
  int nstar;
  char m_dataType, m_Endian;
};

#endif
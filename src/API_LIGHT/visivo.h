// Define C Interface
#ifndef VISIVOC_H
#define VISIVOC_H

#include "visivodef.h"
#include "vstable.h"

#ifdef __cplusplus
extern "C"
{
#endif

	int VV_Init(VisIVOViewer *env);
	int VV_Clean(VisIVOViewer *env);
	int VV_SetXYZ(VisIVOViewer *env, float *X, float *Y, float *Z, char *names, int nRows);
	int VV_SetVect(VisIVOViewer *env, float *VX, float *VY, float *VZ, char *names, int nRows);
	int VV_SetColorScalar(VisIVOViewer *env, float *color, char *name, int nRows);
	int VV_SetRadiusScalar(VisIVOViewer *env, float *radiusscalar, char *name, int nRows);
	int VV_SetHeightScalar(VisIVOViewer *env, float *heightscalar, char *name, int nRows);
	int VV_SetVolume(VisIVOViewer *env, float *vol, char *name, char *geo, char *size);
	int VV_SetGSliceScan(VisIVOViewer *env, float *fgslice, int *iglsice);
	int VV_EnableSplotch(VisIVOViewer *env);
	int VV_SetSplotchLabIntensity(VisIVOViewer *env, float *intensity, char *name, int nRows);
	int VV_SetSplotchLabColor(VisIVOViewer *env, float *color, char *name, int nRows);
	int VV_SetSplotchLabhsml(VisIVOViewer *env, float *hsml, char *name, int nRows);
	int VV_SetAtt(VisIVOViewer *env, int code, char *value);
	int VV_SetTable(VisIVOViewer *env, VSTable *table);
	int VV_SetTableFromImporter(VisIVOViewer* viewer, VisIVOImporter* importer, size_t tableIndex);
	int VV_SetCameraPath(VisIVOViewer *env, int type, float *camera, int zoomend, float *zsf, int framesec,
						 int length, float *campos, float *camfp, float *camroll, int fcycle);
	int VV_View(VisIVOViewer *env);

	int VI_Init(VisIVOImporter *env);
	int VI_Clean(VisIVOImporter *env);
	int VI_SetAtt(VisIVOImporter *env, int code, char *value);
	int VI_Import(VisIVOImporter *env);
	void VI_FreeImporter(VisIVOImporter* env);

	int VF_Init(VisIVOFilter *env);
	int VF_Clean(VisIVOFilter *env);
	int VF_SetAtt(VisIVOFilter *env, int code, char *value);
	int VF_Filter(VisIVOFilter *env);

	int VS_VBTMetaData(VBT *tab, int satistic, char *value);
	int VS_VBTPointers(VBT *tab, unsigned int *colList, unsigned int nOfCols,
					   long unsigned int fromRow, long unsigned int toRow, float **fArray);
	int VS_VBTAllColumn(VBT *tab, int idCol, float *oneCol);
	int VS_VBTPartialColumn(VBT *tab, int idCol, float *oneCol, int from, int to);
	int VS_DirectImage(VisIVOViewer *VVenv, VisIVOImporter *VIenv, float random);

	void VA_Init(VisIVOAsynchId *id);
	int VA_Wait(VisIVOAsynchId *id);
	int VA_GetError(VisIVOAsynchId *id);
	int VA_GetState(VisIVOAsynchId *id);
	int VA_DirectImage(VisIVOViewer *VVenv, VisIVOImporter *VIenv, float random, VisIVOAsynchId *id);
	int VA_Filter(VisIVOFilter *e, VisIVOAsynchId *id);
	int VA_Import(VisIVOImporter *e, VisIVOAsynchId *id);
	int VA_View(VisIVOViewer *e, VisIVOAsynchId *id);

#ifdef __cplusplus
}
#endif

#endif
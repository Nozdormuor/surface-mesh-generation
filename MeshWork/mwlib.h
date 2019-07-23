#ifndef FILE_MESHWORK
#define FILE_MESHWORK

/* *************************************************************************/
/* File:   mwlib.hpp                                                 */
/* Function:  principal function  for users                      */
/* Date:   25. jan. 2016                                                   */
/* Function:                                       */
/* *************************************************************************/
// 
namespace meshwork
{
	
	 void MW_STL_GenerateSurfaceMesh(STLGeometry *pgeom, Mesh *pmesh);

	 void MW_WriteSTLFormat (Mesh *ngmesh, const char *filname);

	
}


#endif
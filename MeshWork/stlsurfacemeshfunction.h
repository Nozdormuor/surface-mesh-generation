#ifndef FILE_STLSURFACEMESHFUNTION
#define FILE_STLSURFACEMESHFUNTION
/* *************************************************************************/
/* File:   stlsurfacemeshfuntion.hpp                                                 */
/* Function:  function for surface mesh generation                       */
/* Date:   26. jan. 2016                                                   */
/* Function:                                                                            */
/* *************************************************************************/
namespace meshwork
{
class STLSurfaceMeshing;
int MW_STLSurfaceMeshing (STLGeometry & geom,class Mesh & mesh);
void TransFromTripsToCGAL(const vector<Point<3> >& ap,const vector<STLTriangle> &tris,const vector<int> trips,STLSurfaceMesh *pSFM);
bool SetBoundPropOnGEdge(TriangulationInKdtree *fvkdta,STLSurfaceMesh *psfma,const vector<Point<3> >& ap,const vector<STLLine> &lines,const vector<int> &facelines);
void MeshingWithNewNodes(STLSurfaceMesh *pSFM,  TriangulationInKdtree *fvt,vector<Point<3>> &meshpoints,vector<int>& meshsurfaces,class Mesh & mesh);
VertexHandle Insert1PointInInteriorOfTriangle(Point<3> &pt,FaceHandle *pFH,STLSurfaceMesh &sfm,TriangulationInKdtree *fvt,bool flag);
VertexHandle Insert1PointOnEdge(Point<3> &pt,FaceHandle *pFH,Edge edge,STLSurfaceMesh &sfm,TriangulationInKdtree *fvt,bool flag);
void LocalImproveAfterInsert( VertexHandle &vt,STLSurfaceMesh &sfm, TriangulationInKdtree *fvt);
}
#endif
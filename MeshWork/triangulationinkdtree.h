#ifndef FILE_TRIANGULATIONINKDTREE
#define FILE_TRIANGULATIONINKDTREE

namespace meshwork
{
	class STLSurfaceMesh;

	class TriangulationInKdtree{

	public:
		TriangulationInKdtree(STLSurfaceMesh *pM);
		~TriangulationInKdtree();
	public:
		VertexHandle * InsertVertexHandle(Triangulation::VertexHandle vh)
		{
			Triangulation::VertexHandle *pVH=new (Triangulation::VertexHandle)(vh);
			fvtree->InsertVert((void *)pVH);
			return pVH;
		}
		FaceHandle * InsertFaceHandle(Triangulation::FaceHandle fh)
		{
			Triangulation::FaceHandle *pFH=new (Triangulation::FaceHandle)(fh);
			fvtree->InsertExinfo((void *)pFH,1);
			return pFH;
		}
		void RemoveFaceHandle(Triangulation::FaceHandle fh)
		{
			fvtree->DeleteExinfo(fh->Info().bkaddress,1);
			delete (Triangulation::VertexHandle *)fh->Info().bkaddress;
			fh->Info().bkaddress=0;
		}
		bool GetFacehandleFromPoint(Point<3> pp,Triangulation::FaceHandle &fh,double *pdist);
		bool GetVertexhandleFromPoint(Point<3> p,Triangulation::VertexHandle &vh);
		Kodtree * GetKdtree(void)
		{
			return fvtree;
		}
		STLSurfaceMesh *GetSurfremesh(void)
		{
			return pMesh;
		}
		static void pofvforcoordnodes3(Point<3> &p,void *pv);
	private:
		static const double epsoverlap;
		Kodtree *fvtree;
		STLSurfaceMesh *pMesh;
		void InsertSurfremesh();
		//	void wrapPointsUpasVerts(void  ** &vti);
		static bool ifexinfooverlapbox(void *info,int infotype,const Box<3> &bd,double eps);
		static bool ifexinfoshouldbeincell(void *info,int infotype,CellNode *cnode);
	};

}
#endif
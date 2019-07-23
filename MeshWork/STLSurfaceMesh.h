#ifndef FILE_STLSURFACEMESH
#define FILE_STLSURFACEMESH

/* *************************************************************************/
/* File:   STLSurfaceMeshing.hpp                                                 */
/* Function:  mesh surface by delaunay algorithm                        */
/* Date:   26. jan. 2016                                                   */
/* Class:  Edge2v,  STLSurfaceMeshing                                                         */
/* *************************************************************************/
//STLSurfaceMesh类的私有公有函数还没处理呢！！！！！！！

namespace meshwork
{

	typedef Triangulation::VertexHandle VertexHandle;
	typedef Triangulation::FaceHandle FaceHandle;

	struct Edge2v
	{
		VertexHandle va,vb;
		Edge2v()
		{
			va=VertexHandle();
			vb=VertexHandle();
		}
		Edge2v(VertexHandle v0,VertexHandle v1):va(v0),vb(v1){};
	};

	class STLSurfaceMesh
	{
	public:
		enum tdsvertexstatus 
		{
			initial=-4,removed,accept,newone
		};//顶点的状态
		Triangulation sft;//CGAL三角片类
		std::list<Edge2v> qedgeforimp;//要反转的边	
		void *ptonbd,*ptinque;
	public:
		STLSurfaceMesh()
		{
			ptonbd=new int;
			ptinque=new int;
		}
		~STLSurfaceMesh()
		{ 
			delete (int *)ptonbd; 
			delete (int *)ptinque;
		}


		void  Get3VertexPtOfTdsTriFromaIndex(const FaceHandle &tri,int i, Point<3> &pa,Point<3> &pb,Point<3> &pc);
		void  Get1VertexPtOfTdsTri(const FaceHandle &tri,int i,Point<3> &p);
		bool SetEdgeProperty(Edge2v &e2v,void * ptprop);
		void SetEdgeProperty(FaceHandle &tri, int index,void * ptprop);
		bool ChangeEdgeFormFrom2vToTn(Edge2v e2v, Edge &e);
		void ChangeEdgeFormFromTnTo2v(Edge e, Edge2v &e2v);
		void InitializeAllVertexFlag(void);
		void  SetFlagOfTdsVertex(const VertexHandle &v, int flag);
		void FindTheNearestEdgeOfTdsTriToPtWithSqdist(Point<3> pt,const FaceHandle &tri,Edge &edge,double &dist);
		double  SqDistPtToTdsEdge(Point<3> pt,const FaceHandle &tri,int i);
		void  Get2VertexPtOfTdsEdge(Edge edge,Point<3> &pa,Point<3> &pb);
		void  Insert1PtInInteriorOfTdsTri(Point<3> &pt,const FaceHandle &tri,VertexHandle &v);
		void MaintainTriPropertyAroundaVertex(VertexHandle &v);
		void  Insert1PtOnTdsEdge(double pt[3],Edge edge,VertexHandle &v);
		bool  IsTdsEdgeOnBoundary(Edge edge);
		void  Get2TdsVertexOfEdge(Edge edge,VertexHandle &va,VertexHandle &vb);
		bool CanTdsEdgeSwap(Edge edge,bool geomapp=false);
		void  GetSpecifiedNeighbVertexPtOfTdsTri(FaceHandle tri,int i,Point<3> &p);
		bool  IsImprovedAfterEdgeSwap(Edge edge);
		void  Get4VertexPtOf2NeighbTri(const FaceHandle &tri,int i,Point<3> p4t[4]);
		void  SelectTheDiagonalOf4Pt(Point<3> p4t[4],int &k,double &cos0,double &cos1);
		double  ComputeCosValueOfAngleInaTriangle(Point<3> pa,Point<3> pb,Point<3> pc);
		void SwapEdge(FaceHandle tri1,int index);
		void MaintainTriProperty(FaceHandle &tri);
		bool IsVertexRemain(const VertexHandle &v);
		bool CanVertexDel(VertexHandle v);
		bool DeleteOneOldNode(const VertexHandle &v, std::queue<VertexHandle> *pqvt=0);
		void  ReduceDegreeBySwap(const VertexHandle &v);
		bool  IsTriPairCorrect(const FaceHandle &tri1,const FaceHandle &tri2);
		void CompPropertyOfVertex(const VertexHandle &v,int &intsegcount,Edge2v &e2v);
		bool  ChooseAMergeTdsVertexFrom4Incidents(const VertexHandle &v,VertexHandle &vm);
		bool Can2TdsVertsJoin(const VertexHandle &va,const VertexHandle &vb);
		void  Merge2TdsVertexOf1TdsEdge(const VertexHandle &v,const VertexHandle &vm,FaceHandle &ta,FaceHandle &tb);
		void LocalImpAfter1VtDel(std::list<Edge2v> &queueofedge, std::queue<VertexHandle> &quevt);
		void PushFreeEdgeInque(FaceHandle tri,int ind,std::list<Edge2v> &queueofedge);
		void	ImpSurfMeshFromEdgesInQue(std::list<Edge2v> &queueofedge);
		void UpdateQueueAndSwapStatusOfEdgesInQueue(list<Edge2v> &queueofedge,FaceHandle &tria, FaceHandle &trib);
		void  ImpSurfaceMesh(std::list<Edge2v> &queueofedge);
	};
	extern bool skipBeInfinite(VertexHandle &vt,Triangulation &sft);
	extern bool skipBeInfinite(FaceHandle &fa,Triangulation &sft);
	extern bool skipBeInfinite(Edge &ed,Triangulation &sft);
}
#endif
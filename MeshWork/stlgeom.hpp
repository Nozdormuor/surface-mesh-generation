#ifndef FILE_STLGEOM
#define FILE_STLGEOM

/**************************************************************************/
/* File:   stlgeom.hpp                                                     */
/* Date:   09. jan. 2016                                                   */
/**************************************************************************/

/**
STL Geometry


Terminology:

Point ... coordinates of STL triangles
Triangle  (short Trig)  STL triangle
TopEdge .... edge in topology, boundary of STL triangles (many)
Edge .... Edges which will occur in the mesh (confirmed edges, less)
*/

namespace meshwork
{

	class STLGeometry : public STLTopology,  public MeshWorkGeometry
	{
		vector <int> vertices;

		vector<STLEdge> edges;

		vector<int> markededges;//if have been used = 1, otherwise=0;

		Table<int> edgesperpoint;
		//stlgeom line
		int linecnt;
		vector<STLLine> lines;
		vector<int> linepoints;//if = 1, otherwise=0;
		vector<int> lineendpoints; //per geometrypoint, 1 = is endpoint; 0 = no endpoint,

		//stlgeom face 
		int facecnt; 
		vector<STLFace> faces;

		//spiralpoints:
		vector<int> spiralpoints;

	public:
		STLGeometry();
		virtual ~STLGeometry();

		void Load (istream & ist);
		void LoadBinary (istream & ist);

		STLEdge &GetEdge(int nr) 
		{
			return edges[nr];
		}
        int GetNE() 
		{
			return edges.size();
		}

		STLLine &GetLine(int nr) 
		{
			return lines[nr];
		}
		int GetNL() 
		{
			return lines.size();
		}
		vector<STLLine> &GetLines()
		{
			return lines;
		}

		STLFace &GetFace(int nr) 
		{
			return faces[nr];
		}
		int GetNF() 
		{
			return facecnt;
		}
		//get NO edges per point
		int GetEPPSize() const 
		{
			return edgesperpoint.Size();
		}
		int GetNEPP(int pn) 
		{
		  if (edgesperpoint.Size() == 0) 
		  {
			  BuildEdgesPerPoint();
		  }
		  return edgesperpoint.RowSize(pn);
		};
		int GetEdgePP(int pn, int vi)
		{
		  if (edgesperpoint.Size() == 0) 
		  {
			  BuildEdgesPerPoint();
		  }
		  return edgesperpoint.Get(pn,vi);
		};
		void AddEdgePP(int pn, int vn) 
		{
			edgesperpoint.Add(pn,vn);
		}
		int IsEdge(int p1, int p2);
		int IsEdgeNum(int p1, int p2);

		virtual void InitSTLGeometry (const vector<STLReadTriangle> & readtrigs);
		void RecoverBody();

		void BuildEdges();
		void FindEdgesFromAngles();
		int AddEdge(const STLEdge& v) 
		{
			edges.push_back(v);
			return edges.size();
		}
		void BuildEdgesPerPoint();

		void BuildSTLLines();
		int AddLine(STLLine line) 
		{
			lines.push_back(line); return lines.size();
		}
		bool IsLineClosed(STLLine &line)
		{
			return line.EndP()==line.StartP();
		}

		void BuildVertex();
		int AddVertex(int p);

		void FindLinePoints();
		void FindLineEndPoints();
		int STLGeometry :: IsLineEndPoint(int pn) ;
		void FindSpiralPoint();
		int AddSpiralPoint(int n);

		void CalcFaceNums();
		void BuildSTLFace();
		void AddTrigToFace(int nface, int trig);
		void AddLineToFace(int nface, int line);

		void FixGeometry();
		void ClearSmallFaces(void);
		double AreaOfFace(STLFace &face);
		 int GetFlatNeighbFace(int face);
		 void UnionFace( int f1 , int f2 );
		 void ClearSmallInteriorEdges();

		void RestrictLocalH(class Mesh & mesh, double gh, double minh);

private:
		int FindAdjacentEdge  (STLEdge &edge, int pnum) const;

	};

}
#endif

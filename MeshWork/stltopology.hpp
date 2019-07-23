#ifndef FILE_STLTOPOLOGY
#define FILE_STLTOPOLOGY

/**************************************************************************/
/* File:   stltopology.hpp                                                */
/* Record STLgeometry topologic information                         */
/* Date:   12. jan. 2016                                                   */
/**************************************************************************/

/*
  The STLTopology contains topologic information as
  triangle->point, point->triangles, triangle->edge, 2-points->edge,...
*/
namespace meshwork
{

	class STLGeometry;

#define STLBASE 1

	class STLPointIndex
	{
		int i;
	public:
		STLPointIndex () { ; }
		STLPointIndex (int ai) : i(ai) { ; }
		STLPointIndex & operator= (const STLPointIndex & ai) { i = ai.i; return *this; }
		STLPointIndex & operator= (int ai) { i = ai; return *this; }
		operator int () const { return i; }
		STLPointIndex operator++ (int) { return i++; }
		STLPointIndex operator-- (int) { return i--; }
	};



	class STLTrigIndex
	{
		int i;
	public:
		STLTrigIndex () { ; }
		STLTrigIndex (int ai) : i(ai) { ; }
		STLTrigIndex & operator= (const STLTrigIndex & ai) { i = ai.i; return *this; }
		STLTrigIndex & operator= (int ai) { i = ai; return *this; }
		operator int () const { return i; }
		STLTrigIndex operator++ (int) { return i++; }
		STLTrigIndex operator-- (int) { return i--; }
	};





	// triangle structure for loading stl files
	class STLReadTriangle
	{
		Vec<3> normal;
		Point<3> pts[3];
	public:
		STLReadTriangle (const Point<3> * apts, const Vec<3> & anormal);
		STLReadTriangle () {};
		const Point<3> & operator[] (int i) const 
		{
			return pts[i]; 
		}
		const Vec<3> & Normal() const 
		{ 
			return normal; 
		}
	};

	/**
	Topology Triangle:
	A Triangle sharing 3 points, 3 edges , 3 or 6 trigs and 2 domains.
	*/
	class STLTriangle
	{
		// topology edges of triangle, edge[i] opposite to point[i]
		int topedges[3];
		// neighbour triangles, trig[i] opposite to point[i]
		int nbtrigs[2][3]; 
		// normalized stored normal vector ??
		Vec<3> normal;
		// point numbers of triangle
		int pts[3];
		// front-side and back-side domains
		int domains[2];


	public:

		Box<3> box;
		Point<3> center;
		double rad;
		int facenum;

		struct 
		{
			unsigned int toperror : 1;
		} flags;




		STLTriangle (const int * apts);
		STLTriangle () 
		{
			pts[0]=0;pts[1]=0;pts[2]=0;
		}

		int operator[] (int i) const 
		{
			return pts[i]; 
		}
		int & operator[] (int i) 
		{ 
			return pts[i]; 
		}

		int EdgeNum(int i) const 
		{ 
			return topedges[i];
		}
		int & EdgeNum(int i) 
		{ 
			return topedges[i]; 
		}

		int PNum(int i) const 
		{
			return pts[(i)]; 
		}
		int & PNum(int i) 
		{ 
			return pts[(i)]; 
		}

		int NBTrigNum(int i) const
		{ 
			return nbtrigs[0][i];
		}
		int & NBTrigNum(int i)
		{ 
			return nbtrigs[0][i]; 
		}
		int NBTrig (bool side, int i) const
		{
			return nbtrigs[side][i]; 
		}
		int & NBTrig (bool side, int i) 
		{ 
			return nbtrigs[side][i]; 
		}

		int Domain (bool side) const
		{ 
			return domains[side]; 
		}
		int & Domain (bool side) 
		{
			return domains[side]; 
		}

		int GetCw(int i)
		{
			return pts[(i+2)%3];
		}
		int GetCcw(int i) 
		{
			return pts[(i+1)%3];
		}

		int PIndex( int pointnum) 
		{
			for (int i=0; i<3; i++)
			{
				if(pointnum==pts[i])
					return i;
			}
			return -1;
		}

		int AnotherPointNum(int pointindex1, int pointindex2);
		int AnotherPointIndex(int pointnum1, int pointnum2);

		int IsCoEdge(STLTriangle trig);
		// consistently oriented neighbour:
		int IsNeighbourFrom(const STLTriangle& t) const;
		// opposite to consistently oriented neighbour:
		int IsWrongNeighbourFrom(const STLTriangle& t) const;

		// Stored normal vector, normalized
		void SetNormal (const Vec<3> & n);
		const Vec<3> & Normal () const 
		{ 
			return normal; 
		}

		int GetFaceNum() 
		{
			return facenum;
		}
		void SetFaceNum(int i)
		{
			facenum = i;
		}

		void GetNeighbourPoints(const STLTriangle& t, int& p1, int& p2) const;
		int GetNeighbourPointsAndOpposite(const STLTriangle& t, int& p1, int& p2, int& po) const;
	};

	ostream& operator<<(ostream& os, const STLTriangle& t);

	/**
	Topology Edge:
	Useful unside a face.
	A edges sharing more than 2 faces: trigs are undefined 
	*/
	class STLTopEdge 
	{
		int pts[2];  
		int trigs[2];  
		double cosangle;
		int status;  // excluded, confirmed, candidate, undefined
	public:
		STLTopEdge ();
		STLTopEdge (int p1, int p2, int trig1, int trig2);

		int operator[] (int i) const
		{ 
			return pts[i]; 
		}
		int & operator[] (int i)
		{
			return pts[i];
		}

		int PNum(int i) const 
		{ 
			return pts[i]; 
		}
		int & PNum(int i) 
		{
			return pts[i]; 
		}

		int TrigNum(int i) const 
		{
			return trigs[i]; 
		}
		int & TrigNum(int i)
		{
			return trigs[i]; 
		}

		void SetCosAngle (double ca)
		{ 
			cosangle = ca;
		}
		double CosAngle () const
		{
			return cosangle;
		}
		double Angle () const 
		{ 
			return acos (cosangle);
		}

		void SetStatus (int stat)
		{
			status = stat; 
		}
		int GetStatus () const
		{ 
			return status; 
		}
	};

	ostream& operator<<(ostream& os, const STLTopEdge& e);

	/**
	Topology geometry:
	*/
	class STLTopology
	{
		static double  geom_tol_fact;
	protected:
		vector<STLTriangle> trias;
		vector<STLTopEdge> topedges;
		vector<Point<3> > points;
		//normals belong to points!
		vector<Vec<3>> normals; 
		// mapping of node to trigs
		Table<int> trigsperpoint; 
		// mapping of node to edges
		Table<int> topedgesperpoint; 

		//  Box3dTree * searchtree; // ADT
		Point3dTree * pointtree;
		// searchtree for trigs and points
		Box<3> boundingbox;
		double pointtol;

	public:
		enum STL_GEOM_STATUS { STL_GOOD, STL_WARNING, STL_ERROR };

	protected:
		STL_GEOM_STATUS status;  
		bool topology_ok;

	public:
		STLTopology();
		virtual ~STLTopology();

		virtual void InitSTLGeometry (const vector<STLReadTriangle> & readtrigs);
		const Box<3> & GetBoundingBox () const
		{ 
			return boundingbox;
		}
		Box<3> & GetBoundingBox ()  
		{
			return boundingbox; 
		}

		int GetNP() const 
		{
			return points.size();
		}
		const Point<3> & GetPoint(int nr) const 
		{ 
			return points[nr]; 
		}
		int GetPointNum (const Point<3> & p);
		void SetPoint(int nr, const Point<3> & p) 
		{ 
			points[nr] = p;
		}
		const vector<Point<3> >& GetPoints() const
		{ 
			return points; 
		}
		const Point<3> & operator[] (STLPointIndex i) const 
		{ 
			return points[i];
		}
		Point<3> & operator[] (STLPointIndex i) 
		{ 
			return points[i]; 
		}
		const Vec<3> & GetNormal(int nr) const
		{
			return normals[nr];
		}
		void SetNormal(int nr, const Vec<3>& n)
		{
			normals[nr] = n;
		}

		int GetNT() const
		{
			return trias.size();
		}
		void AddTriangle(const STLTriangle& t);
		const STLTriangle & GetTriangle (int nr) const 
		{ 
			return trias[nr];
		}
		STLTriangle & GetTriangle (int nr)
		{ 
			return trias[nr]; 
		}
		const vector<STLTriangle>& GetTriangles() const
		{ 
			return trias; 
		}
		const STLTriangle & operator[] (STLTrigIndex i) const
		{
			return trias[i]; 
		}
		STLTriangle & operator[] (STLTrigIndex i) 
		{ 
			return trias[i]; 
		}

		int GetNTE() const
		{ 
			return topedges.size(); 
		}
		const STLTopEdge & GetTopEdge (int nr) const 
		{ 
			return topedges[nr]; 
		}
		STLTopEdge & GetTopEdge (int nr) 
		{ 
			return topedges[nr];
		}
		int GetTopEdgeNum (int pi1, int pi2) const;

		int NOTrigsPerPoint(int pn)
		{ 
			return trigsperpoint.RowSize(pn); 
		}
		int NTopEdgesPerPoint (int pn) const 
		{ 
			return topedgesperpoint.RowSize(pn);
		}
		int NeighbourTrig(int trig, int nr) const 
		{
			return trias[trig].NBTrigNum(nr); 
		}

	protected:
		int AddPoint(const Point<3> & p) 
		{ 
			points.push_back(p); 
			return (points.size()-1);
		}
		void FindTrigsOfPoints();
		void BuildSTLTopEdge();
		void FindEdgesOfPoint();
		void FindNeighbourTrigs();
		void FindEdgesOfTrig();
		void CalculateCosangleOfEdge();
		void CalculateVectorOfPoint();

	};

}
#endif

#ifndef FILE_STLBOUND
#define FILE_STLBOUND


/**************************************************************************/
/* File:   stlbound.hh                                                     */
/* Date:   09. jan. 2016                                                      */
/* Class:  STLEdge  STLLine  STLChart                            */
/**************************************************************************/

namespace meshwork
{

	class STLGeometry;
	class STLTopology;

	/* ******************************* STLEdge ******************************* */
	class STLEdge
	{
	public:
		int pts[2];
		int trigs[2]; //left and right trig

		STLEdge (const int * apts) 
		{
			pts[0] = apts[0]; pts[1] = apts[1];
		}
		STLEdge (int v1, int v2) 
		{
			pts[0] = v1; pts[1] = v2;
		}
		STLEdge () 
		{
			pts[0]=0;pts[1]=0;
		}
		int PNum(int i) const
		{
			return pts[(i)];
		}

		int LeftTrig() const 
		{
			return trigs[0];
		}
		int RightTrig() const
		{
			return trigs[1];
		}
		void SetLeftTrig(int i) 
		{
			trigs[0] = i;
		}
		void SetRightTrig(int i)
		{
			trigs[1] = i;
		}
	};

	enum STL_ED_STATUS { ED_EXCLUDED, ED_CONFIRMED, ED_CANDIDATE, ED_UNDEFINED };

	/* ******************************* STLLine ******************************* */
	//a line defined by several points (polyline)
	class STLLine
	{
	private:
		vector<int> pts;
		vector<double> dists;
		vector<int> lefttrigs;
		vector<int> righttrigs;
		int split;

	public:
		STLLine();
		void AddPoint(int i) 
		{
			pts.push_back(i);
		}
		void AddLeftTrig(int i) 
		{
			lefttrigs.push_back(i);
		}
		void AddRightTrig(int i)
		{
			righttrigs.push_back(i);
		}
		int PNum(int i) const 
		{
			return pts[i];
		}
		int NP() const 
		{
			return pts.size();
		}
		int GetNS() const;
		void GetSeg(int nr, int& p1, int& p2) const;
		double GetSegLen(const vector<Point<3> >& ap, int nr) const;
		double GetDist(int nr) const 
		{ 
			return dists[nr];
		}
		int GetLeftTrig(int nr) const
		{
			return lefttrigs[nr];
		}
		int GetRightTrig(int nr) const
		{
			return righttrigs[nr];
		}
		void GetBoundingBox (const vector<Point<3> > & ap, Box<3> & box) const;
		void AddDist (double dist) 
		{
			dists.push_back(dist); 
		}
		int StartP() const 
		{
			return pts[0];
		}
		int EndP() const 
		{
			return pts[pts.size()-1];
		}
		int ShouldSplit() const
		{
			return split;
		}  
		void Clear()
		{
			pts.clear();
			dists.clear();
			lefttrigs.clear();
			righttrigs.clear();
		}
		double GetLength(const vector<Point<3> >& ap) const;
		Point<3> GetPointInDist(const vector<Point<3> >& ap, double dist, int& index) const;
		STLLine Mesh(const vector<Point<3> >& ap, vector<Point<3>>& mp, double ghi,class Mesh& mesh) const;
	};


	/* ******************************* STLFace ******************************* */
	class STLFace
	{
	private:
		vector<int> facetriangles;
		vector<int> facelines;

	public:

		STLFace();
		vector<int> & Triangles()
		{
			return facetriangles;
		}
		vector<int> & Lines()
		{
			return facelines;
		}
		int GetTriangle(int nr)
		{
			return facetriangles[nr];
		}
		const int GetTriangle(int nr) const 
		{
			return facetriangles[nr];
		}
		int GetLine(int nr) 
		{
			return facelines[nr];
		}
		const int GetLine(int nr) const 
		{
			return facelines[nr];
		}
		int GetFaceLineNum()
		{
			return facelines.size();
		}
		const int GetFaceLineNum()const 
		{
			return facelines.size();
		}
		int GetFaceTriangleNum()
		{
			return facetriangles.size();
		}
		const int GetFaceTriangleNum()const 
		{
			return facetriangles.size();
		}
		int SetNodesOnSurface( Point<3> pa , Point<3> pb ,double ra, double rb,const vector<Point<3> >& ap,const vector<STLTriangle> &tris, Point<3>  p[2]  );
		int IsPointOnSurface( const vector<Point<3>> &pois , const vector<STLTriangle> &trips ,Point<3> p);
		int IntersectionOfCircleAndSurface( const vector<Point<3> >& pois,const vector<STLTriangle> &trips ,Point<3> p1 , Point<3> p2 , double r , Point<3> p[2] );
		vector<int> Mesh(const vector<Point<3> >& ap,const vector<STLTriangle> &tris,vector<Point<3>>& mp,const vector<STLLine> &line, double ghi,class Mesh& mesh) ;
	};

}

#endif

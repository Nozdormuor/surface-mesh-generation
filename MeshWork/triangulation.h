#ifndef FILE_TRIANGULATION
#define FILE_TRIANGULATION

/* *************************************************************************/
/* File:   triangulation.hpp                                                 */
/* Function:  topological operation of triangulation                         */
/* Date:   21. jan. 2016                                                   */
/* Class: Vertex Face Triangulation                                      */
/* *************************************************************************/
namespace meshwork
{

	class Vertex;
	class Face;
	class Triangulation;
	class FaceCirculator
	{
	public:
		list<Face>::iterator fi;
		int index;
	public:
		FaceCirculator(){;}
		FaceCirculator(list<Face>::iterator f,int i){fi=f;index=i;}
		FaceCirculator & operator++();
		FaceCirculator & operator--();
		FaceCirculator  operator++(int);
		FaceCirculator  operator--(int);
	};
	inline bool operator == (const FaceCirculator &lhs, const FaceCirculator &rhs)
	{
		return(lhs.fi==rhs.fi);
	}
	inline bool operator != (const FaceCirculator &lhs, const FaceCirculator &rhs)
	{
		return !(lhs.fi==rhs.fi);
	}
	class VertexCirculator
	{
	public:
		list<Vertex>::iterator  v0;
		list<Vertex>::iterator  v;
		list<Face>::iterator f;
	public:
		VertexCirculator(){;}
		VertexCirculator(list<Vertex>::iterator  v0,list<Vertex>::iterator  v,list<Face>::iterator f){this->v0=v0;this->v=v;this->f=f;}
		VertexCirculator & operator++();
		VertexCirculator & operator--();
		VertexCirculator  operator++(int);
		VertexCirculator  operator--(int);
	};
	inline bool operator == (const VertexCirculator &lhs, const VertexCirculator &rhs)
	{
		return(lhs.v==rhs.v);
	}
	inline bool operator != (const VertexCirculator &lhs, const VertexCirculator &rhs)
	{
		return !(lhs.v==rhs.v);
	}
	//enum ED_STATUS { ED_CANFLIP, ED_CANOTFLIP };

	struct VertexInfomation
	{
		VertexInfomation()
		{
			haveR=false;
			radius=-1;
			index=-1;
			status=-4; 
			pta=ptb=0;
			bkvaddress=0;
		}
		int status;
		double radius;
		void *pta,*ptb,*bkvaddress;
		bool haveR;
		int index;
	};
	typedef pair<list<Face>::iterator, int> Edge;
	struct Edge2v;
	struct TriInfomation
	{
		TriInfomation()
		{
			eproperty[0]=eproperty[1]=eproperty[2]=0;
			tflag=0;
			bkaddress=0;
		}
		void *eproperty[3];
		void *bkaddress;
		int tflag;
	};
	/*************************************Vextex*************************************/
	class Vertex
	{
	private:
		Point<3> point;
		int label;
		VertexInfomation infos;
		list<Face>::iterator incident_face;
	public:
		Vertex();
		Vertex(Point<3> poi);
		Vertex(Point<3> poi,int lab);
		Vertex(double x, double y, double z);
		Vertex(double x, double y, double z,int lab);
		Vertex & operator= (const Vertex & v2)
		{ 
			point=v2.point;
			label=v2.label;
			incident_face=v2.incident_face;
			return *this; 
		}
		void SetPoint(Point<3> pp){point=pp;}
		void SetPoint(double x, double y, double z)
		{
			Point<3> poi(x,y,z);
			point=poi;
		}
		void SetLabel(int lab)
		{
			label=lab;
		}
		int &GetLabel()
		{
			return label;
		}
		Point<3> &GetPoint()
		{
			return point;
		}
		double & operator() (int i) 
		{ 
			return point[i]; 
		}
		const double & operator() (int i) const
		{
			return point[i]; 
		}
		double & operator[] (int i) 
		{ 
			return point[i]; 
		}
		const double & operator[] (int i) const 
		{
			return point[i];
		}

		bool operator==(const Vertex & vv);

		FaceCirculator  IncidentFaces();
		VertexCirculator IncidentVertices();
		list<Face>::iterator GetFace()
		{
			return incident_face;
		}
		int  Degree();
		int SetFace(const list<Face>::iterator fh);
		bool IsIncidentFace(const list<Face>::iterator fh);
		struct VertexInfomation & Info(){return infos;}
	};

	/*************************************Face*************************************/

	class Face
	{
	private:
		list<Vertex>::iterator vertex[3];
		list<Face>::iterator neighbor[3];
		int edstatus[3];
		TriInfomation infos;
	public:
		Face();
		Face(list<Vertex>::iterator v1,list<Vertex>::iterator v2, list<Vertex>::iterator v3);
		Face(list<Face>::iterator f1, list<Face>::iterator f2, list<Face>::iterator f3);
		Face(list<Vertex>::iterator v1,list<Vertex>::iterator v2, list<Vertex>::iterator v3,list<Face>::iterator f1, list<Face>::iterator f2, list<Face>::iterator f3);
		Face & operator= (const Face & f2);

		bool operator==(const Face& f2);

		list<Vertex>::iterator & operator() (int i) { return vertex[i]; }
		const list<Vertex>::iterator & operator() (int i) const { return vertex[i]; }
		list<Vertex>::iterator & operator[] (int i) { return vertex[i]; }
		const list<Vertex>::iterator & operator[] (int i) const { return vertex[i]; }

		void SetVertex(list<Vertex>::iterator vh,int index) {vertex[index]=vh;}
		void SetNeighbor(list<Face>::iterator fh, int index){neighbor[index]=fh;}
		void SetEdgeStatus(int st, int index){edstatus[index]=st;}

		list<Face>::iterator &GetNeighbor(const int index){return neighbor[index];}
		list<Vertex>::iterator &GetVertex(const int index){return vertex[index];}
		int &EdgeStatus(const int index){return edstatus[index];}

		bool IsVertex(const list<Vertex>::iterator vh);
		bool IsNeighbor(const list<Face>::iterator fh);
		int  Index_f(const list<Face>::iterator fh);
		int  Index_v(const list<Vertex>::iterator vh);
		int Index(const list<Face>::iterator fh);
		int Index(const list<Vertex>::iterator vh);
		int Index(Vertex v);
		struct TriInfomation & Info(){return infos;}

	};

	/*************************************Triangulation*************************************/
	class Triangulation
	{
	public:
		typedef list<Vertex>::iterator VertexHandle;
		typedef list<Face>::iterator FaceHandle;
		typedef list<Vertex>::iterator AllVerticesIterator;
		typedef list<Face>::iterator AllFacesIterator;
		//typedef pair<FaceHandle, int> Edge;
		
	private:
		list<Vertex> vertexs;
		list<Face> faces;
		VertexHandle infinitevertex;
		bool hasinfinite;
	public:
		Triangulation();

		VertexHandle CreateVertex();
		VertexHandle CreateVertex(Point<3> poi);
		VertexHandle CreateVertex(double x, double y, double z);
		VertexHandle CreateVertex(Vertex v);
		FaceHandle CreateFace();
		FaceHandle CreateFace(Face fa);
		FaceHandle CreateFace(VertexHandle v1, VertexHandle v2, VertexHandle v3,FaceHandle f1, FaceHandle f2, FaceHandle f3);

		Vertex & GetVertex(VertexHandle vh) {return *vh;}
		Face & GetFace(FaceHandle fh) { return *fh;}

		list<Face>& GetFaces(){return faces;}
		list<Vertex>& GetVertices(){return vertexs;}

		AllVerticesIterator  AllVerticesBegin()
		{
			return vertexs.begin();	
		}
		AllVerticesIterator  AllVerticesEnd()
		{
			return vertexs.end();
		}
		AllFacesIterator  AllFacesBegin()
		{
			return faces.begin();
		}
		AllFacesIterator  AllFacesEnd()
		{
			return faces.end();
		}

		void Flip(FaceHandle f, int i);
		VertexHandle InsertInFace(FaceHandle f);
		VertexHandle InsertInEdge(FaceHandle f, int i);

		int NumberOfVertices()
		{
			return vertexs.size();
		}
		int NumberOfFaces()
		{
			return faces.size();
		}
		void SetInfiniteVertex(const VertexHandle& v) 
		{
			infinitevertex=v;
			hasinfinite = true;
		}
		bool IsInfinite(FaceHandle f) const ;
		bool IsInfinite(VertexHandle v) const 
		{
			return hasinfinite&&v == InfiniteVertex();
		}
		VertexHandle InfiniteVertex() const
		{
			return  infinitevertex;
		}
		bool HaveInfinite()
		{
			return hasinfinite;
		}
		void Clear();
		bool  IsEdge(VertexHandle va, VertexHandle vb, FaceHandle &fr,  int & i) const;
		bool IsEdge(VertexHandle va, VertexHandle vb) const;
		FaceCirculator  Incident_faces(VertexHandle v)
		{
			return GetVertex(v).IncidentFaces();
		}
		VertexCirculator  Incident_vertices(VertexHandle v)
		{
			return GetVertex(v).IncidentVertices();
		}
		void RemoveDegree_3(VertexHandle vh, FaceHandle fh);
		void DeleteFace(FaceHandle f);
		void DeleteVertex(VertexHandle v);
		VertexHandle JoinVertices(FaceHandle f, int i, VertexHandle v);
		void SetAdjacency(FaceHandle f0, int i0, FaceHandle f1, int i1) const;

		//test
		void OutputDXF(char *filename);
		void OutputSTL( char *filename );
	private:
		int MirrorIndex(FaceHandle f, int i) const;

	};


}

#endif
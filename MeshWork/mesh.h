#ifndef FILE_MESH
#define FILE_MESH
/**************************************************************************/
/* File:   mesh.hpp                                                */
/* Record mesh information                         */
/* Date:   26. jan. 2016                                                   */
/**************************************************************************/
namespace meshwork
{
	class Mesh
	{
	private:
	/**
       Representation of local mesh-size h
    */
    LocalH * lochfunc;
    ///全局最大尺寸
    double hglob;
    ///全局最小尺寸
    double hmin;

	vector<MeshPoint> nodes;
	//vector<double> nodesradius;
	//vector<int> nodesflag;
	vector<Element2d> surfelements;
	vector<FaceDescriptor> facedecoding;

	public:
		Mesh();
		void  SetGlobalH (double h);
		void  SetLocalH (const Point<3> & pmin, const Point<3> & pmax, double grading);
		void  RestrictLocalH (const Point<3> & p, double hloc);
		void  RestrictLocalHLine (const Point<3> & p1, const Point<3> & p2,double hloc);
		double GetGloabalH(){return hglob;}
		double  GetH (const Point<3> & p) const;
		LocalH & LocalHFunction () { return * lochfunc; }

	    int AddPoint(MeshPoint mp);
		int AddElement2d(Element2d ele);
		int AddFaceDescriptor(FaceDescriptor fd);

		//MeshPoint& GetPoint(int i){return nodes[i];}
		Element2d& GetElement2d(int i)
		{
			return surfelements[i];
		}
		FaceDescriptor& GetFaceDescriptor(int i) 
		{
			return facedecoding[i];
		}

		int GetNP () const 
		{ 
			return nodes.size(); 
		}
		MeshPoint & GetPoint(int i)
		{ 
			return nodes[i]; 
		}
		const MeshPoint & GetPoint(int i) const 
		{ 
			return nodes[i]; 
		}
		const vector<MeshPoint> & GetPoints() const 
		{ 
			return nodes;
		}
		vector<MeshPoint> & GetPoints() 
		{ 
			return nodes; 
		}

		int GetNSE () const 
		{ 
			return surfelements.size();
		}
		Element2d & GetSurfaceElement(int i)
		{ 
			return surfelements[i]; 
		}
		const Element2d & GetSurfaceElement(int i) const
		{ 
			return surfelements[i]; 
		}
		const vector<Element2d> & GetSurfaceElements() const
		{
			return surfelements;
		}
		vector<Element2d> & GetSurfaceElements() 
		{ 
			return surfelements; 
		}


	};

}
#endif
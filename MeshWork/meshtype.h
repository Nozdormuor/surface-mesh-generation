#ifndef FILE_MESHTYPE
#define FILE_MESHTYPE
/**************************************************************************/
/* File:   meshtype.hpp                                                */
/* fundamental mesh class                         */
/* Date:   19. feb. 2016                                                   */
/**************************************************************************/
namespace meshwork
{
	enum POINTTYPE { FIXEDPOINT = 1, EDGEPOINT = 2, SURFACEPOINT = 3, INNERPOINT = 4 };
	class MeshPoint : public Point<3>
	{
		
		POINTTYPE type;

	public:
		MeshPoint () 
		{ 
			;
		}

		MeshPoint (const Point<3> & ap, POINTTYPE apt = INNERPOINT): Point<3> (ap), type(apt) 
		{ 
			;
		}

		void SetPoint (const Point<3> & ap)
		{ 
			Point<3>::operator= (ap); 
		}

		POINTTYPE Type() const { return type; }

		void SetType(POINTTYPE at) { type = at; }

	};

	inline ostream & operator<<(ostream  & s, const MeshPoint & pt)
	{ 
		return (s << Point<3> (pt)); 
	}

	  /**
     Triangle element for surface mesh generation.
  */
  class Element2d
  { 
    /// point numbers
    int   pnum[3];

	int edgetype[3];//0:free,1:edge;
    /// surface nr
    int findex;

  public:
    ///
    Element2d ();
    ///
    Element2d (int pi1, int pi2, int pi3);
    ///
    Element2d (int pi1, int pi2, int pi3, int face);

    int & operator[] (int i) 
	{ 
		return pnum[i]; 
	}
    ///
    const int & operator[] (int i) const 
	{ 
		return pnum[i]; 
	}

    int & PNum (int i) 
	{ 
		return pnum[i]; 
	}
    ///
    const int & PNum (int i) const 
	{ 
		return pnum[i]; 
	}
    ///
    void SetIndex (int si) 
	{ 
		findex = si;
	}
    ///
    int GetIndex () const 
	{ 
		return findex; 
	}
	///
	void SetEdgeType(int index, int type) 
	{
		edgetype[index]=type;
	}
	///
	int GetEgdeType(int index) const 
	{
		return edgetype[index];
	}
    ///
    void GetBox (const vector<MeshPoint> & points, Box<3> & box) const;
    /// invert orientation
    inline void Invert ();
    /// first point number is smallest
    inline void NormalizeNumbering ();

	Element2d & operator= (const Element2d & p2)
	{
		for (int i = 0; i <3; i++) pnum[i] = p2.pnum[i]; 
		for (int i = 0; i <3; i++) edgetype[i] = p2.edgetype[i]; 
		return *this;
	}

    // friend ostream & operator<<(ostream  & s, const Element2d & el);
    friend class Mesh;
  };

  ostream & operator<<(ostream  & s, const Element2d & el);

  class FaceDescriptor
  {
	  /// which surface, 0 if not available
	  int surfnr;
	  /// domain nr inside
	  int domin;
	  /// domain nr outside
	  int domout;
	  /// root of linked list 
public:
	  vector<int> element;

  public:
	 FaceDescriptor();
	 FaceDescriptor(int surfnri);
	 FaceDescriptor(const FaceDescriptor& other);

	  int SurfNr () const 
	  { 
		  return surfnr;
	  }
	  int DomainIn () const 
	  { 
		  return domin; 
	  }
	  int DomainOut () const 
	  {
		  return domout;
	  }

	  void SetSurfNr (int sn) 
	  { 
		  surfnr = sn;
	  }
	  void SetDomainIn (int di) 
	  { 
		  domin = di;
	  }
	  void SetDomainOut (int dom) 
	  { 
		  domout = dom;
	  }

	  void AddElement(int i) 
	  {
		  element.push_back(i);
	  }
	  friend class Mesh;
  };

  ostream & operator<< (ostream  & s, const FaceDescriptor & fd);

}

#endif
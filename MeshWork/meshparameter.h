#ifndef FILE_MESHPARAMETER
#define FILE_MESHPARAMETER

/* *************************************************************************/
/* File:   meshparameter.hpp                                                 */
/* Function:  mesh parameter                        */
/* Date:   25. jan. 2016                                                   */
/* Class:  meshparameter                                     */
/* *************************************************************************/

namespace meshwork
{

  class  MeshingParameters
  {
  public:
    /**
       3d optimization strategy:
       // m .. move nodes
       // M .. move nodes, cheap functional
       // s .. swap faces
       // c .. combine elements
       // d .. divide elements
       // p .. plot, no pause
       // P .. plot, Pause
       // h .. Histogramm, no pause
       // H .. Histogramm, pause
       */
    const char * optimize3d;
    /// number of 3d optimization steps
    int optsteps3d;
    /**
       2d optimization strategy:
       // s .. swap, opt 6 lines/node
       // S .. swap, optimal elements
       // m .. move nodes
       // p .. plot, no pause
       // P .. plot, pause
       // c .. combine
       **/
    const char * optimize2d;
    /// number of 2d optimization steps
    int optsteps2d;
    /// power of error (to approximate max err optimization)
    double opterrpow;
    /// do block filling ?  
    int blockfill;
    /// block filling up to distance
    double filldist;
    /// radius of local environment (times h)
    double safety;
    /// radius of active environment (times h)
    double relinnersafety;
    /// use local h ?
    int uselocalh;
    /// grading for local h
    double grading;
    /// use delaunay meshing
    int delaunay;
    /// maximal mesh size
    double maxh;
    /// minimal mesh size
    double minh;
    /// file for meshsize
    const char * meshsizefilename;
    /// start surfacemeshing from everywhere in surface
    int startinsurface;
    /// check overlapping surfaces (debug)
    int checkoverlap;
    /// check overlapping surface mesh before volume meshing
    int checkoverlappingboundary;
    /// check chart boundary (sometimes too restrictive)
    int checkchartboundary;
    /// safty factor for curvatures (elemetns per radius)
    double curvaturesafety;
    /// minimal number of segments per edge
    double segmentsperedge;
    /// use parallel threads
    int parthread;
    /// weight of element size w.r.t element shape
    double elsizeweight;
    /// init with default values


    /// from mp3:
    /// give up quality class, 2d meshing
    int giveuptol2d;
    /// give up quality class, 3d meshing
    int giveuptol;
    /// maximal outer steps
    int maxoutersteps;
    /// class starting star-shape filling
    int starshapeclass;
    /// if non-zero, baseelement must have baseelnp points
    int baseelnp;        
    /// quality tolerances are handled less careful
    int sloppy;
  
    /// limit for max element angle (150-180)
    double badellimit;

    bool check_impossible;
  
    ///
    int secondorder;
    /// high order element curvature
    int elementorder;
    /// quad-dominated surface meshing
    int quad;
    ///
    int inverttets;
    ///
    int inverttrigs;
    ///
    int autozrefine;
    ///
    MeshingParameters ();
    ///
    void Print (ostream & ost) const;

    void CopyFrom(const MeshingParameters & other);
  };

  class STLParameters
  {
  public:
	  /// angle for edge detection
	  double yangle;
	  double contyangle; //edges continued with contyangle
	  /// angle of geometry edge at which the mesher should set a point
	  double edgecornerangle;
	  /// angle inside on chart
	  double chartangle;
	  /// angle for overlapping parts of char
	  double outerchartangle;
	  /// 0 .. no, 1 .. local, (2 .. global)
	  int usesearchtree;
	  ///
	  double resthatlasfac; 
	  int resthatlasenable;
	  double atlasminh;

	  double resthsurfcurvfac; 
	  int resthsurfcurvenable;

	  double resthchartdistfac;
	  int resthchartdistenable;

	  double resthcloseedgefac;
	  int resthcloseedgeenable;

	  double resthedgeanglefac;
	  int resthedgeangleenable;

	  double resthsurfmeshcurvfac;
	  int resthsurfmeshcurvenable;

	  double resthlinelengthfac;
	  int resthlinelengthenable;

	  ///
	  int recalc_h_opt;
	  ///
	  STLParameters();
	  ///
	  void Print (ostream & ost) const;
  };

}


#endif
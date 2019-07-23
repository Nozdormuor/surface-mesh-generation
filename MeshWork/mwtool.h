#ifndef FILE_MWTOOL
#define FILE_MWTOOL

/**************************************************************************/
/* File:   mwtool.hpp                                                */
/* principal function for MeshWork                         */
/* Date:   26. jan. 2016                                                   */
/**************************************************************************/

namespace meshwork
{
	double GetDistFromLine(const Point<3> & lp1, const Point<3> & lp2, Point<3> & p);
	double GetDistFromInfiniteLine(const Point<3>& lp1, const Point<3>& lp2, const Point<3>& p);
	double ComputeCylinderRadius (const Vec<3> & n1, const Vec<3> & n2,double h1, double h2);
	int AddPointIfNotExists(vector<Point<3>>& ap, const Point<3> p, double eps);
	int IntersectionOfCircleAndSegment3d( Point<3> p , double rc ,Point<3> plp1 , Point<3> plp2 , Point<3>  pip[2] );
	void GetNearNodes( int node ,  double maxr, double grad,vector<Point<3>>& snodesf ,vector<double> & nodesf_r ,vector<int> &nearnodes);
	int CalculateCircleFromTwoPoints3d( Point<3> p1 , Point<3> p2 , double r ,Point<3> &pcb , double *prb );
	int IsBoxAndTriangleIntersect( double x1 , double y1 ,double z1 , double x2 , double y2
		,double z2 , Point<3> p1 , Point<3> p2 , Point<3> p3 );
	int IsPointOnTriangle( Point<3> p ,  Point<3> p1 ,  Point<3> p2 , Point<3> p3 );
	int CalculateBoxOfCircle3d( Point<3> pcent , Point<3> p1 ,double r , double *x1 ,
		double *y1, double *z1 , double *x2 , double *y2 , double *z2 );
	void EquationOfPlane( Point<3> p1 , Point<3> p2 , double norm[4] );
	int IntersectionOfPlaneAndTriangle( Point<3> p1 , Point<3> p2 , Point<3> p3 , double a ,
		double b , double c , double d , Point<3>  pip[2] );
	int IsOnePointNearOthers( Point<3> p , double r, vector<Point<3>> &snodesf ,vector<double> &nodesf_r, vector<int> &nearnodes,int nj );
	 bool IsBoxContainPoint( Point<3> p,const Box<3> &bound,const Box<3> &rootbound);
	 bool IsTwoBoxOverlap(const Box<3> &a,const Box<3> &b);
	 bool IsTwoBoxNeighber(const Box<3> &a,const Box<3> &b);
	 bool IsPointOverlapWithBox(const Point<3> &p,const Box<3> &bd,const Box<3> &rootbound,double eps);
	 void BoxOfVertices( void **v , int num ,Box<3> &b ,void (*Funcpointofv)(Point<3> &p,void *v));
	 double SquareDistanceOfPointToBox(Point<3> p,const Box<3> &bd);
	 void GetTheLongestDistanceOfBox(const Box<3> &b,int &di, double *pdist=0);
	 void Copy3DPoint(const Point<3> &pfr,Point<3> &pto); //&?
	 double SquareDistanceOfInnerPointToBoxBound(Point<3> p,const Box<3> &bd);
	  bool IsTriangleAndBoxOverlap(Point<3> p1 ,Point<3> p2 ,Point<3> p3 ,const double bd[6],double eps );
	  bool IsTriangleAndRayIntersect(double p0[3],double p1[3],double p2[3],double orig[3],double dir[3]);//, double &t){
	 void jf_error(char *ch);
	 void vec_2p(double * , double * , double *) ;
	 int vec_unit(double *) ;
	 void vec_negation(double *vector);
	 double vec_value( double * ) ;
	 double vec_squarevalue( double *vec );
	 double vec_dotproduct(double * ,double *) ;
	 void vec_crossproduct(double * ,double * ,double * ) ;
	 double Distance3D( double * ,double *) ;
	 double SquareDistance3D( double * ,double * ) ;
	 double vec_blendproduct( double * , double * , double * ) ;
	 void normal_3p( double * , double * ,double * ,double * ) ;
	 double Volume_4p(Point<3> ,Point<3> ,Point<3> ,Point<3> ) ;

}


#endif
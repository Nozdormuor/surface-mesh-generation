#include"externalcommon.h"
#include"internalcommon.h"
namespace meshwork
{
	double GetDistFromLine(const Point<3> & lp1, const Point<3> & lp2, Point<3> & p)
	{
		Vec<3> vn = lp2 - lp1;
		Vec<3> v1 = p - lp1;
		Vec<3> v2 = lp2 - p;

		Point<3> pold = p;
		Point<3>  pp;
		if (v2 * vn <= 0) {pp = lp2; return (pold - pp).Length();}
		if (v1 * vn <= 0) {pp = lp1; return (pold - pp).Length();}

		double vnl = vn.Length();
		if (vnl == 0) {return Dist(lp1,p);}

		vn /= vnl;
		pp = lp1 + (v1 * vn) * vn;
		return (pold - pp).Length();
	}
	double GetDistFromInfiniteLine(const Point<3>& lp1, const Point<3>& lp2, const Point<3>& p)
	{
		Vec<3> vn= lp2-lp1;
		Vec<3> v1= p-lp1;

		double vnl = vn.Length();

		if (vnl == 0)
		{
			return Dist (lp1, p);
		}
		else
		{
			return Cross (vn, v1).Length() / vnl;
		}
	}
	double ComputeCylinderRadius (const Vec<3> & n1, const Vec<3> & n2,
		double h1, double h2)
	{
		Vec<3> t1, t2;
		double n11 = n1 * n1;
		double n12 = n1 * n2;
		double n22 = n2 * n2;
		double det = n11 * n22 - n12 * n12;

		if (fabs (det) < 1e-14 * n11 * n22)
			return 1e20;

		// a biorthogonal bases   (ti * nj) = delta_ij:
		t1 = (n22/det) * n1 + (-n12/det) * n2;
		t2 = (-n12/det) * n1 + (n11/det) * n2;

		// normalize:
		t1 /= t1.Length();
		t2 /= t2.Length();

		/*
		vector to center point has form
		v = lam1 n1 + lam2 n2
		and fulfills
		t2 v = h1/2
		t1 v = h2/2
		*/

		double lam1 = 0.5 * h2 / (n1 * t1);
		double lam2 = 0.5 * h1 / (n2 * t2);

		double rad = (lam1 * n1 + lam2 * n2).Length();

		return rad;
	}
	//add a point into a pointlist, return pointnumber
	int AddPointIfNotExists(vector<Point<3>>& ap, const Point<3> p, double eps)
	{
		double eps2 = eps*eps;
		for (int i = 0; i < ap.size(); i++)
			if (Dist(ap[i],p) <= eps2 ) 
				return i;
		ap.push_back(p);
		return ap.size()-1;
	}

	int IntersectionOfCircleAndSegment3d(  Point<3> p , double rc ,
		Point<3> plp1 , Point<3> plp2 , Point<3>  pip[2] )
	{
		double xc=p[0];
		double yc=p[1];
		double zc=p[2];

		double dp1p2 , p1m , dltsq ,dlt ,blt, v1 ;
		Point<3>  unitvecp1p2 , vecm ;
		int it = 0 ;

		dp1p2 = Dist( plp1 , plp2 ) ;  //if dp1p2 very small,error
		unitvecp1p2[0] = ( plp2[0] - plp1[0] )/dp1p2 ;
		unitvecp1p2[1] = ( plp2[1] - plp1[1] )/dp1p2 ;
		unitvecp1p2[2] = ( plp2[2] - plp1[2] )/dp1p2 ;
		p1m = ( xc - plp1[0] ) * unitvecp1p2[0]
		+ ( yc - plp1[1] ) * unitvecp1p2[1]
		+ ( zc - plp1[2] ) * unitvecp1p2[2] ;
		vecm[0] = plp1[0] + p1m * unitvecp1p2[0] ;
		vecm[1] = plp1[1] + p1m * unitvecp1p2[1] ;
		vecm[2] = plp1[2] + p1m * unitvecp1p2[2] ;
		dltsq = rc * rc-( (vecm[0]-xc)*( vecm[0]-xc) + (vecm[1]-yc)*( vecm[1]-yc)
			+(vecm[2]-zc) * ( vecm[2]-zc) ) ;
		if( dltsq < 0 ) return(0) ;
		dlt = sqrt( dltsq ) ;
		blt=0. ; //0.01 * min(dp1p2,rc) ;
		if( (v1 = p1m + dlt) >=-blt && v1<=dp1p2+blt ){
			pip[it][0] = plp1[0] + v1 * unitvecp1p2[0] ;
			pip[it][1] = plp1[1] + v1 * unitvecp1p2[1] ;
			pip[it++][2] = plp1[2] + v1 * unitvecp1p2[2] ;
		}
		if( (v1 = p1m - dlt) >=-blt && v1<=dp1p2+blt ){
			pip[it][0] = plp1[0] + v1 * unitvecp1p2[0] ;
			pip[it][1] = plp1[1] + v1 * unitvecp1p2[1] ;
			pip[it++][2] = plp1[2] + v1 * unitvecp1p2[2] ;
		}
		return(it) ;
	}

	void GetNearNodes( int node ,  double maxr, double grad,vector<Point<3>>& snodesf ,vector<double> & nodesf_r ,vector<int> &nearnodes)
	{
		int i ;
		double dist0  , dlx , drx , dly , dry , dlz , drz ;

		double a=nodesf_r[node];
		double g_maxincresize=grad+1.;
		dist0=a+2.01*a*g_maxincresize+a*g_maxincresize*g_maxincresize ;
		if( dist0>4.*maxr ) dist0=4.*maxr ;
		dlx=snodesf[node][0]-dist0 ;
		drx=snodesf[node][0]+dist0 ;
		dly=snodesf[node][1]-dist0 ;
		dry=snodesf[node][1]+dist0 ;
		dlz=snodesf[node][2]-dist0 ;
		drz=snodesf[node][2]+dist0 ;

		for( i=0 ; i<snodesf.size() ; i++ ){
			if( snodesf[i][0]<dlx||snodesf[i][0]>drx||
				snodesf[i][1]<dly||snodesf[i][1]>dry||
				snodesf[i][2]<dlz||snodesf[i][2]>drz )continue ;
			if( 0.9*Dist( snodesf[i],snodesf[node] ) < nodesf_r[i]+
				nodesf_r[node]+2.*sqrt(nodesf_r[i]*nodesf_r[node])&&(i!= node) ){
					nearnodes.push_back(i);
					if( nearnodes.size() >= 500 )                  /* test every time or here only 1 */
						cout<<( " error-near_snodef\n") ;
			}
		}
	}
	int CalculateCircleFromTwoPoints3d( Point<3> p1 , Point<3> p2 , double r ,Point<3> &pcb , double *prb )
	{
		double d ;
		pcb[0] = ( p1[0] + p2[0] ) /2. ;
		pcb[1] = ( p1[1] + p2[1] ) /2. ;
		pcb[2] = ( p1[2] + p2[2] ) /2. ;
		d = Dist( p1 , p2 ) /2. ;
		if( r <= 1.0000001*d ) return(0) ; //rev.12.11
		*prb = sqrt( r * r - d * d ) ;
		return(1) ;
	}
	int IsBoxAndTriangleIntersect( double x1 , double y1 ,double z1 , double x2 , double y2
		,double z2 , Point<3> p1 , Point<3> p2 , Point<3> p3 )
	{
		if( p1[0]<x1 && p2[0]<x1 && p3[0]<x1 || p1[0]>x2 && p2[0]>x2 && p3[0]>x2
			|| p1[1]<y1 && p2[1]<y1 && p3[1]<y1 || p1[1]>y2 && p2[1]>y2 && p3[1]>y2
			|| p1[2]<z1 && p2[2]<z1 && p3[2]<z1 || p1[2]>z2 && p2[2]>z2 && p3[2]>z2
			) return(0) ;
		return(1) ;
	}
	int IsPointOnTriangle( Point<3> p ,  Point<3> p1 ,  Point<3> p2 , Point<3> p3 )
	{
		Vec<3> vpp1,vpp2,vpp3,vp1p2,vp1p3,
			varea1,varea2,varea3,vare ;
		//vec_2p( p , p1 , vpp1 ) ;
		vpp1=p1-p;
		//vec_2p( p , p2 , vpp2 ) ;
		vpp2=p2-p;
		//vec_2p( p , p3 , vpp3 ) ;
		vpp3=p3-p;
		//vec_2p( p1 , p2 , vp1p2 ) ;
		vp1p2=p2-p1;
		//vec_2p( p1 , p3 , vp1p3 ) ;
		vp1p3=p3-p1;
		//vec_crop( vpp1 , vpp2 , varea1 ) ;
		varea1=Cross(vpp1 , vpp2 );
		//vec_crop( vpp2 , vpp3 , varea2 ) ;
		varea2=Cross( vpp2 , vpp3);
		//vec_crop( vpp3 , vpp1 , varea3 ) ;
		varea3=Cross(vpp3 , vpp1);
		//vec_crop( vp1p2 , vp1p3 , vare ) ;
		vare=Cross(vp1p2 , vp1p3);
#define eps4 1.01  // 1.1->1.01
		if( varea1.Length() +varea2.Length() +varea3.Length()<
			eps4 *   vare.Length() ) return(1) ;
		return(0) ;
	}
	int CalculateBoxOfCircle3d( Point<3> pcent , Point<3> p1 ,double r , double *x1 ,
		double *y1, double *z1 , double *x2 , double *y2 , double *z2 )
	{
		double mx , my , mz , dx , dy , dz , ds ;
		mx = p1[0] - pcent[0] ;
		my = p1[1] - pcent[1] ;
		mz = p1[2] - pcent[2] ;
		ds = sqrt( mx * mx + my * my + mz * mz ) ;
		if( ds < 0.0000001 ) return(0) ;               //rev.12.11
		dx = r * sqrt( my * my +mz * mz ) /ds * 1.01 ;
		dy = r * sqrt( mx * mx +mz * mz ) /ds * 1.01 ;
		dz = r * sqrt( my * my +mx * mx ) /ds * 1.01 ;
		*x1 = pcent[0] - dx ; *y1 = pcent[1] - dy ; *z1 = pcent[2] - dz ;
		*x2 = pcent[0] + dx ; *y2 = pcent[1] + dy ; *z2 = pcent[2] + dz ;
		return(1) ; /* change to box of sphere good or not */
	}
	void EquationOfPlane( Point<3> p1 , Point<3> p2 , double norm[4] )
	{
		norm[0] = p2[0] - p1[0] ;
		norm[1] = p2[1] - p1[1] ;
		norm[2] = p2[2] - p1[2] ;
		norm[3] = -p1[0] * norm[0] - p1[1] * norm[1] - p1[2] *norm[2] ;
	}
	int IntersectionOfPlaneAndTriangle( Point<3> p1 , Point<3> p2 , Point<3> p3 , double a ,
		double b , double c , double d , Point<3>  pip[2] )
	{
		double h1 , h2 , h3 ;
		h1 = a * p1[0] + b * p1[1] + c * p1[2] + d ;
		h2 = a * p2[0] + b * p2[1] + c * p2[2] + d ;
		h3 = a * p3[0] + b * p3[1] + c * p3[2] + d ;
		if( h1 >=0&&h2 >= 0&&h3 >= 0 || h1 <= 0&&h2 <= 0&&h3 <= 0 ) return(0) ;
		if( h1 * h2 >= 0&&fabs(h1-h3)>=0.0000001&&fabs(h2-h3)>=0.0000001 ){
			pip[0][0] = ( h1 * p3[0] - h3 * p1[0] ) / ( h1 - h3 ) ;
			pip[0][1] = ( h1 * p3[1] - h3 * p1[1] ) / ( h1 - h3 ) ;
			pip[0][2] = ( h1 * p3[2] - h3 * p1[2] ) / ( h1 - h3 ) ;
			pip[1][0] = ( h2 * p3[0] - h3 * p2[0] ) / ( h2 - h3 ) ;
			pip[1][1] = ( h2 * p3[1] - h3 * p2[1] ) / ( h2 - h3 ) ;
			pip[1][2] = ( h2 * p3[2] - h3 * p2[2] ) / ( h2 - h3 ) ;
			return 1 ;
		}else if( h3 * h2 >= 0&&fabs(h1-h3)>=0.0000001&&fabs(h2-h1)>=0.0000001 ){
			pip[0][0] = ( h1 * p3[0] - h3 * p1[0] ) / ( h1 - h3 ) ;
			pip[0][1] = ( h1 * p3[1] - h3 * p1[1] ) / ( h1 - h3 ) ;
			pip[0][2] = ( h1 * p3[2] - h3 * p1[2] ) / ( h1 - h3 ) ;
			pip[1][0] = ( h2 * p1[0] - h1 * p2[0] ) / ( h2 - h1 ) ;
			pip[1][1] = ( h2 * p1[1] - h1 * p2[1] ) / ( h2 - h1 ) ;
			pip[1][2] = ( h2 * p1[2] - h1 * p2[2] ) / ( h2 - h1 ) ;
			return 1 ;
		}else if( h1 * h3 >= 0&&fabs(h1-h2)>=0.0000001&&fabs(h2-h3)>=0.0000001 ){
			pip[0][0] = ( h1 * p2[0] - h2 * p1[0] ) / ( h1 - h2 ) ;
			pip[0][1] = ( h1 * p2[1] - h2 * p1[1] ) / ( h1 - h2 ) ;
			pip[0][2] = ( h1 * p2[2] - h2 * p1[2] ) / ( h1 - h2 ) ;
			pip[1][0] = ( h2 * p3[0] - h3 * p2[0] ) / ( h2 - h3 ) ;
			pip[1][1] = ( h2 * p3[1] - h3 * p2[1] ) / ( h2 - h3 ) ;
			pip[1][2] = ( h2 * p3[2] - h3 * p2[2] ) / ( h2 - h3 ) ;
			return 1 ;
		}
		return(0) ;   // rev.12.11
	}
	int IsOnePointNearOthers( Point<3> p , double r, vector<Point<3>> &snodesf ,vector<double> &nodesf_r, vector<int> &nearnodes,int nj )
	{
		int i ;    double wp ;
		for(i=0 ;i<nearnodes.size(); i++ ){
			if(nearnodes[i]==nj)continue;
			wp = Dist( p , snodesf[ nearnodes[i] ]) ;
			if( 1.3 * wp < nodesf_r[ nearnodes[i] ] + r) return(1) ; //1.6=>1.4 2002 9 11
		}
		return(0) ;
	}

	double SquareDistanceOfPointToBox( Point<3> p,const Box<3> &bd)
	{

		double a[3];
		double q=0;
		for(int i=0; i<3; i++){
			if(p[i]>bd[i+3]) a[i]=p[i]-bd[i+3];
			else if(p[i]<bd[i]) a[i]=bd[i]-p[i];
			else a[i]=0;
			q+=a[i]*a[i];
		}
		return q;
	}
	double SquareDistanceOfInnerPointToBoxBound( Point<3> p,const Box<3> &bd)
	{

		double a=min(p[0]-bd[0],bd[3]-p[0]);
		double b=min(p[1]-bd[1],bd[4]-p[1]);
		double c=min(p[2]-bd[2],bd[5]-p[2]);
		double d= min(c,min(a,b));
		return d*d;
	}
	void GetTheLongestDistanceOfBox(const Box<3> &b,int &di, double *pdist)
	{

		di=0;
		double dist=0.;
		for(int i=0; i<3; i++)
			if(b[i+3]-b[i]>dist){
				dist=b[i+3]-b[i];
				di=i;
			}
			if(pdist!=0) *pdist=dist;
	}
	bool IsTwoBoxOverlap(const Box<3> &a,const Box<3> &b)
	{

		if(a[0]>b[3]||a[1]>b[4]||a[2]>b[5]||a[3]<b[0]||a[4]<b[1]||a[5]<b[2])
			return false;
		return true;
	}


	void Copy3DPoint(const Point<3> &pfr,Point<3> &pto)
	{ //const or not?

		pto[0]=pfr[0];
		pto[1]=pfr[1];
		pto[2]=pfr[2];
	}

	void jf_error(char *ch)
	{

		printf("%s\n",ch);
		exit(1);
	}


	bool IsBoxContainPoint( Point<3> p,const Box<3> &bound,const Box<3> &rootbound)
	{

		if(p[0]<bound[0]||p[1]<bound[1]||p[2]<bound[2]||p[0]>bound[3]||p[1]>bound[4]||p[2]>bound[5])
			return false;
		else if(bound[0]!=rootbound[0]&&p[0]==bound[0]||
			bound[1]!=rootbound[1]&&p[1]==bound[1]||
			bound[2]!=rootbound[2]&&p[2]==bound[2] ) //need or not to follow the convention?
			return false;
		else 
			return true;
	}

	bool IsPointOverlapWithBox(const Point<3> &p,const Box<3> &bd,const Box<3> &rootbound,double eps)
	{

		Box<3> bound;

		double a[3];
		for(int i=0; i<3; i++)
			a[i]=bd[i+3]-bd[i];
		for(int i=0; i<3; i++){
			bound[i]=bd[i]-eps*a[i];
			bound[i+3]=bd[i+3]+eps*a[i];
		}//convention
		if(p[0]<bound[0]||p[1]<bound[1]||p[2]<bound[2]||p[0]>bound[3]||p[1]>bound[4]||p[2]>bound[5])
			return false;
		else if(bound[0]!=rootbound[0]&&p[0]==bound[0]||
			bound[1]!=rootbound[1]&&p[1]==bound[1]||
			bound[2]!=rootbound[2]&&p[2]==bound[2] )
			return false;
		else 
			return true;
	}

	bool IsTwoBoxNeighber(const Box<3> &a,const Box<3> &b)
	{

		if(a[0]>b[3]||a[1]>b[4]||a[2]>b[5]||a[3]<b[0]||a[4]<b[1]||a[5]<b[2])
			return false;
		return true;
	}

	void BoxOfVertices( void **v, int num ,Box<3> &box,void (*pofv)(Point<3> &p,void *v) )
	{

		int i , j ;
		double a ;

		Point<3> p;

		pofv(p,v[0]);
		for(i=0 ; i<3 ; i++ )
		{
			box[i]=box[i+3]=p[i] ;
		}
		for(j=1 ; j<num ; j++ )
		{
			pofv(p,v[j]);
			for( i=0 ; i<3 ; i++ )
			{
				if( p[i]<box[i] ) box[i]=p[i] ;
				if( p[i]>box[i+3] ) box[i+3]=p[i] ;
			}
		}
		a=max( box[3]-box[0] ,max(box[4]-box[1],box[5]-box[2]) ) ;
		for( i=0 ; i<3 ; i++ ){
			box[i] -= 0.01*a ;
			box[i+3] += 0.01*a ; //keep unchanged for a undegenerate 3D box or use a,b and c?
		}
	}

	/*
	void vec_2p(double *pointa , double *pointb , double *vector)
	{
	int i;
	for(i=0; i<3; i++)
	vector[i]=pointb[i]-pointa[i];
	}
	void
	vec_neg(double *vector)
	{
	vector[0]= -vector[0];
	vector[1]= -vector[1];
	vector[2]= -vector[2];
	}
	int
	vec_uni(double *vector)
	{
	double len;  int i;
	len=(double)sqrt(vector[0]*vector[0]+vector[1]*vector[1]
	+vector[2]*vector[2]);
	if(len<=0.00000000001) return(0);        // eps_2 zero vector bound :min_siz/2. 
	for(i=0; i<3; i++)   vector[i]/=len;
	return(1);
	}
	double vec_val( double *vec )
	{
	return( sqrt(vec[0] * vec[0] +vec[1] * vec[1] + vec[2] *vec[2] )) ;
	}
	double vec_sqval( double *vec )
	{
	return vec[0] * vec[0] +vec[1] * vec[1] + vec[2] *vec[2]  ;
	}
	double
	vec_dotp(double *vector1,double *vector2)
	{
	return(vector1[0]*vector2[0]+vector1[1]*vector2[1]+vector1[2]*vector2[2]);
	}
	void
	vec_crop(double *vector1,double *vector2,double *vector3)
	{
	vector3[0]=vector1[1]*vector2[2]-vector2[1]*vector1[2];
	vector3[1]=vector2[0]*vector1[2]-vector1[0]*vector2[2];
	vector3[2]=vector1[0]*vector2[1]-vector2[0]*vector1[1];
	}

	double
	Distance3D( double *v1,double *v2){
	return( sqrt( (v1[0]-v2[0])*(v1[0]-v2[0]) +
	(v1[1]-v2[1])*(v1[1]-v2[1]) +
	(v1[2]-v2[2])*(v1[2]-v2[2]) )
	) ;
	}
	double
	SqDistance3D( double *v1,double *v2){
	return( (v1[0]-v2[0])*(v1[0]-v2[0]) +
	(v1[1]-v2[1])*(v1[1]-v2[1]) +
	(v1[2]-v2[2])*(v1[2]-v2[2])
	) ;
	}


	double vec_blep( double *vec1 , double *vec2 , double *vec3 )
	{
	double vecp[3] ;
	vec_crop( vec2 , vec3 , vecp ) ;
	return( vec_dotp(vec1 , vecp ) ) ;
	}
	void norm_3p( double *p1 , double *p2 ,double *p3 ,double *normal )
	{
	double v12[3] , v13[3] ;
	vec_2p( p1 , p2 , v12 ) ;
	vec_2p( p1 , p3 , v13 ) ; // the direction of crop has changed ?!
	vec_crop( v12 , v13 , normal ) ;
	}

	double VolumOf4p(jf_point p0,jf_point p1,jf_point p2,jf_point p3)
	{
	jf_point vec ,p01,p02,p03 ;
	p01[0] = p1[0] -p0[0] ;   p01[1] = p1[1] -p0[1] ;   p01[2] = p1[2] -p0[2] ;
	p02[0] = p2[0] -p0[0] ;   p02[1] = p2[1] -p0[1] ;   p02[2] = p2[2] -p0[2] ;
	p03[0] = p3[0] -p0[0] ;   p03[1] = p3[1] -p0[1] ;   p03[2] = p3[2] -p0[2] ;
	vec_crop(p01,p02,vec) ;                            // valu  or parem ? 
	return vec_dotp(vec,p03) ;
	}

	int triBoxOverlap(double boxcenter[3],double boxhalfsize[3],double triverts[3][3]);

	bool isTriangleBoxOver(jf_point p1 ,jf_point p2 ,jf_point p3 ,const double bd[6],double eps ){

	//  int i ;

	double bound[6];
	double a[3];
	for(int i=0; i<3; i++)
	a[i]=bd[i+3]-bd[i];
	for(int i=0; i<3; i++){
	bound[i]=bd[i]-eps*a[i];
	bound[i+3]=bd[i+3]+eps*a[i];
	}
	double boxcenter[3],boxhalfsize[3],triverts[3][3];
	for(int i=0; i<3; i++){
	boxcenter[i]=(bound[i]+bound[i+3])/2.;
	boxhalfsize[i]=(bound[i+3]-bound[i])/2;
	triverts[0][i]=p1[i];
	triverts[1][i]=p2[i];
	triverts[2][i]=p3[i];
	}
	if(triBoxOverlap(boxcenter,boxhalfsize,triverts)==0) return false;
	else return true;
	//  for(i=0 ; i<3 ; i++ )
	//  if( (p1[i]<box[i]&&p2[i]<box[i]&&p3[i]<box[i])||
	//   (p1[i]>box[i+3]&&p2[i]>box[i+3]&&p3[i]>box[i+3]) )return false ;
	// return true ;
	}*/
	bool IsTriangleAndRayIntersect(double p0[3],double p1[3],double p2[3],double orig[3],double dir[3])
	{//, double &t){
		double tvec[3], pvec[3], qvec[3]; //revised from gel
		double det,inv_det;
		double edge[3][3];
		//	return true;
		vec_2p(p0,p1,edge[0]);
		vec_2p(p0,p2,edge[2]);
		/* begin calculating determinant - also used to calculate U parameter */
		vec_crossproduct(dir,edge[2],pvec);
		/* if determinant is near zero, ray lies in plane of triangle */
		det = vec_dotproduct(edge[0], pvec);

		if (det > -0.00000000001 && det < 0.00000000001)
			return false;
		inv_det = 1.0 / det;

		/* calculate distance from v0 to ray origin */
		vec_2p(p0,orig,tvec);

		/* calculate U parameter and test bounds */
		double u = vec_dotproduct(tvec, pvec) * inv_det;
		if (u < 0.0 || u > 1.0)
			return false;

		/* prepare to test V parameter */
		vec_crossproduct(tvec, edge[0],qvec);

		/* calculate V parameter and test bounds */
		double v = vec_dotproduct(dir, qvec) * inv_det;
		if (v < 0.0 || u + v > 1.0)
			return false;

		/* calculate t, ray intersects triangle */
		return  vec_dotproduct(edge[2], qvec) * inv_det>=0;

		//   return true;
	}


	/********************************************************/

	/* AABB-triangle overlap test code                      */

	/* by Tomas Akenine-Moller                              */

	/* Function: int triBoxOverlap(double boxcenter[3],      */

	/*          double boxhalfsize[3],double triverts[3][3]); */

	/* History:                                             */

	/*   2001-03-05: released the code in its first version */

	/*   2001-06-18: changed the order of the tests, faster */

	/*                                                      */

	/* Acknowledgement: Many thanks to Pierre Terdiman for  */

	/* suggestions and discussions on how to optimize code. */

	/* Thanks to David Hunt for finding a ">="-bug!         */

	/********************************************************/

#include <math.h>

#include <stdio.h>



#define X 0

#define Y 1

#define Z 2



#define CROSS(dest,v1,v2) \
	dest[0]=v1[1]*v2[2]-v1[2]*v2[1]; \
	dest[1]=v1[2]*v2[0]-v1[0]*v2[2]; \
	dest[2]=v1[0]*v2[1]-v1[1]*v2[0]; 



#define DOT(v1,v2) (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])



#define SUB(dest,v1,v2) \
	dest[0]=v1[0]-v2[0]; \
	dest[1]=v1[1]-v2[1]; \
	dest[2]=v1[2]-v2[2]; 



#define FINDMINMAX(x0,x1,x2,min,max) \
	min = max = x0;   \
	if(x1<min) min=x1;\
	if(x1>max) max=x1;\
	if(x2<min) min=x2;\
	if(x2>max) max=x2;



	int planeBoxOverlap(double normal[3], double vert[3], double maxbox[3])	// -NJMP-

	{

		int q;

		double vmin[3],vmax[3],v;

		for(q=X;q<=Z;q++)

		{

			v=vert[q];					// -NJMP-

			if(normal[q]>0.0f)

			{

				vmin[q]=-maxbox[q] - v;	// -NJMP-

				vmax[q]= maxbox[q] - v;	// -NJMP-

			}

			else

			{

				vmin[q]= maxbox[q] - v;	// -NJMP-

				vmax[q]=-maxbox[q] - v;	// -NJMP-

			}

		}

		if(DOT(normal,vmin)>0.0f) return 0;	// -NJMP-

		if(DOT(normal,vmax)>=0.0f) return 1;	// -NJMP-



		return 0;

	}





	/*======================== X-tests ========================*/

#define AXISTEST_X01(a, b, fa, fb)			   \
	p0 = a*v0[Y] - b*v0[Z];			       	   \
	p2 = a*v2[Y] - b*v2[Z];			       	   \
	if(p0<p2) {min=p0; max=p2;} else {min=p2; max=p0;} \
	rad = fa * boxhalfsize[Y] + fb * boxhalfsize[Z];   \
	if(min>rad || max<-rad) return 0;



#define AXISTEST_X2(a, b, fa, fb)			   \
	p0 = a*v0[Y] - b*v0[Z];			           \
	p1 = a*v1[Y] - b*v1[Z];			       	   \
	if(p0<p1) {min=p0; max=p1;} else {min=p1; max=p0;} \
	rad = fa * boxhalfsize[Y] + fb * boxhalfsize[Z];   \
	if(min>rad || max<-rad) return 0;



	/*======================== Y-tests ========================*/

#define AXISTEST_Y02(a, b, fa, fb)			   \
	p0 = -a*v0[X] + b*v0[Z];		      	   \
	p2 = -a*v2[X] + b*v2[Z];	       	       	   \
	if(p0<p2) {min=p0; max=p2;} else {min=p2; max=p0;} \
	rad = fa * boxhalfsize[X] + fb * boxhalfsize[Z];   \
	if(min>rad || max<-rad) return 0;



#define AXISTEST_Y1(a, b, fa, fb)			   \
	p0 = -a*v0[X] + b*v0[Z];		      	   \
	p1 = -a*v1[X] + b*v1[Z];	     	       	   \
	if(p0<p1) {min=p0; max=p1;} else {min=p1; max=p0;} \
	rad = fa * boxhalfsize[X] + fb * boxhalfsize[Z];   \
	if(min>rad || max<-rad) return 0;



	/*======================== Z-tests ========================*/



#define AXISTEST_Z12(a, b, fa, fb)			   \
	p1 = a*v1[X] - b*v1[Y];			           \
	p2 = a*v2[X] - b*v2[Y];			       	   \
	if(p2<p1) {min=p2; max=p1;} else {min=p1; max=p2;} \
	rad = fa * boxhalfsize[X] + fb * boxhalfsize[Y];   \
	if(min>rad || max<-rad) return 0;



#define AXISTEST_Z0(a, b, fa, fb)			   \
	p0 = a*v0[X] - b*v0[Y];				   \
	p1 = a*v1[X] - b*v1[Y];			           \
	if(p0<p1) {min=p0; max=p1;} else {min=p1; max=p0;} \
	rad = fa * boxhalfsize[X] + fb * boxhalfsize[Y];   \
	if(min>rad || max<-rad) return 0;



	int triBoxOverlap(double boxcenter[3],double boxhalfsize[3],double triverts[3][3])

	{



		/*    use separating axis theorem to test overlap between triangle and box */

		/*    need to test for overlap in these directions: */

		/*    1) the {x,y,z}-directions (actually, since we use the AABB of the triangle */

		/*       we do not even need to test these) */

		/*    2) normal of the triangle */

		/*    3) crossproduct(edge from tri, {x,y,z}-directin) */

		/*       this gives 3x3=9 more tests */

		double v0[3],v1[3],v2[3];

		//   double axis[3];

		double min,max,p0,p1,p2,rad,fex,fey,fez;		// -NJMP- "d" local variable removed

		double normal[3],e0[3],e1[3],e2[3];



		/* This is the fastest branch on Sun */

		/* move everything so that the boxcenter is in (0,0,0) */

		SUB(v0,triverts[0],boxcenter);

		SUB(v1,triverts[1],boxcenter);

		SUB(v2,triverts[2],boxcenter);



		/* compute triangle edges */

		SUB(e0,v1,v0);      /* tri edge 0 */

		SUB(e1,v2,v1);      /* tri edge 1 */

		SUB(e2,v0,v2);      /* tri edge 2 */



		/* Bullet 3:  */

		/*  test the 9 tests first (this was faster) */

		fex = fabs(e0[X]);

		fey = fabs(e0[Y]);

		fez = fabs(e0[Z]);

		AXISTEST_X01(e0[Z], e0[Y], fez, fey);

		AXISTEST_Y02(e0[Z], e0[X], fez, fex);

		AXISTEST_Z12(e0[Y], e0[X], fey, fex);



		fex = fabs(e1[X]);

		fey = fabs(e1[Y]);

		fez = fabs(e1[Z]);

		AXISTEST_X01(e1[Z], e1[Y], fez, fey);

		AXISTEST_Y02(e1[Z], e1[X], fez, fex);

		AXISTEST_Z0(e1[Y], e1[X], fey, fex);



		fex = fabs(e2[X]);

		fey = fabs(e2[Y]);

		fez = fabs(e2[Z]);

		AXISTEST_X2(e2[Z], e2[Y], fez, fey);

		AXISTEST_Y1(e2[Z], e2[X], fez, fex);

		AXISTEST_Z12(e2[Y], e2[X], fey, fex);



		/* Bullet 1: */

		/*  first test overlap in the {x,y,z}-directions */

		/*  find min, max of the triangle each direction, and test for overlap in */

		/*  that direction -- this is equivalent to testing a minimal AABB around */

		/*  the triangle against the AABB */



		/* test in X-direction */

		FINDMINMAX(v0[X],v1[X],v2[X],min,max);

		if(min>boxhalfsize[X] || max<-boxhalfsize[X]) return 0;



		/* test in Y-direction */

		FINDMINMAX(v0[Y],v1[Y],v2[Y],min,max);

		if(min>boxhalfsize[Y] || max<-boxhalfsize[Y]) return 0;



		/* test in Z-direction */

		FINDMINMAX(v0[Z],v1[Z],v2[Z],min,max);

		if(min>boxhalfsize[Z] || max<-boxhalfsize[Z]) return 0;



		/* Bullet 2: */

		/*  test if the box intersects the plane of the triangle */

		/*  compute plane equation of triangle: normal*x+d=0 */

		CROSS(normal,e0,e1);

		// -NJMP- (line removed here)

		if(!planeBoxOverlap(normal,v0,boxhalfsize)) return 0;	// -NJMP-



		return 1;   /* box and triangle overlaps */

	}

	//vector operator
	void vec_2p(double *pointa , double *pointb , double *vector)
	{
		int i;
		for(i=0; i<3; i++)
			vector[i]=pointb[i]-pointa[i];
	}

	int  vec_unit(double *vector)
	{
		double len;  int i;
		len=(double)sqrt(vector[0]*vector[0]+vector[1]*vector[1]
		+vector[2]*vector[2]);
		if(len<=0.0000000000000000001) return(0);        /* eps_2 zero vector bound :min_siz/2. */
		for(i=0; i<3; i++)   vector[i]/=len;
		return(1);
	}

	void  vec_negation(double *vector)
	{
		vector[0]= -vector[0];
		vector[1]= -vector[1];
		vector[2]= -vector[2];
	}

	double vec_value( double *vec )
	{
		return( sqrt(vec[0] * vec[0] +vec[1] * vec[1] + vec[2] *vec[2] )) ;
	}
	double
		vec_dotproduct(double *vector1,double *vector2)
	{
		return(vector1[0]*vector2[0]+vector1[1]*vector2[1]+vector1[2]*vector2[2]);
	}
	void
		vec_crossproduct(double *vector1,double *vector2,double *vector3)
	{
		vector3[0]=vector1[1]*vector2[2]-vector2[1]*vector1[2];
		vector3[1]=vector2[0]*vector1[2]-vector1[0]*vector2[2];
		vector3[2]=vector1[0]*vector2[1]-vector2[0]*vector1[1];
	}

	double vec_squarevalue( double *vec )
	{
		return vec[0] * vec[0] +vec[1] * vec[1] + vec[2] *vec[2]  ;
	}

	double
		Distance3D( double *v1,double *v2){
			return( sqrt( (v1[0]-v2[0])*(v1[0]-v2[0]) +
				(v1[1]-v2[1])*(v1[1]-v2[1]) +
				(v1[2]-v2[2])*(v1[2]-v2[2]) )
				) ;
	}

	double
		SquareDistance3D( double *v1,double *v2){
			return ( (v1[0]-v2[0])*(v1[0]-v2[0]) +
				(v1[1]-v2[1])*(v1[1]-v2[1]) +
				(v1[2]-v2[2])*(v1[2]-v2[2])
				) ;
	}

	double vec_blendproduct( double *vec1 , double *vec2 , double *vec3 )
	{
		double vecp[3] ;
		vec_crossproduct( vec2 , vec3 , vecp ) ;
		return( vec_dotproduct(vec1 , vecp ) ) ;
	}

	void normal_3p( double *p1 , double *p2 ,double *p3 ,double *normal )
	{
		double v12[3] , v13[3] ;
		vec_2p( p1 , p2 , v12 ) ;
		vec_2p( p1 , p3 , v13 ) ; // the direction of crop has changed ?!
		vec_crossproduct( v12 , v13 , normal ) ;
	}

	double Volume_4p(Point<3> p0,Point<3> p1,Point<3> p2,Point<3> p3)
	{
		double vec[3] ,p01[3],p02[3],p03[3] ;
		p01[0] = p1[0] -p0[0] ;   p01[1] = p1[1] -p0[1] ;   p01[2] = p1[2] -p0[2] ;
		p02[0] = p2[0] -p0[0] ;   p02[1] = p2[1] -p0[1] ;   p02[2] = p2[2] -p0[2] ;
		p03[0] = p3[0] -p0[0] ;   p03[1] = p3[1] -p0[1] ;   p03[2] = p3[2] -p0[2] ;
		vec_crossproduct(p01,p02,vec) ;                            /* valu  or parem ? */
		return vec_dotproduct(vec,p03) ;
	}

	double modevects( double p0[3] , double p1[3] )
	{
		double vec[3] ;
		vec_crossproduct( p0 , p1  , vec ) ;
		return( sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]) ) ;
	}

	double area4p(double p0[3] , double p1[3] ,double p2[3] ,double p3[3] )
	{
		double p01[3],p02[3],p03[3], v1[3] , v2[3] ;
		p01[0] = p1[0] -p0[0] ;   p01[1] = p1[1] -p0[1] ;   p01[2] = p1[2] -p0[2] ;
		p02[0] = p2[0] -p0[0] ;   p02[1] = p2[1] -p0[1] ;   p02[2] = p2[2] -p0[2] ;
		p03[0] = p3[0] -p0[0] ;   p03[1] = p3[1] -p0[1] ;   p03[2] = p3[2] -p0[2] ;
		v1[0] = p1[0] -p3[0] ;   v1[1] = p1[1] -p3[1] ;   v1[2] = p1[2] -p3[2] ;
		v2[0] = p2[0] -p3[0] ;   v2[1] = p2[1] -p3[1] ;   v2[2] = p2[2] -p3[2] ;
		return( modevects(p01,p02)+modevects(p02,p03)+modevects(p01,p03)+
			modevects(v1,v2)) ;
	}


	bool IsTriangleAndBoxOverlap(Point<3> p1 ,Point<3> p2 ,Point<3> p3 ,const double bd[6],double eps ){

		//  int i ;

		double bound[6];
		double a[3];
		for(int i=0; i<3; i++)
			a[i]=bd[i+3]-bd[i];
		for(int i=0; i<3; i++)
		{
			bound[i]=bd[i]-eps*a[i];
			bound[i+3]=bd[i+3]+eps*a[i];
		}
		double boxcenter[3],boxhalfsize[3],triverts[3][3];
		for(int i=0; i<3; i++)
		{
			boxcenter[i]=(bound[i]+bound[i+3])/2.;
			boxhalfsize[i]=(bound[i+3]-bound[i])/2;
			triverts[0][i]=p1[i];
			triverts[1][i]=p2[i];
			triverts[2][i]=p3[i];
		}
		if(triBoxOverlap(boxcenter,boxhalfsize,triverts)==0) return false;
		else return true;
		//  for(i=0 ; i<3 ; i++ )
		//  if( (p1[i]<box[i]&&p2[i]<box[i]&&p3[i]<box[i])||
		//   (p1[i]>box[i+3]&&p2[i]>box[i+3]&&p3[i]>box[i+3]) )return false ;
		// return true ;
	}


	void CompRatioPoint( double p0[3] ,double p1[3] ,double dlt,double p[3] )
	{

		int i ;
		for( i=0 ; i<3 ; i++ )
			p[i]=p0[i]+dlt*(p1[i]-p0[i]) ;

	}
	double sqDistPointToRatioPoint(double p[3],double p0[3],double p1[3],double dlt)
	{
		double prt[3];
		CompRatioPoint(p0,p1,dlt,prt);
		return SquareDistance3D(p,prt);
	}

	double sqDistPointToTri(double p[3],double p0[3],double p1[3],double p2[3])
	{//,double sqdist0){
		//return signed_distance(p,p0,p1,p2);

		double v0p[3],v20[3],v01[3];
		vec_2p(p0,p,v0p);
		vec_2p(p2,p0,v20);
		vec_2p(p0,p1,v01);
		//double nm012[3];
		//vec_crop(v20,v01,nm012);
		//double sqdptoinner=vec_dotp(nm012,v0p);
		//sqdptoinner*=sqdptoinner/vec_sqval(nm012);
		//if(sqdptoinner>=sqdist0) return sqdptoinner;

		double d0p20=vec_dotproduct(v0p,v20);
		double d0p01=vec_dotproduct(v0p,v01);
		if(d0p20>=0&&d0p01<=0) return SquareDistance3D(p,p0);

		double v1p[3],v12[3];
		vec_2p(p1,p,v1p);
		vec_2p(p1,p2,v12);
		double d1p01=vec_dotproduct(v1p,v01);
		double d1p12=vec_dotproduct(v1p,v12);
		if(d1p01>=0&&d1p12<=0) return SquareDistance3D(p,p1);

		double v2p[3];
		vec_2p(p2,p,v2p);
		double d2p12=vec_dotproduct(v2p,v12);
		double d2p20=vec_dotproduct(v2p,v20);
		if(d2p12>=0&&d2p20<=0) return SquareDistance3D(p,p2);

		double 
			nm012[3],
			nm01p[3],nm12p[3],nm20p[3];
		vec_crossproduct(v20,v01,nm012);

		vec_crossproduct(v01,v0p,nm01p);
		double dt01=vec_dotproduct(nm012,nm01p);
		if(dt01<=0&&d0p01>=0&&d1p01<=0)
			return sqDistPointToRatioPoint(p,p0,p1,d0p01/(d0p01-d1p01)); 
		//		return sqDistPointToSeg3D(p,p0,p1); // rt==0;

		vec_crossproduct(v12,v1p,nm12p);
		double dt12=vec_dotproduct(nm012,nm12p);
		if(dt12<=0&&d1p12>=0&&d2p12<=0)
			return sqDistPointToRatioPoint(p,p1,p2,d1p12/(d1p12-d2p12)); 
		//		return sqDistPointToSeg3D(p,p1,p2);

		vec_crossproduct(v20,v2p,nm20p);
		double dt20=vec_dotproduct(nm012,nm20p);
		if(dt20<=0&&d2p20>=0&&d0p20<=0)
			return sqDistPointToRatioPoint(p,p2,p0,d2p20/(d2p20-d0p20)); 
		//		return sqDistPointToSeg3D(p,p2,p0);

		if(dt01>=0&&dt12>=0&&dt20>=0){
			double a=vec_dotproduct(nm012,v0p);
			return a*a/vec_squarevalue(nm012);
			//		return sqdptoinner;
		}else{
			double a=vec_dotproduct(nm012,v0p);
			return a*a/vec_squarevalue(nm012); //2014-5-27 previous version 2011-4-13
		}
	}

	int is_point_same( Point<3>p1 , Point<3> p2 )
	{

		if( p1[0]==p2[0]&&p1[1]==p2[1]&&p1[2]==p2[2] ) return(1) ;
		return(0) ;
	}

	double DistPointToSegm( double p[3] , double ps[3] , double pe[3] ,double pnear[3],int *ip)
	{

		double vsp[3] , vse[3] /*, vsep*/,pj[3] ;
		double dsp , dse , dep ,dist , dsptose ;

		vec_2p(ps , p , vsp ) ;
		vec_2p(ps , pe , vse ) ;
		dsp=Distance3D(ps , p) ;
		dse=Distance3D(ps , pe) ;
		if(vec_unit(vse)==0)
		{
			if(pnear!=NULL) memcpy(pnear,ps,3*sizeof(double));
			if(ip!=NULL) *ip=0;
			return Distance3D(p,ps);
		}
		if( (dsptose=vec_dotproduct(vsp , vse))>0.000001*dse&&dsptose<0.999999*dse )
		{
			for(int i=0 ; i<3 ; i++ )    pj[i]=ps[i]+dsptose*vse[i] ;
			if(pnear!=NULL) memcpy(pnear,pj,3*sizeof(double));
			if(ip!=NULL) *ip=1;
			return Distance3D(p,pj);

		}else
		{
			dep=Distance3D(pe,p) ;
			dist=min( dsp , dep ) ;
			if(dist==dsp){	 
				if(pnear!=NULL) memcpy(pnear,ps,3*sizeof(double));
				if(ip!=NULL) *ip=0;
			}else{	 
				if(pnear!=NULL) memcpy(pnear,pe,3*sizeof(double));
				if(ip!=NULL) *ip=2;
			}
		}
		return dist ;
	}

	void PointProjectTo2PLine(double p[3] ,double p0[3] ,double p1[3] ,double pb[3])
	{

		double v01[3] , v0p[3] ;
		double dlt ;
		int i ;

		vec_2p(p0,p1,v01) ;
		vec_2p(p0,p,v0p) ;
		if(vec_unit( v01 )==0) jf_error("norm err ppt2pl" ) ;
		dlt= vec_dotproduct(v01,v0p) ;
		for( i=0 ; i<3 ; i++ )
			pb[i]=p0[i]+dlt*v01[i] ;
	}

	double Comp2fAngleTestDegenerate( double pa[3] ,double pb[3] ,double pc[3] ,double pd[3] )
	{

		double norm1[3] , norm2[3] ;

		normal_3p( pa , pb , pc , norm1 ) ;
		normal_3p( pa , pb , pd , norm2 ) ;
		if(vec_unit( norm1 )==0)
			return -2.;
		if(vec_unit( norm2 )==0)
			return -2.;
		return( vec_dotproduct(norm1 , norm2 ) ) ;
	}

	double compMinCos2Faceangle1pTo4p(double p[3],double p4t[4][3])
	{

		double a1,a2,a3,a4;

		if((a1=Comp2fAngleTestDegenerate(p4t[0],p4t[1],p4t[2],p))<=0) return -1;
		if((a2=Comp2fAngleTestDegenerate(p4t[1],p4t[2],p4t[0],p))<=0) return -1;
		if((a3=Comp2fAngleTestDegenerate(p4t[2],p4t[3],p4t[0],p))<=0) return -1;
		if((a4=Comp2fAngleTestDegenerate(p4t[3],p4t[0],p4t[2],p))<=0) return -1;
		return min(min(a1,a2),min(a3,a4));
	}

	void rotatePoints(double (*pts)[3],int numpt)
	{

		double tmpt[3];
		memcpy(tmpt,pts[0],3*sizeof(double));
		for(int i =0; i<numpt-1; i++)
			memcpy(pts[i],pts[i+1],3*sizeof(double));
		memcpy(pts[numpt-1],tmpt,3*sizeof(double));
	}


}

#include "externalcommon.h"
#include "internalcommon.h"
namespace meshwork
{
	/* ******************************* STLLine ******************************* */
	STLLine :: STLLine()
	{
		split = 0;
	};

	int STLLine :: GetNS() const
	{
		if (pts.size() <= 1) {return 0;}
		return pts.size()-1;
	}

	void STLLine :: GetSeg(int nr, int& p1, int& p2) const
	{
		p1 = pts[nr];
		p2 = pts[nr+1];
	}

	double STLLine :: GetSegLen(const vector<Point<3> >& ap, int nr) const
	{
		return Dist(ap[PNum(nr)],ap[PNum(nr+1)]);
	}

	double STLLine :: GetLength(const vector<Point<3> >& ap) const
	{
		double len = 0;
		for (int i = 1; i < pts.size(); i++)
		{
			len += (ap[pts[(i)]] - ap[pts[i-1]]).Length();
		}
		return len;
	}

	Point<3> STLLine :: GetPointInDist(const vector<Point<3> >& ap, double dist, int& index) const
	{
		if (dist <= 0)
		{
			index = 0;
			return ap[StartP()];
		}

		double len = 0;
		int i;
		for (i = 0; i < pts.size()-1; i++)
		{
			double seglen = Dist (ap[pts[i]],
				ap[pts[i+1]]);

			if (len + seglen > dist)
			{
				index = i;
				double relval = (dist - len) / (seglen + 1e-16);
				Vec<3> v = ap[pts[i+1]]-ap[pts[i]];
				return ap[pts[i]] + relval * v;
			}

			len += seglen;
		}

		index = pts.size() - 2;
		return ap[EndP()];
	}

	void STLLine :: GetBoundingBox (const vector<Point<3> > & ap, Box<3> & box) const
	{
		if(pts.size()==0)
			return;
		box.Set (ap[pts[0]]);
		for (int i = 1; i < pts.size(); i++)
			box.Add (ap[pts[i]]);
	}

	STLLine STLLine :: Mesh(const vector<Point<3> >& ap, vector<Point<3>>& mp, double ghi,class Mesh& mesh) const
	{

		STLLine line;

		//stlgh = ghi; //uebergangsloesung!!!!

		double len = GetLength(ap);
		if(len==0.)
			return line;
		double inthl = 0; //integral of 1/h
		double dist = 0;
		double h;
		int ind;
		Point<3> p;

		Box<3> bbox;
		GetBoundingBox (ap, bbox);
		double diam = bbox.Diam();

		double minh = mesh.LocalHFunction().GetMinH (bbox.PMin(), bbox.PMax());

		double maxseglen = 0;
		for (int i = 0; i < GetNS(); i++)
			maxseglen = max (maxseglen, GetSegLen (ap, i));

		int nph = 10+int(maxseglen / minh); 

		vector<double> inthi(GetNS()*nph);
		vector<double> curvelen(GetNS()*nph);


		for (int i =0; i < GetNS(); i++)
		{
			//double seglen = GetSegLen(ap,i);
			for (int j = 0; j < nph; j++)
			{
				p = GetPointInDist(ap,dist,ind);
				//h = GetH(p,dist/len);
				h = mesh.GetH(p);


				dist += GetSegLen(ap,i)/(double)nph;

				inthl += GetSegLen(ap,i)/nph/(h);
				inthi[i*nph+j] = GetSegLen(ap,i)/nph/h;
				curvelen[i*nph+j] = GetSegLen(ap,i)/nph;
			}
		}


		int inthlint = int(inthl+1);

		if ( (inthlint < 3) && (StartP() == EndP()))
		{
			inthlint = 3;
		}
		if ( (inthlint == 1) && ShouldSplit())
		{
			inthlint = 2; 
		}

		double fact = inthl/(double)inthlint;
		dist = 0;
		int j = 0;


		p = ap[StartP()];
		int pn = AddPointIfNotExists(mp, p, 1e-10*diam);

		line.AddPoint(pn);

		inthl = 0; //restart each meshseg
		for (int i = 1; i <= inthlint; i++)
		{
			while (inthl < 1.000000001 && j < inthi.size())
			{
				inthl += inthi[j]/fact;
				dist += curvelen[j];
				j++;
			}

			//went too far:
			j--;
			double tofar = (inthl - 1)/inthi[j];
			inthl -= tofar*inthi[j];
			dist -= tofar*curvelen[j]*fact;

			if (i == inthlint && fabs(dist - len) >= 1E-8) 
			{
				cout<<("meshline failed!!!"); 
			}

			if (i != inthlint) 
			{
				p = GetPointInDist(ap,dist,ind);
				pn = AddPointIfNotExists(mp, p, 1e-10*diam);
				line.AddPoint(pn);

			}

			inthl = tofar*inthi[j];
			dist += tofar*curvelen[j]*fact;
			j++;
		}
		p = ap[EndP()];
		pn = AddPointIfNotExists(mp, p, 1e-10*diam);
		line.AddPoint(pn);


		return line;
	}
	/* ******************************* STLFace ******************************* */

	STLFace :: STLFace()
	{
		;
	}
	//set node on face
	int STLFace ::SetNodesOnSurface( Point<3> pa , Point<3> pb ,double ra, double rb,const vector<Point<3> >& ap,const vector<STLTriangle> &tris, Point<3>  p[2]  )
	{
		double r , rcir ;
		Point<3>  pcent  ;
		//if(sortpointrv(snodesf[ia],snodesf[ja])>0) swap(&ia,&ja);
		//r = ( nodesf_r[ia] + nodesf_r[ja] ) ;
		r=ra+rb;
		if( CalculateCircleFromTwoPoints3d( pa , pb , r , pcent , &rcir ) == 0 ) return(-1) ;
		//cout<<rcir<<endl;
		if( rcir < 0.001 * r ){
			if( IsPointOnSurface(ap,tris,pcent)==0 )
				return(-1) ;
			p[0] = pcent  ;
			return(0) ;
		}
		return( IntersectionOfCircleAndSurface( ap ,tris,pcent ,pa , rcir , p  ) ) ;
	}

	int STLFace ::IntersectionOfCircleAndSurface( const vector<Point<3> >& pois,const vector<STLTriangle> &trips ,
		Point<3> p1 , Point<3> p2 , double r , Point<3> p[2] ){

			int i ,it =-1 ,ret ;
			double x1 , y1 ,z1 ,x2 , y2 , z2 ,v1 ,norm_p[4] ;
			Point<3> pip[2] , pjp[2] ;

			if(CalculateBoxOfCircle3d( p1 , p2 , r , &x1 , &y1 , &z1 , &x2 ,&y2 , &z2 )==0){
				CalculateBoxOfCircle3d( p1 , p2 , r , &x1 , &y1 , &z1 , &x2 ,&y2 , &z2 );
				//     drawline( p1,p2,1, 3  ) ;
				p2[0]=p1[0]+0.2; p2[1]=p1[1]+0.2;	 p2[2]=p1[2]+0.2;
				//     drawline( p1,p2,1, 3  ) ;
				cout<<( " maybe local size too small ,clear and redefines it " ) ;
			}
			EquationOfPlane( p1 , p2 , norm_p ) ;
			for( i = 0 ; i< GetFaceTriangleNum() ; i++ ){
				if( IsBoxAndTriangleIntersect( x1,y1,z1,x2,y2,z2,pois[trips[GetTriangle(i)][0]] ,
					pois[trips[GetTriangle(i)][1]] ,pois[trips[GetTriangle(i)][2]]) == 0 ) continue ;
				if( IntersectionOfPlaneAndTriangle( pois[trips[GetTriangle(i)][0]] , pois[trips[GetTriangle(i)][1]]
				,pois[trips[GetTriangle(i)][2]],norm_p[0],norm_p[1],norm_p[2],norm_p[3],pip)
					== 0 ) continue ;
				if( Dist( pip[0] , pip[1] ) < 0.0005 * r ){
					if( (v1 = Dist( pip[0] , p1 ) )> 0.9995 * r&& v1 < 1.0005 * r ){
						/* or v2 is in this scope */
						p[++it][0] = ( pip[0][0] + pip[1][0] ) / 2. ;
						p[it][1] = ( pip[0][1] + pip[1][1] ) / 2. ;
						p[it][2] = ( pip[0][2] + pip[1][2] ) / 2. ;
						//trinum[it] = i ;
						if( it == 1 ){
							if( Dist( p[0] , p[1] ) > 0.05 * r ) return(1) ;
							it = 0 ;
						}
					}
				}else{
					if( (ret = IntersectionOfCircleAndSegment3d(p1 ,r ,
						pip[0] , pip[1] , pjp ) ) == 0 ) continue ;
					if( ret == 2 ){
						//        if( it == 0 ) prompt( "error : set_node_zigf ",1 ) ;
						p[0] = pjp[0];
						p[1] = pjp[1];
						return(1) ;
					}
					++it ;
					p[it] =pjp[0] ;
					if( it == 1 ){
						if( Dist( p[0] , p[1] ) > 0.05 * r ) return(1) ;
						it = 0 ;
					}
				}
			}
			return( it ) ;
	}

	int STLFace ::IsPointOnSurface( const vector<Point<3>> &pois , const vector<STLTriangle> &trips ,Point<3> p){

		int i ;

		for( i = 0 ; i < GetFaceTriangleNum() ; i++ ){
#define Tmin_size 0.1
			if( IsBoxAndTriangleIntersect(p[0]-Tmin_size,p[1]-Tmin_size,p[2]-Tmin_size,
				p[0]+Tmin_size,p[1]+Tmin_size,p[2]+Tmin_size,pois[trips[GetTriangle(i)][0]] ,
				pois[trips[GetTriangle(i)][1]] ,pois[trips[GetTriangle(i)][2]])== 0 ) continue ;
			if( IsPointOnTriangle( p , pois[trips[GetTriangle(i)][0]] ,
				pois[trips[GetTriangle(i)][1]] , pois[trips[GetTriangle(i)][2]] ) ) {
					return(1) ;
			}
		}
		return(0) ; /* maybe err when x=c,y=c or z=c for some triangle */
	}

	vector<int> STLFace ::Mesh(const vector<Point<3> >& ap,const vector<STLTriangle> &tris,vector<Point<3>>& mp,const vector<STLLine> &mline, double ghi,class Mesh& mesh) 
	{
		if(facetriangles.size()==0)
			return vector<int>();
		vector<Point<3>> snodesf;
		vector<double> nodesf_r;
		vector<int> facenodes;
		int temp;
		int ind;
		//find edges node
		for(int i = 0; i< GetFaceLineNum();i++)
		{
			Box<3> bbox;
			if(mline[GetLine(i)].NP()==0)
				continue;
			mline[GetLine(i)].GetBoundingBox (mp, bbox);
			double diam = bbox.Diam();
			temp=snodesf.size();
			ind=AddPointIfNotExists(snodesf,mp[mline[GetLine(i)].StartP()],1e-10*diam);
			if (ind==temp) 
				facenodes.push_back(mline[GetLine(i)].StartP());
			for(int j = 1; j < mline[GetLine(i)].NP()-1;j++)
			{
				snodesf.push_back(mp[mline[GetLine(i)].PNum(j)]);
				facenodes.push_back(mline[GetLine(i)].PNum(j));
			}
			temp=snodesf.size();
			ind=AddPointIfNotExists(snodesf,mp[mline[GetLine(i)].EndP()],1e-10*diam);
			if (ind==temp) 
				facenodes.push_back(mline[GetLine(i)].EndP());
		}
		//if no edge, create 2 nodes
		if(snodesf.size()==0)
		{
			snodesf.push_back(ap[tris[GetTriangle(0)].PNum(0)]);
			mp.push_back(ap[tris[GetTriangle(0)].PNum(0)]);
			facenodes.push_back(mp.size()-1);
			double h= mesh.GetH(ap[tris[GetTriangle(0)].PNum(0)]);
			Point<3> pjp[2] ;
			for(int i = 0; i < GetFaceTriangleNum(); i++)
			{
				if( Dist(ap[tris[GetTriangle(i)].PNum(0) ],ap[ tris[GetTriangle(i)].PNum(1) ])>
					0.00001&&IntersectionOfCircleAndSegment3d( snodesf[0],2.*h,
					ap[tris[GetTriangle(i)].PNum(0) ],ap[ tris[GetTriangle(i)].PNum(1)] , pjp)){
						snodesf.push_back( pjp[0]);
						mp.push_back(pjp[0]);
						facenodes.push_back(mp.size()-1);
						break;

				}
				if( Dist(ap[tris[GetTriangle(i)].PNum(1) ],ap[ tris[GetTriangle(i)].PNum(2) ])>
					0.00001&&IntersectionOfCircleAndSegment3d( snodesf[0],2.*h,
					ap[tris[GetTriangle(i)].PNum(1) ],ap[ tris[GetTriangle(i)].PNum(2)] , pjp))
				{
					snodesf.push_back( pjp[0]);
					mp.push_back(pjp[0]);
					facenodes.push_back(mp.size()-1);
					break;

				}
				if( Dist(ap[tris[GetTriangle(i)].PNum(2) ],ap[ tris[GetTriangle(i)].PNum(0) ])>
					0.00001&&IntersectionOfCircleAndSegment3d( snodesf[0],2.*h,
					ap[tris[GetTriangle(i)].PNum(2) ],ap[ tris[GetTriangle(i)].PNum(0)] , pjp)){
						snodesf.push_back( pjp[0]);
						mp.push_back(pjp[0]);
						facenodes.push_back(mp.size()-1);
						break;

				}
			}
			if(snodesf.size()==1)
				cout<<"err"<<endl;
		}
		//calculate node radius of edge
		for(int i = 0; i < snodesf.size(); i++)
		{
			nodesf_r.push_back(mesh.GetH(snodesf[i])/2.);
		}
		//cout<<snodesf.size()<<endl;
		// fill node in STLface
		{
			double r ; Point<3>  p[2] ;
			int i,j,nj,itw,jtw ;
			vector<int> nearnodes;
			if(snodesf.size()<=0) 
				cout<<("err zidfill");

			//for(int i = 0;i< snodesf.size() ; i++)
			//{
			//	for(int j = 0; j <snodesf.size() ; j++)
			//	{
			//		if(i!=j&&Dist(snodesf[i],snodesf[j])<1)
			//		{
			//			cout<<Dist(snodesf[i],snodesf[j])<<endl;
			//		}
			//	}
			//}
			for( i=0 ; i< snodesf.size() ; i++ )
			{
				nearnodes.clear();
				GetNearNodes( i,ghi,mparam.grading,snodesf,nodesf_r,nearnodes ) ;
				for (j=0 ; j<nearnodes.size() ; j++ ){
					if( (nj=nearnodes[j])<i )continue ;
					if( (jtw=SetNodesOnSurface( snodesf[i], snodesf[nj],nodesf_r[i],nodesf_r[nj], ap, tris, p)) < 0 ) continue ;
					for ( itw=0 ; itw<= jtw ; itw++ ){
						r=mesh.GetH( p[itw])/2. ;
						if( IsOnePointNearOthers( p[itw],r ,snodesf,nodesf_r,nearnodes,nj ) )  continue ;
						//addnodef( jf , p[itw] , trinum[itw] , r ) ;
						snodesf.push_back(p[itw]);
						nodesf_r.push_back(r);
						nearnodes.push_back(snodesf.size()-1)  ;
						mp.push_back(p[itw]);
						facenodes.push_back(mp.size()-1);
						if(nearnodes.size()>=500) cout<< "error nearnode\n " ;
					}
				}
			}
		}
		//cout<<facenodes.size()<<endl;
		return facenodes;
	}



}

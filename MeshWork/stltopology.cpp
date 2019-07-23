#include"externalcommon.h"
#include"internalcommon.h"


namespace meshwork
{

	//+++++++++++++++++++++++++++++STLReadTriangle+++++++++++++++++++++++++++++++++++++

	STLReadTriangle :: STLReadTriangle (const Point<3> * apts,
		const Vec<3> & anormal)
	{
		pts[0] = apts[0];
		pts[1] = apts[1];
		pts[2] = apts[2]; 
		normal = anormal;
	}


	//+++++++++++++++++++++++++++++STLTriangle+++++++++++++++++++++++++++++++++++++
	STLTriangle :: STLTriangle(const int * apts)
	{
		pts[0] = apts[0];
		pts[1] = apts[1];
		pts[2] = apts[2];

		facenum = 0;
	}

	void STLTriangle :: SetNormal (const Vec<3> & n)
	{
		double len = n.Length();
		if (len > 0)
		{
			normal = n;
			normal.Normalize();
		}
		else
		{
			normal = Vec<3> (1, 0, 0);
		}
	}

	int STLTriangle :: AnotherPointNum(int pointindex1, int pointindex2)
	{
		for (int i=0; i<3; i++)
		{
			if(i!=pointindex1&&i!=pointindex2)
				return pts[i];
		}
		return -1;
	}

	int STLTriangle :: AnotherPointIndex(int pointnum1, int pointnum2)
	{
		for (int i=0; i<3; i++)
		{
			if(pts[i]!=pointnum1&&pts[i]!=pointnum2)
				return i;
		}
		return -1;
	}

	int STLTriangle ::IsCoEdge(STLTriangle trig)
	{
		if(pts[0]==trig[0]&&pts[1]==trig[1]||pts[0]==trig[0]&&pts[1]==trig[2]||pts[0]==trig[0]&&pts[2]==trig[1]||pts[0]==trig[0]&&pts[2]==trig[2]) return 1;
		if(pts[0]==trig[1]&&pts[1]==trig[0]||pts[0]==trig[1]&&pts[1]==trig[2]||pts[0]==trig[1]&&pts[2]==trig[0]||pts[0]==trig[1]&&pts[2]==trig[2]) return 1;
		if(pts[0]==trig[2]&&pts[1]==trig[0]||pts[0]==trig[2]&&pts[1]==trig[1]||pts[0]==trig[2]&&pts[2]==trig[0]||pts[0]==trig[2]&&pts[2]==trig[1]) return 1;
		return 0;
	}
	int STLTriangle :: IsNeighbourFrom(const STLTriangle& t) const
	{
		//triangles must have same orientation!!!

		for(int i = 0; i <= 2; i++)
			for(int j = 0; j <= 2; j++)
				if (t.pts[(i+1)%3] == pts[j] && 
					t.pts[i] == pts[(j+1)%3])

					return 1;

		return 0;      
	}

	int STLTriangle :: IsWrongNeighbourFrom(const STLTriangle& t) const
	{
		//triangles have not same orientation!!!
		for(int i = 0; i <= 2; i++)
			for(int j = 0; j <= 2; j++)
				if (t.pts[(i+1)%3] == pts[(j+1)%3] &&
					t.pts[i] == pts[j])

					return 1;

		return 0;      
	}

	void STLTriangle :: GetNeighbourPoints(const STLTriangle& t, int& p1, int& p2) const
	{
		vector<int> p;
		if(PNum(0)==t.PNum(0)||PNum(0)==t.PNum(1)||PNum(0)==t.PNum(2))
			p.push_back(PNum(0));
		if(PNum(1)==t.PNum(0)||PNum(1)==t.PNum(1)||PNum(1)==t.PNum(2))
			p.push_back(PNum(1));
		if(PNum(2)==t.PNum(0)||PNum(2)==t.PNum(1)||PNum(2)==t.PNum(2))
			p.push_back(PNum(2));
		p1=p.front(); p2 = p.back();
		return;
	}
	int STLTriangle :: GetNeighbourPointsAndOpposite(const STLTriangle& t, int& p1, int& p2, int& po) const
	{
		for(int i = 0; i <3; i++)
			for(int j = 0; j < 3; j++)
				if (t.PNum((i+1) %3)== PNum(j%3) &&
					t.PNum(i%3) == PNum((j+1)%3))
				{
					p1 = PNum(j%3); 
					p2 = PNum((j+1)%3); 
					po = PNum((j+2)%3); 
					return 1;
				}

				return 0;
	}
	ostream& operator<<(ostream& os, const STLTriangle& t)
	{
		os << "points:" <<"[";
		os << t[0] << ",";
		os << t[1] << ",";
		os << t[2] << "]" << endl;

		os << "topedge:" <<"[";
		os << t.EdgeNum(0) << ",";
		os << t.EdgeNum(1) << ",";
		os << t.EdgeNum(2) << "]" << endl;

		os << "neighbour triangles:" <<"[";
		os << t.NBTrigNum(0) << ",";
		os << t.NBTrigNum(1) << ",";
		os << t.NBTrigNum(2) << "]" << endl;

		return os;
	}

	//+++++++++++++++++++++++++++++STLTopEdge+++++++++++++++++++++++++++++++++++++
	STLTopEdge :: STLTopEdge ()
	{
		pts[0] = pts[1] = 0;
		trigs[0] = trigs[1] = 0;
		cosangle = 1;
		status = ED_UNDEFINED;
	}

	STLTopEdge :: STLTopEdge (int p1, int p2, int trig1, int trig2)
	{ 
		pts[0] = p1; 
		pts[1] = p2; 
		trigs[0] = trig1; 
		trigs[1] = trig2; 
		cosangle = 1;
		status = ED_UNDEFINED;
	}

	ostream& operator<<(ostream& os, const STLTopEdge& e)
	{
		os <<  "points: " <<"[";
		os << e[0] << ",";
		os << e[1] << "]";
		os <<  "trigs: " <<"[";
		os << e.TrigNum(0)  << ",";
		os << e.TrigNum(1) << "]";
		return os;
	}

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//++++++++++++++ ++++++++++++++++++STLTopology++++++++++++++++++++++++
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	double  STLTopology::geom_tol_fact = 1E-6;
	STLTopology :: STLTopology()
	{
		;
	}

	STLTopology :: ~STLTopology()
	{
		;
	}



	void STLTopology :: InitSTLGeometry(const vector<STLReadTriangle> & readtrias)
	{

		cout << "number of triangles = " << readtrias.size() << endl;
		if (!readtrias.size()) return;

		boundingbox.Set (readtrias[0][0]);
		for (int i = 0; i < readtrias.size(); i++)
		{
			for (int k = 0; k < 3; k++)
				boundingbox.Add (readtrias[i][k]);
		}
		//		cout<<boundingbox<<endl;
		Box<3> bb = boundingbox;
		bb.Increase (1);

		pointtree = new Point3dTree (bb.PMin(), bb.PMax());
		vector<int> pintersect;
		pointtol = boundingbox.Diam() * geom_tol_fact;

		cout << "point tolerance = "<< pointtol<<endl;

		for(int i = 0; i < readtrias.size(); i++)
		{
			const STLReadTriangle & t = readtrias[i];

			STLTriangle st;
			st.SetNormal (t.Normal());

			for (int k = 0; k < 3; k++)
			{
				Point<3> p = t[k];
				Point<3> pmin = p - Vec<3> (pointtol, pointtol, pointtol);
				Point<3> pmax = p + Vec<3> (pointtol, pointtol, pointtol);

				pointtree->GetIntersecting (pmin, pmax, pintersect);

				if (pintersect.size() > 1)
					cout<<"too many close points" << endl;
				int foundpos = -1;
				if (pintersect.size())
					foundpos = pintersect[0];
				if (foundpos == -1)
				{
					foundpos = AddPoint(p);
					pointtree->Insert (p, foundpos);
				}
				if (Dist(p, points[foundpos]) > 1e-10)
					cout << "identify close points: " << p << " " << points[foundpos] 
				<< ", dist = " << Dist(p, points[foundpos])
					<< endl;
				st[k] = foundpos;
			}

			if ( (st[0] == st[1]) ||
				(st[0] == st[2]) || 
				(st[1] == st[2]) )
			{
				//PrintError("STL Triangle degenerated");
				cout << "STL Triangle degenerated" << endl;
			}
			else
			{
				AddTriangle(st);
			}

		} 
		//FindNeighbourTrigs();
		FindTrigsOfPoints();
		BuildSTLTopEdge();
		FindEdgesOfPoint();
		FindNeighbourTrigs();
		FindEdgesOfTrig();
		CalculateCosangleOfEdge();
		CalculateVectorOfPoint();
	}

	int STLTopology :: GetPointNum (const Point<3> & p)
	{
		Point<3> pmin = p - Vec<3> (pointtol, pointtol, pointtol);
		Point<3> pmax = p + Vec<3> (pointtol, pointtol, pointtol);

		vector<int> pintersect;

		pointtree->GetIntersecting (pmin, pmax, pintersect);
		if (pintersect.size() == 1)
			return pintersect[0];
		else 
			return -1;
	}

	void STLTopology :: AddTriangle(const STLTriangle& t)
	{
		trias.push_back(t);

		const Point<3> & p1 = GetPoint (t.PNum(0));
		const Point<3> & p2 = GetPoint (t.PNum(1));
		const Point<3> & p3 = GetPoint (t.PNum(2));

		Box<3> box;
		box.Set (p1);
		box.Add (p2);
		box.Add (p3);
		/*
		//  Point<3> pmin(p1), pmax(p1);
		pmin.SetToMin (p2);
		pmin.SetToMin (p3);
		pmax.SetToMax (p2);
		pmax.SetToMax (p3);
		*/

		trias.back().box = box; 
		trias.back().center = Center (p1, p2, p3);
		double r1 = Dist (p1, trias.back().center);
		double r2 = Dist (p2, trias.back().center);
		double r3 = Dist (p3, trias.back().center);
		trias.back().rad = max (max (r1, r2), r3);

	}

	int STLTopology :: GetTopEdgeNum (int pi1, int pi2) const
	{
		if (! topedgesperpoint.Size()) return -1;
		for(int i=0; i<topedgesperpoint.RowSize(pi1);i++)
		{
			if(topedges[topedgesperpoint.Get(pi1,i)][0]==pi1&&topedges[topedgesperpoint.Get(pi1,i)][1]==pi2||
				topedges[topedgesperpoint.Get(pi1,i)][0]==pi2&&topedges[topedgesperpoint.Get(pi1,i)][1]==pi1)
				return topedgesperpoint.Get(pi1,i);
		}
		return -1;
	}

	void STLTopology :: FindTrigsOfPoints()
	{
		trigsperpoint.SetSize(points.size());
		for ( int i =0; i < trias.size(); i++)
		{
			for( int j = 0; j<3; j++)
			{
				trigsperpoint.Add(trias[i][j],i);
			}
		}
	}

	void STLTopology ::BuildSTLTopEdge()
	{
		Table<int> pointsaround;
		for(int i = 0; i<points.size(); i++)
		{
			pointsaround.SetSize(0);
			for(int j=0; j <trigsperpoint.RowSize(i); j++)
			{

				int pointIndex=(trias[trigsperpoint.Get(i,j)]).PIndex(i);
				int cw=(trias[trigsperpoint.Get(i,j)]).GetCw(pointIndex);
				int ccw=(trias[trigsperpoint.Get(i,j)]).GetCcw(pointIndex);
				if(cw>i)
				{
					int k;
					for( k =0; k< pointsaround.Size(); k++)
					{
						if(pointsaround.Get(k,0)==cw)
						{
							pointsaround.Add(k,trigsperpoint.Get(i,j));
							break;
						}
					}
					if(k==pointsaround.Size())
					{
						pointsaround.Add(k,cw);
						pointsaround.Add(k,trigsperpoint.Get(i,j));
					}
				}
				if(ccw>i)
				{
					int k;
					for( k =0; k< pointsaround.Size(); k++)
					{
						if(pointsaround.Get(k,0)==ccw)
						{
							pointsaround.Add(k,trigsperpoint.Get(i,j));
							break;
						}
					}
					if(k==pointsaround.Size())
					{
						pointsaround.Add(k,ccw);
						pointsaround.Add(k,trigsperpoint.Get(i,j));
					}
				}
			}
			for(int l=0; l < pointsaround.Size(); l++)
			{
				STLTopEdge edge(i, pointsaround.Get(l,0),pointsaround.Get(l,1),pointsaround.Get(l,2));
				topedges.push_back(edge);
			}
		}
	}

	void STLTopology :: FindEdgesOfPoint()
	{
		topedgesperpoint.SetSize(points.size());
		for ( int i =0; i < topedges.size(); i++)
		{
			for( int j = 0; j<2; j++)
			{
				topedgesperpoint.Add(topedges[i][j],i);
			}
		}
	}

	void STLTopology :: FindNeighbourTrigs()
	{
		int index;
		for (int i = 0;i < topedges.size(); i++)
		{
			index=trias[topedges[i].TrigNum(0)].AnotherPointIndex(topedges[i][0],topedges[i][1]);
			trias[topedges[i].TrigNum(0)].NBTrigNum(index)= topedges[i].TrigNum(1);  
			index=trias[topedges[i].TrigNum(1)].AnotherPointIndex(topedges[i][0],topedges[i][1]);
			trias[topedges[i].TrigNum(1)].NBTrigNum(index)= topedges[i].TrigNum(0); 

		}
	}


	void STLTopology :: FindEdgesOfTrig()
	{
		for (int i = 0;i < topedges.size(); i++)
		{

			trias[topedges[i].TrigNum(0)].EdgeNum(trias[topedges[i].TrigNum(0)].AnotherPointIndex(topedges[i][0],topedges[i][1]))
				= i;  
			trias[topedges[i].TrigNum(1)].EdgeNum(trias[topedges[i].TrigNum(1)].AnotherPointIndex(topedges[i][0],topedges[i][1]))
				= i; 

		}
	}

	void STLTopology :: CalculateCosangleOfEdge()
	{
		for (int i = 0; i < topedges.size(); i++)
		{
			topedges[i].SetCosAngle((trias[topedges[i].TrigNum(0)].Normal()*trias[topedges[i].TrigNum(1)].Normal()));
		}
	}

	void STLTopology ::CalculateVectorOfPoint()
	{
		int i,k;
		int np = GetNP();
		normals.resize(GetNP());
		vector<int> normal_cnt(GetNP()); // counts number of added normals in a point

		for (i = 0; i < np; i++)
		{
			normal_cnt[i] = 0;
			normals[i] = Vec<3> (0,0,0);
		}

		for(i = 0; i < GetNT(); i++)
		{
			//      STLReadTriangle t = GetReadTriangle(i);
			//      STLTriangle st;

			Vec<3> n = GetTriangle(i).Normal ();

			for (k = 0; k < 3; k++)
			{
				int pi = GetTriangle(i).PNum(k);

				normal_cnt[pi]++;
				SetNormal(pi, GetNormal(pi) + n);
			}
		} 

		//normalize the normals
		for (i = 0 ;i < GetNP(); i++)
		{
			SetNormal(i,1./(double)normal_cnt[i]*GetNormal(i));
		}
	}

}


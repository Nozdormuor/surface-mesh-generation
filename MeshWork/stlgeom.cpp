#include "externalcommon.h"
#include "internalcommon.h"

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++   STL GEOMETRY   ++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

namespace meshwork
{

	STLGeometry :: STLGeometry()
	{
		;
	}

	STLGeometry :: ~STLGeometry()
	{
		;
	}
	void   STLGeometry ::Load (istream & ist)
	{
		//STLGeometry * geom = new STLGeometry();

		vector<STLReadTriangle> readtrigs;

		char buf[100];
		Point<3> pts[3];
		Vec<3> normal;

		int cntface = 0;
		int vertex = 0;
		bool badnormals = false;

		while (ist.good())
		{
			ist >> buf;

			int n = strlen (buf);
			for (int i = 0; i < n; i++)
				buf[i] = tolower (buf[i]);

			if (strcmp (buf, "facet") == 0)
			{
				cntface++;
			}

			if (strcmp (buf, "normal") == 0)
			{
				ist >> normal(0)
					>> normal(1)
					>> normal(2);
				normal.Normalize();
			}

			if (strcmp (buf, "vertex") == 0)
			{
				ist >> pts[vertex](0)
					>> pts[vertex](1)
					>> pts[vertex](2);

				vertex++;

				if (vertex == 3)
				{
					if (normal.Length() <= 1e-5)

					{
						normal = Cross (pts[1]-pts[0], pts[2]-pts[0]);
						normal.Normalize();
					}

					else

					{
						Vec<3> hnormal = Cross (pts[1]-pts[0], pts[2]-pts[0]);
						hnormal.Normalize();

						if (normal * hnormal < 0.5)
							badnormals = true;
					}

					vertex = 0;

					if ( (Dist2 (pts[0], pts[1]) > 1e-16) &&
						(Dist2 (pts[0], pts[2]) > 1e-16) &&
						(Dist2 (pts[1], pts[2]) > 1e-16) )

					{
						readtrigs.push_back (STLReadTriangle (pts, normal));

						if (readtrigs.size() % 100000 == 0)
							cout << readtrigs.size() << " triangles loaded" << endl;
					}
					else
					{
						cout << "Skipping flat triangle " 
							<< "l1 = " << Dist(pts[0], pts[1])
							<< ", l2 = " << Dist(pts[0], pts[2])
							<< ", l3 = " << Dist(pts[2], pts[1]) << endl;
					}

				}
			}
		}
		cout << readtrigs.size() << " triangles loaded" << endl;
		InitSTLGeometry(readtrigs);

	}
	void STLGeometry :: InitSTLGeometry(const vector<STLReadTriangle> & readtrias)
	{

		STLTopology::InitSTLGeometry(readtrias);
		RecoverBody();
	}

	void STLGeometry ::RecoverBody()
	{
		BuildEdges();
		//cout<<lines.size()<<endl;
		//for(int i = 0; i < lines.size(); i++)
		//{
		// for(int j = 0; j < 2; j++)
		// {
		//  cout<<lines[i].PNum(j)<<" , ";
		// }
		// cout<<endl;
		//}
		BuildSTLLines();
		//cout<<lines.size()<<endl;
		//for(int i = 0; i < lines.size(); i++)
		//{
		// for(int j = 0; j < lines[i].NP(); j++)
		// {
		//  cout<<lines[i].PNum(j)<<" , ";
		// }
		// cout<<endl;
		//}
		//BuildVertex();
		//cout<<vertices.size()<<endl;
		//for(int i = 0; i < vertices.size(); i++)
		//{
		//  cout<<vertices[i]<<endl;
		//}
		FindLinePoints();
		FindSpiralPoint();
		FindLineEndPoints();
		//cout<<spiralpoints.size()<<")))))))))))))"<<endl;
		//for(int i = 0; i < spiralpoints.size(); i++)
		//{
		// if(spiralpoints[i]==1)
		//  cout<<i<<endl;
		//}
		CalcFaceNums();
		BuildSTLFace();
		//cout<<faces.size()<<endl;
		//for(int i = 0; i < faces.size(); i++)
		//{
		//	for(int j = 0; j < faces[i].Triangles().size(); j++)
		//	{
		//		if(faces[i].GetTriangle(j)==78||faces[i].GetTriangle(j)==743||faces[i].GetTriangle(j)==822||faces[i].GetTriangle(j)==1487)
		//		{
		//			cout<<"trei"<<":";
		//		cout<<points[trias[faces[i].GetTriangle(j)][0]][0]<<" "<<points[trias[faces[i].GetTriangle(j)][0]][1]<<" "<<points[trias[faces[i].GetTriangle(j)][0]][2]<<endl;
		//		cout<<points[trias[faces[i].GetTriangle(j)][1]][0]<<" "<<points[trias[faces[i].GetTriangle(j)][1]][1]<<" "<<points[trias[faces[i].GetTriangle(j)][1]][2]<<endl;
		//		cout<<points[trias[faces[i].GetTriangle(j)][2]][0]<<" "<<points[trias[faces[i].GetTriangle(j)][2]][1]<<" "<<points[trias[faces[i].GetTriangle(j)][2]][2]<<endl;
		//		}

		//	}
		// cout<<"lines:"<<endl;
		// for(int j = 0; j < faces[i].Lines().size(); j++)
		// {
		//  cout<<faces[i].GetLine(j)<<" , ";
		// }
		// cout<<endl;
		//}
		FixGeometry();
		BuildVertex();
	}

	void STLGeometry ::BuildEdges()
	{
		FindEdgesFromAngles();
		BuildEdgesPerPoint();
	}

	void STLGeometry :: FindEdgesFromAngles()
	{
		for(int i = 0; i<topedges.size(); i++)
		{
			STLTopEdge &ed=topedges[i];
			if(ed.CosAngle()<0.866)
			{
				STLEdge se(ed.PNum(0),ed.PNum(1));
				se.SetLeftTrig(ed.TrigNum(0));
				se.SetRightTrig(ed.TrigNum(1));
				AddEdge(se);
			}
		}
	}

	void STLGeometry :: BuildEdgesPerPoint()
	{
		edgesperpoint.SetSize(GetNP());
		//add edges to points
		for (int i = 0; i < GetNE(); i++)
		{
			//cout << "EDGE " << GetEdge(i).PNum(1) << " - " << GetEdge(i).PNum(2) << endl;
			for (int j = 0; j < 2; j++)
			{
				AddEdgePP(GetEdge(i).PNum(j),i);
			}
		}
	}

	int STLGeometry :: IsEdge(int ap1, int ap2)
	{
		int i,j;
		for (i = 0; i < GetNEPP(ap1); i++)
		{
			for (j = 0; j < GetNEPP(ap2); j++)
			{
				if (GetEdgePP(ap1,i) == GetEdgePP(ap2,j)) {return 1;}
			}
		}
		return 0;
	}
	void STLGeometry :: BuildSTLLines()
	{
		linecnt=0;
		int chosededge(0);
		int startedge;
		int nextedge;
		int startpoint;
		int nextpoint;
		list<int> linepoints;
		list<double> linedists;
		list<int> lefttrigs;
		list<int> righttrigs;
		markededges.resize(GetNE());
		for(int i = 0; i < GetNE(); i++)
			markededges[i]=0;

		while(chosededge<GetNE())
		{
			linepoints.clear();
			linedists.clear();
			lefttrigs.clear();
			righttrigs.clear();
			for( chosededge; chosededge<GetNE(); chosededge++)
			{
				if(markededges[chosededge]==0)
				{
					break;
				}			
			}
			if(chosededge==GetNE())  break;
			linepoints.push_back(edges[chosededge].PNum(0));
			linepoints.push_back(edges[chosededge].PNum(1));
			linedists.push_back(Dist(points[edges[chosededge].PNum(0)],points[edges[chosededge].PNum(1)]));
			lefttrigs.push_back(edges[chosededge].LeftTrig());
			righttrigs.push_back(edges[chosededge].RightTrig());
			markededges[chosededge]=1;


			//find edge from front	
			startedge=chosededge;
			startpoint=edges[startedge].PNum(0);

			for(;;)
			{
				if(edgesperpoint.RowSize(startpoint)!=2)
					break;

				nextedge = FindAdjacentEdge(edges[startedge],startpoint);

				if(markededges[nextedge]!=0)
					break;

				nextpoint = edges[nextedge].PNum(1)+edges[nextedge].PNum(0)-startpoint;

				Vec<3> vec1=(points[edges[startedge].PNum(1)+edges[startedge].PNum(0)-startpoint]-points[startpoint]);
				vec1.Normalize();
				Vec<3> vec2=(points[edges[nextedge].PNum(1)+edges[nextedge].PNum(0)-startpoint]-points[startpoint]);
				vec2.Normalize();

				if(vec1*vec2>-0.35) 
					break;

				linepoints.push_front(nextpoint);
				linedists.push_front(Dist(points[nextpoint],points[startpoint]));
				lefttrigs.push_front(edges[nextedge].LeftTrig());
				righttrigs.push_front(edges[nextedge].RightTrig());
				markededges[nextedge]=1;

				startpoint=nextpoint;	
				startedge=nextedge;
			}

			//find edge from back	
			startedge=chosededge;
			startpoint=edges[startedge].PNum(1);

			for(;;)
			{
				if(edgesperpoint.RowSize(startpoint)!=2)
					break;

				nextedge = FindAdjacentEdge(edges[startedge],startpoint);

				if(markededges[nextedge]!=0)
					break;

				nextpoint = edges[nextedge].PNum(1)+edges[nextedge].PNum(0)-startpoint;

				Vec<3> vec1=(points[edges[startedge].PNum(1)+edges[startedge].PNum(0)-startpoint]-points[startpoint]);
				vec1.Normalize();
				Vec<3> vec2=(points[edges[nextedge].PNum(1)+edges[nextedge].PNum(0)-startpoint]-points[startpoint]);
				vec2.Normalize();

				if(vec1*vec2>-0.35) 
					break;

				linepoints.push_back(nextpoint);
				linedists.push_back(Dist(points[nextpoint],points[startpoint]));
				lefttrigs.push_back(edges[nextedge].LeftTrig());
				righttrigs.push_back(edges[nextedge].RightTrig());
				markededges[nextedge]=1;

				startpoint=nextpoint;		
				startedge=nextedge;
			}
			STLLine stlline;
			linecnt++;
			int size = linepoints.size();
			for(int i = 0; i < size; i++)
			{
				stlline.AddPoint(linepoints.front());
				linepoints.pop_front();
			}
			for(int i = 0; i < size-1; i++)
			{
				stlline.AddDist(linedists.front());
				linedists.pop_front();
			}
			for(int i = 0; i < size-1; i++)
			{
				stlline.AddLeftTrig(lefttrigs.front());
				lefttrigs.pop_front();
				stlline.AddRightTrig(righttrigs.front());
				righttrigs.pop_front();
			}
			AddLine(stlline);
			chosededge++;
		}
	}

	void STLGeometry :: BuildVertex()
	{
		for(int i = 0; i < lines.size(); i++)
		{
			if(lines[i].NP()==0)
				continue;
			if(IsLineClosed(lines[i]))
				continue;
			AddVertex(lines[i].EndP());
			AddVertex(lines[i].StartP());
		}
	}

	int STLGeometry ::AddVertex(int p)
	{
		for(int i= 0; i < vertices.size(); i++)
			if(vertices[i]==p)
				return -1;
		vertices.push_back((p));
		return vertices.size();
	}

	int STLGeometry :: FindAdjacentEdge(STLEdge &edge, int pnum) const
	{
		for(int i=0; i < edgesperpoint.RowSize(pnum); i++)
		{
			if(edges[edgesperpoint.Get(pnum,i)].PNum(0)!=edge.PNum(0)&&edges[edgesperpoint.Get(pnum,i)].PNum(0)!=edge.PNum(1))
				return edgesperpoint.Get(pnum,i);
			if(edges[edgesperpoint.Get(pnum,i)].PNum(1)!=edge.PNum(0)&&edges[edgesperpoint.Get(pnum,i)].PNum(1)!=edge.PNum(1))
				return edgesperpoint.Get(pnum,i);
		}
		return -1;
	}

	void STLGeometry :: CalcFaceNums()
	{
		int markedtrigs1 = 0;
		int starttrig(0);
		int laststarttrig = 0;
		int i, k, nnt;
		facecnt = 0;


		for (i = 0; i <GetNT(); i++)
			GetTriangle(i).SetFaceNum(-1);


		while (markedtrigs1 < GetNT())
		{
			for (i = laststarttrig; i < GetNT(); i++)
			{
				if (GetTriangle(i).GetFaceNum()==-1) 
				{
					starttrig = i;
					laststarttrig = i;
					break;
				}
			} 
			//add all triangles around starttriangle, which is reachable without going over an edge
			vector<int> todolist;
			vector<int> nextlist;
			facecnt++;
			markedtrigs1++;
			GetTriangle(starttrig).SetFaceNum(facecnt-1);
			todolist.push_back(starttrig);
			int ap1, ap2;

			while(todolist.size())
			{
				for (i = 0; i < todolist.size(); i++)
				{
					const STLTriangle& tt = GetTriangle(todolist[i]);
					for (k = 0; k < 3; k++)
					{
						nnt = NeighbourTrig(todolist[i],k);
						STLTriangle& nt = GetTriangle(nnt);
						if (nt.GetFaceNum()==-1)
						{
							tt.GetNeighbourPoints(nt,ap1,ap2);
							if (!IsEdge(ap1,ap2))
							{
								nextlist.push_back(nnt);
								nt.SetFaceNum(facecnt-1);
								markedtrigs1++;
							}
						}
					}
				}

				todolist.resize(0);
				for (i = 0; i < nextlist.size(); i++)
				{
					todolist.push_back(nextlist[i]);
				}
				nextlist.resize(0);	  
			}
		}
	}

	void STLGeometry ::AddTrigToFace(int nface, int trig)
	{
		faces[nface].Triangles().push_back(trig);
	}
	void STLGeometry ::AddLineToFace(int nface, int line)
	{
		faces[nface].Lines().push_back(line);
	}
	void STLGeometry :: BuildSTLFace()
	{
		faces.resize(facecnt);
		for(int i = 0; i < GetNT(); i++)
		{
			AddTrigToFace(GetTriangle(i).GetFaceNum(),i);
		}
		for(int i=0; i < GetNL(); i++)
		{
			AddLineToFace(GetTriangle(GetLine(i).GetLeftTrig(0)).GetFaceNum(),i);
			if(GetTriangle(GetLine(i).GetLeftTrig(0)).GetFaceNum()!=GetTriangle(GetLine(i).GetRightTrig(0)).GetFaceNum())
				AddLineToFace(GetTriangle(GetLine(i).GetRightTrig(0)).GetFaceNum(),i);
		}
	}
	void STLGeometry ::FindLinePoints()
	{
		linepoints.resize(GetNP());
		for(int i = 0;i<GetNP(); i++)
		{
			linepoints[i]=0;
		}
		for(int i = 0; i<GetNL(); i++)
		{
			for(int j = 0; j <lines[i].NP(); j++)
			{
				linepoints[lines[i].PNum(j)]=1;
			}
		}
	}
	void STLGeometry ::FindLineEndPoints()
	{
		lineendpoints.resize(GetNP());
		for(int i = 0;i<GetNP(); i++)
		{
			lineendpoints[i]=0;
		}
		for(int i = 0; i<GetNL(); i++)
		{
			if(lines[i].StartP()==lines[i].EndP())
				continue;
			lineendpoints[lines[i].StartP()]=1;
			lineendpoints[lines[i].EndP()]=1;
		}
	}
	int STLGeometry :: IsLineEndPoint(int pn) 
	{
		//  return 0;
		if (pn <0 || pn >= lineendpoints.size()) 
		{cout<<("Illegal pnum in IsLineEndPoint!!!"); return 0;}
		return lineendpoints[pn];
	}
	void STLGeometry ::FindSpiralPoint()
	{

		spiralpoints.resize(GetNP());
		for(int i = 0; i<GetNP();i++)
			spiralpoints[i]=0;
		for(int i = 0; i < GetNP(); i++)
		{
			for(int j = 0; j < trigsperpoint.RowSize(i); j++)
			{
				if(normals[i]*GetTriangle(trigsperpoint.Get(i,j)).Normal()<0.81)
				{
					AddSpiralPoint(i);
					continue;
				}
			}
		}
	}

	int STLGeometry ::AddSpiralPoint(int n)
	{
		for(int i = 0; i<linepoints.size();i++)
		{
			if(linepoints[n]==1)
				return -1;
		}
		spiralpoints[n]=1;
		return 1;
	}

	void STLGeometry ::FixGeometry()
	{
		ClearSmallFaces();
		//ClearSmallInteriorEdges();
	}

	void STLGeometry ::ClearSmallFaces(void)
	{
		vector<int> smallfaces;
		double area0=3.0*mparam.maxh*mparam.maxh;
		for(int i = 0; i < faces.size();i++)
		{
			if(AreaOfFace(faces[i]) < area0)
				smallfaces.push_back(i);
		}
		int size  = smallfaces.size();
		if(size == faces.size())
			size--;
		for(int i = 0; i < size; i++)
		{
			int f1 = GetFlatNeighbFace(smallfaces[i]);
			UnionFace(f1,smallfaces[i]);
		}

	}

	double STLGeometry ::AreaOfFace(STLFace &face)
	{
		double area = 0.;
		for(int i = 0; i < face.GetFaceTriangleNum();i++)
		{
			Vec<3> vec1 = points[trias[face.GetTriangle(i)][0]] - points[trias[face.GetTriangle(i)][1]];
			Vec<3> vec2 = points[trias[face.GetTriangle(i)][2]] - points[trias[face.GetTriangle(i)][1]];
			Vec<3> vec3 = Cross(vec1,vec2);
			area = area +vec3.Length ()/2.;
		}
		return area;
	}

	int STLGeometry ::GetFlatNeighbFace(int face)
	{
		double lenght = 0.;
		int num;
		for(int i = 0; i < faces[face].GetFaceLineNum();i++)
		{
			if(lines[faces[face].GetLine(i)].GetLength(points) > lenght)
			{
				for(int j = 0 ; j <faces.size(); j++)
				{
					for(int k = 0; k < faces[j].GetFaceLineNum();k++)
					{
						if(faces[face].GetLine(i)==faces[j].GetLine(k)&&face!=j)
						{
							lenght = lines[faces[face].GetLine(i)].GetLength(points);
							num = j;
						}
					}
				}
			}
		}
		return num;
	}

	void STLGeometry ::UnionFace( int f1 , int f2 )
	{
		for(int i = 0; i < faces[f2].GetFaceTriangleNum(); i++)
		{
			trias[faces[f2].GetTriangle(i)].SetFaceNum(f1);
			faces[f1].Triangles().push_back(faces[f2].GetTriangle(i));//Ð§ÂÊµÍÏÂ
		}
		for(int i = 0; i < faces[f2].GetFaceLineNum(); i++ )
		{
			int temp = 0;
			for(temp = 0; temp < faces[f1].GetFaceLineNum(); temp++)
			{
				if(faces[f2].GetLine(i)==faces[f1].GetLine(temp))
					break;
			}
			if(temp == faces[f1].GetFaceLineNum())
				faces[f1].Lines().push_back(faces[f2].GetLine(i));
		}
		faces[f2].Triangles().clear();
		faces[f2].Lines().clear();
	}

	void STLGeometry ::ClearSmallInteriorEdges()
	{
		for(int i=0; i < GetNL(); i++)
		{
			if(GetTriangle(GetLine(i).GetLeftTrig(0)).GetFaceNum()==GetTriangle(GetLine(i).GetRightTrig(0)).GetFaceNum()&&GetLine(i).GetLength(GetPoints())<1.5*mparam.maxh)
				GetLine(i).Clear();
		}
		for(int i = 0; i < lines.size();i++);

	}

	void STLGeometry :: RestrictLocalH(class Mesh & mesh, double gh , double gminh)
	{
		int i,j;

		int ap1,ap2,p3,p4;
		Point<3> p1p, p2p, p3p, p4p;
		Vec<3> n, ntn;
		double rzyl, localh;

		//  double localhfact = 0.5;
		// double geometryignorelength = 1E-4;

		Box<3> bb = GetBoundingBox();
		mesh.SetLocalH(bb.PMin() - Vec<3>(10, 10, 10),bb.PMax() + Vec<3>(10, 10, 10), mparam.grading);

		//mesh.SetGlobalH(gh);

		double mincalch = 1E10;
		double maxcalch = -1E10;

		double objectsize = bb.Diam();
		double geometryignoreedgelength = objectsize * 1e-5;

		if (stlparam.resthsurfcurvenable)
		{
			cout<<("Restrict H due to surface curvature")<<endl;

			vector<double> minh; 
			minh.resize(GetNP());
			for (i = 0; i < GetNP(); i++)
			{
				minh[i] = gh;
			}

			for (i = 0; i < GetNT(); i++)
			{
				const STLTriangle& trig = GetTriangle(i);
				//cout<<points[trig.PNum(0)]<<endl;
				//cout<<points[trig.PNum(1)]<<endl;
				//cout<<points[trig.PNum(2)]<<endl;
				n = GetTriangle(i).Normal();
				for (j = 0; j < 3; j++)
				{
					const STLTriangle& nt = GetTriangle(NeighbourTrig(i,j));
					//cout<<points[nt.PNum(0)]<<endl;
				    //cout<<points[nt.PNum(1)]<<endl;
				    //cout<<points[nt.PNum(2)]<<endl;
					trig.GetNeighbourPointsAndOpposite(nt,ap1,ap2,p3);	 

					//for(int i = 0; i < 3;i++)
					//{
					//	cout<<points[trias[nt.NBTrigNum(i)].PNum(0)]<<endl;
				    //    cout<<points[trias[nt.NBTrigNum(i)].PNum(1)]<<endl;
				    //    cout<<points[trias[nt.NBTrigNum(i)].PNum(2)]<<endl;
					//}

					if (IsEdge(ap1,ap2)) 
					{
						continue;
					}

					p4 = trig.PNum(0) + trig.PNum(1) + trig.PNum(2) - ap1 - ap2;

					p1p = GetPoint(ap1); p2p = GetPoint(ap2); 
					p3p = GetPoint(p3); p4p = GetPoint(p4);

					double h1 = GetDistFromInfiniteLine(p1p,p2p, p4p);
					double h2 = GetDistFromInfiniteLine(p1p,p2p, p3p);
					double diaglen = Dist (p1p, p2p);

					if (diaglen < geometryignoreedgelength)
						continue;
					rzyl = ComputeCylinderRadius 
						(n, GetTriangle (NeighbourTrig(i,j)).Normal(), 
						h1, h2);


					if (h1 < 1e-3 * diaglen && h2 < 1e-3 * diaglen)
						continue;

					if (h1 < 1e-5 * objectsize && h2 < 1e-5 * objectsize)
						continue;


					//	      rzyl = mindist/(2*sinang);
					localh = rzyl / stlparam.resthsurfcurvfac;
					if (localh < mincalch) {mincalch = localh;}
					if (localh > maxcalch) {maxcalch = localh;}
					if (localh < gh) 
					{
						minh[ap1] = min(minh[ap1],localh);
						minh[ap2] = min(minh[ap2],localh);
					}

					if (localh < gminh) {localh = gminh;}

					if(localh < objectsize)
						mesh.RestrictLocalHLine(p1p, p2p, localh);
					//cout << "restrict h along " << p1p << " - " << p2p << " to " << localh << endl;

					if (localh < 0.1)
					{
						localh = 0.1;
					}

				}
			}
		}


		if (stlparam.resthlinelengthenable)
		{
			//restrict h due to short lines
			cout<<("Restrict H due to line-length")<<endl;

			double minhl = 1E50;
			double linefact = 1./stlparam.resthlinelengthfac;
			double l;
			for (i = 0; i < GetNL(); i++)
			{

				l = GetLine(i).GetLength(points);

				if(l==0.)
					continue;
				const Point<3>& pp1 = GetPoint(GetLine(i).StartP());
				const Point<3>& pp2 = GetPoint(GetLine(i).EndP());

				if (l != 0.)
				{
					minhl = min(minhl,l*linefact);
					localh  = l*linefact;
					if (localh < gminh) {localh = gminh;}
					mesh.RestrictLocalH(pp1, localh);
					mesh.RestrictLocalH(pp2, localh);      
				}
			}
		}

		if (stlparam.resthedgeangleenable)
		{
			cout<<("Restrict H due to edges angle")<<endl;

			int lp1, lp2;
			Vec<3> v1,v2;
			mincalch = 1E50;
			maxcalch = -1E50;

			for (i = 0; i <GetNP(); i++)
			{

				if (GetNEPP(i) == 2 && !IsLineEndPoint(i))
				{
					if (GetEdge(GetEdgePP(i,0)).PNum(1) == GetEdge(GetEdgePP(i,1)).PNum(0) ||
						GetEdge(GetEdgePP(i,0)).PNum(0) == GetEdge(GetEdgePP(i,1)).PNum(1))
					{
						lp1 = 0; lp2 = 1;
					}
					else
					{
						lp1 = 1; lp2 = 0;
					}

					v1 = GetPoint(GetEdge(GetEdgePP(i,0)).PNum(1))-GetPoint(GetEdge(GetEdgePP(i,0)).PNum(0));
					v2 = GetPoint(GetEdge(GetEdgePP(i,1)).PNum(lp2))-GetPoint(GetEdge(GetEdgePP(i,1)).PNum(lp1));

					rzyl = ComputeCylinderRadius(v1, v2, v1.Length(), v2.Length());

					localh = rzyl / stlparam.resthedgeanglefac;
					if (localh < mincalch) {mincalch = localh;}
					if (localh > maxcalch) {maxcalch = localh;}

					if (localh < gminh) {localh = gminh;}

					if (localh != 0)
						mesh.RestrictLocalH(GetPoint(i), localh);
				}	  
			}
		}
		//need to be improved, can't get the perfect result!!!!
		if (stlparam.resthcloseedgeenable)
		{
			cout<<("Restrict H due to close edges")<<endl;

			double disttohfact = (1.0 / stlparam.resthcloseedgefac);
			int k,l;
			double h1, h2, h3,dist;
			int rc = 0;
			Point<3> p3p1;
			double mindist = 1E50;
			Vec<3> vec;

			Box3dTree* lsearchtree = new Box3dTree (GetBoundingBox().PMin() - Vec<3>(gh,gh,gh),
				GetBoundingBox().PMax() + Vec<3>(gh,gh,gh));

			vector<Point<3>> pmins(GetNL());
			vector<Point<3>> pmaxs(GetNL());

			double maxhline;
			//cout<<GetNL()<<endl;
			for (i = 0; i < GetNL(); i++)
			{
				maxhline = 0;
				STLLine l1 = GetLine(i);
				Point<3> pmin(GetPoint(l1.StartP())), pmax(GetPoint(l1.StartP())), px;

				for (j = 1; j < l1.NP(); j++)
				{
					px = GetPoint(l1.PNum(j));
					pmin.SetToMin (px);
					pmax.SetToMax (px);
				}
				maxhline = gh;
				Box<3> box(pmin,pmax);
				box.Increase(maxhline);

				lsearchtree->Insert (box.PMin(), box.PMax(), i);
				pmins[i] = box.PMin();
				pmaxs[i] = box.PMax();
			}

			vector<int> linenums;
			int k2;

			for (i = 0; i < GetNL(); i++)
			{

				linenums.resize(0);
				lsearchtree->GetIntersecting(pmins[i],pmaxs[i],linenums);

				STLLine l1 = GetLine(i);
				for (j = 0; j < l1.NP()-1; j++)
				{
					int nn = int (Dist(GetPoint(l1.PNum(j)),GetPoint(l1.PNum(j+1)))/gminh);
					//cout<<nn<<endl;
					for(int m = 0; m <=nn ; m++)
					{
						vec = (GetPoint(l1.PNum(j+1))-GetPoint(l1.PNum(j)));
						vec.Normalize();
						p3p1 = GetPoint(l1.PNum(j))+((double) m)*gminh*vec;
						h1 = mesh.GetH(p3p1);

						for (k2 = 0; k2 < linenums.size(); k2++)
						{
							k = linenums[k2];
							if (k == i) 
							{
								continue;
							} 
							/*  
							//old, without searchtrees
							for (k = i+1; k <= GetNLines(); k++)
							{
							*/
							STLLine l2 = GetLine(k);
							for (l = 1; l < l2.NP(); l++)
							{
								 Point<3> p3p2 = GetPoint(l2.PNum(l-1));
								 Point<3> p3p3 = GetPoint(l2.PNum(l));
								//h2 = (mesh.GetH(p3p2)*mesh.GetH(p3p2));
								//h3 = (mesh.GetH(p3p3)*mesh.GetH(p3p3));
								dist =GetDistFromLine(p3p2,p3p3,p3p1)*disttohfact;
								//cout<<dist<<endl;
								if (dist > 1E-12)
								{
									if (dist < h1) 
									{
										if (dist <gminh) 
										{
											dist = gminh;
										}
										mesh.RestrictLocalH(p3p1,dist); 
										rc++;
										mindist = min(mindist,dist);
									}
									// if (dist < h2) 
									//{
									//  mesh.RestrictLocalH(p3p2,sqrt(dist)); 
									//  rc++;
									//  mindist = min(mindist,sqrt(dist));
									//}
									// if (dist < h3) 
									// {
									//  mesh.RestrictLocalH(p3p3,sqrt(dist)); 
									//  rc++;
									//  mindist = min(mindist,sqrt(dist));
									// }
								}
							}
						}	  
					}
				}
			}
		}
	}


}
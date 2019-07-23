#include"externalcommon.h"
#include"internalcommon.h"


namespace meshwork
{
	int cw(const int i) 
	{
		return(i+2)%3;
	}

	int ccw(const int i) 
	{
		return(i+1)%3;
	}

	FaceCirculator & FaceCirculator::operator++()
	{
		list<Face>::iterator f=fi->GetNeighbor((index+1)%3);
	   index= f->Index(fi->GetVertex(index));
		fi=f;
		return *this;

	}

	FaceCirculator & FaceCirculator::operator--()
	{
		list<Face>::iterator f=fi->GetNeighbor((index+2)%3);
	   index= f->Index(fi->GetVertex(index));
		fi=f;
		return *this;
	}

	FaceCirculator  FaceCirculator::operator++(int)
	{
		FaceCirculator fc(fi,index);
		list<Face>::iterator f=fi->GetNeighbor((index+1)%3);
		 index= f->Index(fi->GetVertex(index));
		fi=f;
		return fc;
	}

	FaceCirculator  FaceCirculator::operator--(int)
	{
		FaceCirculator fc(fi,index);
		list<Face>::iterator f=fi->GetNeighbor((index+2)%3);
		 index= f->Index(fi->GetVertex(index));
		fi=f;
		return fc;
	}

	VertexCirculator & VertexCirculator::operator++()
	{
		int index1= f->Index(v0);
		int index2=f->Index(v);
		int index = 3 - index1-index2;
		v=f->GetVertex(index);
		f=f->GetNeighbor(index2);
		return *this;
	}

	VertexCirculator & VertexCirculator::operator--()
	{
		f=f->GetNeighbor(3-f->Index(v0)-f->Index(v));
		int index1= f->Index(v0);
		int index2=f->Index(v);
		int index = 3 - index1-index2;
		v=f->GetVertex(index);
		return *this;
	}

	VertexCirculator  VertexCirculator::operator++(int)
	{
		VertexCirculator vc(v0,v,f);
		int index1= f->Index(v0);
		int index2=f->Index(v);
		int index = 3 - index1-index2;
		v=f->GetVertex(index);
		f=f->GetNeighbor(index2);
		return vc;
	}

	VertexCirculator  VertexCirculator::operator--(int)
	{
		VertexCirculator vc(v0,v,f);
		f=f->GetNeighbor(3-f->Index(v0)-f->Index(v));
		int index1= f->Index(v0);
		int index2=f->Index(v);
		int index = 3 - index1-index2;
		v=f->GetVertex(index);
		return vc;
	}

	/* ************************************Vextex*************************************/
	Vertex::Vertex()
	{
		point=0.0;
		label=-1;
	}

	Vertex::Vertex(Point<3> poi)
	{
		point=poi;
		label=-1;
	}

	Vertex::Vertex(Point<3> poi,int lab)
	{
		point=poi;
		label=lab;
	}

	Vertex::Vertex(double x, double y, double z)
	{
		Point<3> poi(x,y,z);
		point=poi;
		label=-1;
	}

	Vertex::Vertex(double x, double y, double z,int lab)
	{
		Point<3> poi(x,y,z);
		point=poi;
		label=lab;
	}

	bool Vertex::operator==(const Vertex & vv)
	{
		return(point[0]==vv.point[0]&&point[1]==vv.point[1]&&point[2]==vv.point[2]&&label==vv.label);
	}

	int Vertex::SetFace(const Triangulation::FaceHandle fh)
	{
		incident_face =fh;
		return 1;
	}

	FaceCirculator  Vertex::IncidentFaces()
	{
		list<Face>::iterator f =incident_face;
		int i = f->Index(*this);
		return FaceCirculator(f,i);
	}

	VertexCirculator Vertex::IncidentVertices()
	{
		list<Face>::iterator f =incident_face;
		list<Vertex>::iterator v0= f->GetVertex( f->Index(*this));
		list<Vertex>::iterator v= f->GetVertex( cw(f->Index(*this)));
		return VertexCirculator(v0,v,f);
	}

	int  Vertex::Degree()
	{
		int count = 0;
		FaceCirculator vc = IncidentFaces(), done(vc);
		
			do
			{ 
				count += 1;
			} while ((++vc) != done);
		
		return count;
	}

	/* ************************************Face*************************************/
	Face::Face()
	{
		for(int i = 0; i < 3; i++)
		{
			edstatus[i]=0;
		}
	}

	Face:: Face( list<Vertex>::iterator v1,list<Vertex>::iterator v2, list<Vertex>::iterator v3)
	{
		vertex[0]=v1; vertex[1]=v2;vertex[2]=v3;
		for(int i = 0; i < 3; i++)
		{
			edstatus[i]=0;
		}
	}

	Face:: Face(list<Face>::iterator f1, list<Face>::iterator f2, list<Face>::iterator f3)
	{
		neighbor[0]=f1 ;neighbor[1]=f2; neighbor[2]=f3;
		for(int i = 0; i < 3; i++)
		{
			edstatus[i]=0;
		}
	}

	Face::Face(list<Vertex>::iterator v1,list<Vertex>::iterator v2, list<Vertex>::iterator v3,list<Face>::iterator f1, list<Face>::iterator f2, list<Face>::iterator f3)
	{
		vertex[0]=v1; vertex[1]=v2;vertex[2]=v3;
		neighbor[0]=f1 ;neighbor[1]=f2; neighbor[2]=f3;
		for(int i = 0; i < 3; i++)
		{
			edstatus[i]=0;
		}
	}
	//Face & Face::operator= (const Face & f2)
	//{
	//	vertex[0]=f2.vertex[0]; vertex[1]=f2.vertex[1];vertex[2]=f2.vertex[2];
	//	neighbor[0]=f2.neighbor[0] ;neighbor[1]=f2.neighbor[0]; neighbor[2]=f2.neighbor[2];
	//	edstatus[0]=f2.edstatus[0]; edstatus[1]=f2.edstatus[1]; edstatus[2]=f2.edstatus[2];
	//	infos=f2.info();
	//	return *this;
	//}

	bool Face::operator==(const Face& f2)
	{
		return(vertex[0]==f2.vertex[0]&& vertex[1]==f2.vertex[1]&&vertex[2]==f2.vertex[2]&&
			neighbor[0]==f2.neighbor[0] &&neighbor[1]==f2.neighbor[0]&& neighbor[2]==f2.neighbor[2]&&
			edstatus[0]==f2.edstatus[0]&& edstatus[1]==f2.edstatus[1]&& edstatus[2]==f2.edstatus[2]);
	}

	bool Face::IsVertex(const list<Vertex>::iterator vh)
	{
		for(int i= 0; i < 3; i++)
		{
			if(vertex[i]==vh)
				return true;
		}
		return false;
	}

	bool Face::IsNeighbor(const list<Face>::iterator fh)
	{
		for(int i = 0; i < 3;i++)
		{
			if(neighbor[i]==fh)
				return true;
		}
		return false;
	}

	int  Face::Index_f(const list<Face>::iterator fh)
	{
		for(int i = 0; i < 3;i++)
		{
			if(neighbor[i]==fh)
				return i;
		}
		return -1;
	}

	int  Face::Index_v(const list<Vertex>::iterator vh)
	{
		for(int i = 0; i < 3;i++)
		{
			if(vertex[i]==vh)
				return i;
		}
		return -1;
	}

	int Face::Index(const list<Face>::iterator fh)
	{
		for(int i = 0; i < 3;i++)
		{
			if(neighbor[i]==fh)
				return i;
		}
		return -1;
	}

	int Face::Index(const list<Vertex>::iterator vh)
	{
		for(int i = 0; i < 3;i++)
		{
			if(vertex[i]==vh)
				return i;
		}
		return -1;
	}

	int Face::Index(Vertex v)
	{
		for(int i = 0; i < 3;i++)
		{
			if(*vertex[i]==v)
				return i;
		}
		return -1;
	}

	/* ************************************Triangulation*************************************/
	Triangulation::Triangulation()
	{
		vertexs.clear();
		faces.clear();
		hasinfinite = false;
	}

	Triangulation::VertexHandle Triangulation::CreateVertex()
	{
		Vertex v;
		vertexs.push_back(v);
		return --vertexs.end();
	}

	Triangulation::VertexHandle Triangulation::CreateVertex(Point<3> poi)
	{
		Vertex v(poi);
		vertexs.push_back(v);
		return --vertexs.end();
	}

	Triangulation::VertexHandle Triangulation::CreateVertex(double x, double y, double z)
	{
		Vertex v(x,y,z);
		vertexs.push_back(v);
		return --vertexs.end();
	}

	Triangulation::VertexHandle Triangulation::CreateVertex(Vertex v)
	{
		vertexs.push_back(v);
		return --vertexs.end();
	}

	Triangulation::FaceHandle Triangulation::CreateFace()
	{
		Face fa;
		faces.push_back(fa);
		return --faces.end();
	}

	Triangulation::FaceHandle Triangulation::CreateFace(Face fa)
	{
		Face f=fa;
		faces.push_back(f);
		return --faces.end();
	}

	Triangulation::FaceHandle Triangulation::CreateFace(VertexHandle v1, VertexHandle v2, VertexHandle v3,FaceHandle f1, FaceHandle f2, FaceHandle f3)
	{
		Face f(v1,v2,v3,f1,f2,f3);
		faces.push_back(f);
		return --faces.end();
	}

	void Triangulation::DeleteFace(FaceHandle f)
	{
		faces.erase(f);
	}

	void Triangulation::DeleteVertex(VertexHandle v)
	{
		vertexs.erase(v);
	}

	void Triangulation::Flip(FaceHandle f, int i)
	{
		FaceHandle n  = f->GetNeighbor(i);
		int ni = MirrorIndex(f,i); 
		
		VertexHandle  v_cw = f->GetVertex(cw(i));
		VertexHandle  v_ccw = f->GetVertex(ccw(i));

		FaceHandle tr = f->GetNeighbor(ccw(i));
		int tri =  MirrorIndex(f,ccw(i));  
		FaceHandle bl = n->GetNeighbor(ccw(ni));
		int bli =  MirrorIndex(n,ccw(ni)); 

		f->SetVertex( n->GetVertex(ni),cw(i));
		n->SetVertex( f->GetVertex(i),cw(ni));

		SetAdjacency(f, i, bl, bli);
		SetAdjacency(f, ccw(i), n, ccw(ni));
		SetAdjacency(n, ni, tr, tri);

		v_cw->SetFace(n);
		v_ccw->SetFace(f);
	}

	Triangulation::VertexHandle Triangulation::InsertInFace(FaceHandle f)
	{
		VertexHandle  v = CreateVertex();

		VertexHandle v0 = f->GetVertex(0);
		VertexHandle v2 = f->GetVertex(2);
		VertexHandle v1 = f->GetVertex(1);

		FaceHandle n1 = f->GetNeighbor(1);
		FaceHandle n2 = f->GetNeighbor(2);

		FaceHandle f1 = CreateFace(v0, v, v2, f, n1, FaceHandle());
		FaceHandle f2 = CreateFace(v0, v1, v, f, FaceHandle(), n2);
        SetAdjacency(f1, 2, f2, 1);

			int i1 = MirrorIndex(f,1); //int i1 = n1->index(f);
			n1->SetNeighbor(f1,i1);
		
			int i2 = MirrorIndex(f,2);//int i2 = n2->index(f);
			n2->SetNeighbor(f2,i2);

		f->SetVertex(v, 0);
		f->SetNeighbor(f1, 1);
		f->SetNeighbor(f2, 2);

		v0->SetFace(f2); 
		v->SetFace(f);

		return v;

	}

	Triangulation::VertexHandle Triangulation::InsertInEdge(FaceHandle f, int i)
	{
		VertexHandle v;
		FaceHandle n = f->GetNeighbor(i);
		int in = MirrorIndex(f,i); //n->index(f);
		v = InsertInFace(f);
		v->Degree();
		Flip(n,in); 
		//v->Degree();
		return v;
	}

	void Triangulation::Clear()
	{
		//cout<<"start"<<endl;
		//for(FaceHandle f=faces.begin();f!=faces.end();f++)
		//{
		//	cout<<f->GetVertex(0)->GetPoint()<<endl;
		//	cout<<f->GetVertex(1)->GetPoint()<<endl;
		//	cout<<f->GetVertex(2)->GetPoint()<<endl;
		//}
		faces.clear();
		vertexs.clear();
		return;
	}

	bool Triangulation::IsInfinite(FaceHandle f) const
	{
		return hasinfinite&&f->IsVertex(infinitevertex);
	}

	bool Triangulation:: IsEdge(VertexHandle va, VertexHandle vb, FaceHandle &fr,  int & i) const
		// assume va is a vertex of t
		// returns true (false) if the line segment ab is (is not) an edge of t
		// if true is returned (fr,i) is the edge ab
		// with face fr on the right of a->b
	{
		FaceHandle fc = va->GetFace(); 
		FaceHandle start = fc;
	//	if (fc == 0) return false;
		int inda, indb;
		do 
		{
			inda=fc->Index(va);
			indb =  cw(inda);
			if(fc->GetVertex(indb) == vb)
			{
				fr=fc;
				i = 3 - inda - indb; //works in dim 1 or 2
				return true;
			}
			fc=fc->GetNeighbor(indb); //turns ccw around va
		} while (fc != start);
		return false;
	}

	bool Triangulation::IsEdge(VertexHandle va, VertexHandle vb) const
	{
		VertexCirculator vc = va->IncidentVertices(), done(vc);
		do {
			if( vb == vc.v ) {return true;} 
		} while (++vc != done);
		return false;
	}

	void Triangulation::RemoveDegree_3(VertexHandle v, FaceHandle f)
		// remove a vertex of degree 3
	{

		int i = f->Index(v);
		FaceHandle left = f->GetNeighbor(cw(i));
		int li = MirrorIndex(f,cw(i)); 
		FaceHandle right = f->GetNeighbor(ccw(i));
		int ri = MirrorIndex(f,ccw(i)); 

		FaceHandle ll, rr;
		VertexHandle q = left->GetVertex(li);
		if( left->GetVertex(li) != right->GetVertex(ri))
			cout<<"err"<<endl;

		ll = left->GetNeighbor(cw(li));
			int lli = MirrorIndex(left,cw(li)); 
			ll->SetNeighbor(f,lli);
		f->SetNeighbor( ll,cw(i));
		if (f->GetVertex(ccw(i))->GetFace() == left) f->GetVertex(ccw(i))->SetFace(f);    

		rr = right->GetNeighbor(ccw(ri));
			int rri =  MirrorIndex(right,ccw(ri)); //rr->index(right);
			rr->SetNeighbor(f,rri);	 
		f->SetNeighbor( rr,ccw(i));
		if (f->GetVertex(cw(i))->GetFace() == right) f->GetVertex(cw(i))->SetFace(f);  

		f->SetVertex(q,i);
		if (q->GetFace() == right || q->GetFace() == left) {
			q->SetFace(f);
		}
		DeleteFace(right);
		DeleteFace(left);

		DeleteVertex(v);
	} 

	VertexHandle Triangulation::JoinVertices(FaceHandle f, int i, VertexHandle v)
	{
  // this methods does the "join"-operation and preserves
  // the vertex v among the two vertices that define the edge (f, i) 

  VertexHandle v1 = f->GetVertex( ccw(i) );
  VertexHandle v2 = f->GetVertex( cw(i)  );

  if( v != v1 && v != v2 ) cout<<"err"<<endl;

  if ( v == v2 ) 
  {
    return JoinVertices(f->GetNeighbor(i), MirrorIndex(f,i), v);
  }

  int deg2 = v2->Degree();

  if( deg2 < 3 ) cout<<"err"<<endl;

  if ( deg2 == 3 )
  {
    RemoveDegree_3(v2, f->GetNeighbor(ccw(i)));
    return v1;
  }
  
  /*
  // The following drawing corrsponds to the variables
  // used in this part...
  // The vertex v1 is returned...
  //
  //      itl       i=v3      itr
  //       *---------*---------*
  //        \       / \       /
  //         \  tl /   \  tr /
  //          \   /  f  \   /
  //           \ /       \ /
  //  v1=ccw(i) *---------*  cw(i)=v2
  //           / \       / \
  //          /   \  g  /   \
  //         /  bl \   /  br \
  //        /       \ /	      \
  //       *---------*---------*
  //      ibl       j=v4      ibr
  //                                                           
  // The situation after the "join"-operation is as follows:
  //
  //                 i
  //           *-----*-----*
  //            \    |    /
  //             \ tl|tr /
  //              \  |  /
  //               \ | /
  //                \|/
  //                 *  v1
  //                /|\
  //               / | \
  //              /  |	\
  //             / bl|br \
  //            /    |	  \
  //           *-----*-----*
  //
  */

  // first we register all the needed info
  FaceHandle g = f->GetNeighbor(i);
  int j = MirrorIndex(f,i);

  FaceHandle tl = f->GetNeighbor( cw(i)  );
  FaceHandle tr = f->GetNeighbor( ccw(i) );

  int itl = MirrorIndex(f, cw(i)  );
  int itr = MirrorIndex(f, ccw(i) );

  FaceHandle bl = g->GetNeighbor( ccw(j) );
  FaceHandle br = g->GetNeighbor( cw(j)  );

  int ibl = MirrorIndex(g, ccw(j) );
  int ibr = MirrorIndex(g, cw(j)  );

  // we need to store the faces adjacent to v2 as well as the
  // indices of v2 w.r.t. these faces, so that afterwards we can set 
  // v1 to be the vertex for these faces
  std::vector<FaceHandle> star_faces_of_v2;
  std::vector<int> star_indices_of_v2;
  FaceCirculator fc_start(v2->IncidentFaces());
  FaceCirculator fc = fc_start;

  do 
  {
    FaceHandle ff(fc.fi);
    star_faces_of_v2.push_back(ff);
    star_indices_of_v2.push_back(ff->Index(v2));
    ++fc;
  } while ( fc != fc_start );

  if( int(star_faces_of_v2.size()) != deg2 ) cout<<"err"<<endl;

  // from this point and on we modify the values

  // first set the neighbors
  SetAdjacency(tl, itl, tr, itr);
  SetAdjacency(bl, ibl, br, ibr);

  // make sure that all the faces containing v2 as a vertex, now
  // contain v1
  for (unsigned int k = 0; k < star_faces_of_v2.size(); k++)
  {
    int id = star_indices_of_v2[k];
    if( star_faces_of_v2[k]->GetVertex(id) != v2 )  cout<<"err"<<endl;
    star_faces_of_v2[k]->SetVertex(  v1,id );
  }

  // then make sure that all the vertices have correct pointers to 
  // faces
  VertexHandle v3 = f->GetVertex(i);
  VertexHandle v4 = g->GetVertex(j);
  if ( v3->GetFace() == f )  v3->SetFace(tr);
  if ( v4->GetFace() == g )  v4->SetFace(br);
  if ( v1->GetFace() == f || v1->GetFace() == g ) v1->SetFace(tl);


  // memory management
  star_faces_of_v2.clear();
  star_indices_of_v2.clear();

  DeleteFace(f);
  DeleteFace(g);

  DeleteVertex(v2);
  return v1;
}

	void Triangulation::SetAdjacency(FaceHandle f0, int i0, FaceHandle f1, int i1) const
	{
		f0->SetNeighbor(f1,i0);
		f1->SetNeighbor(f0,i1);
	}

	int Triangulation::MirrorIndex(FaceHandle f, int i) const
	{
		return ccw( f->GetNeighbor(i)->Index(f->GetVertex(ccw(i))));
	}

	//test
	void Triangulation::OutputDXF( char *filename )
	{
		ofstream outfile(filename);
		outfile<<"0\nSECTION\n";
		outfile<<"2\nENTITIES\n";
		for(AllFacesIterator i=AllFacesBegin();i!=AllFacesEnd();i++)
		{
			if(IsInfinite(i)) continue;
			outfile<<"0\n3DFACE\n";
			outfile<<"8\n0\n";
			for(int j=0;j<3;j++)
			{
				Point<3> p = i->GetVertex(j)->GetPoint();
				outfile<<10+j<<"\n"<<p[0]<<"\n"<<20+j<<"\n"<<p[1]<<"\n"<<30+j<<"\n"<<p[2]<<"\n";
				if(j==2)
					outfile<<13<<"\n"<<p[0]<<"\n"<<23<<"\n"<<p[1]<<"\n"<<33<<"\n"<<p[2]<<"\n";
			}
		}
		outfile<<"0\nENDSEC\n";
		outfile<<"0\nEOF";
		outfile.close();
	}

	void Triangulation::OutputSTL( char *filename )
	{
		//ofstream outfile(filename);
		//outfile<<"0\nSECTION\n";
		//outfile<<"2\nENTITIES\n";
		//for(AllFacesIterator i=AllFacesBegin();i!=AllFacesEnd();i++)
		//{
		//	if(IsInfinite(i)) continue;
		//	outfile<<"0\n3DFACE\n";
		//	outfile<<"8\n0\n";
		//	for(int j=0;j<3;j++)
		//	{
		//		Point<3> p = i->GetVertex(j)->GetPoint();
		//		outfile<<10+j<<"\n"<<p[0]<<"\n"<<20+j<<"\n"<<p[1]<<"\n"<<30+j<<"\n"<<p[2]<<"\n";
		//		if(j==2)
		//			outfile<<13<<"\n"<<p[0]<<"\n"<<23<<"\n"<<p[1]<<"\n"<<33<<"\n"<<p[2]<<"\n";
		//	}
		//}
		//outfile<<"0\nENDSEC\n";
		//outfile<<"0\nEOF";
		//outfile.close();
		ofstream  outfile(filename);

		int i;

		outfile.precision(10);

		outfile << "solid" << endl;

		for(AllFacesIterator i=AllFacesBegin();i!=AllFacesEnd();i++)
		{
			if(IsInfinite(i)) continue;
			outfile << "facet normal ";
			const Point<3>& p1 = i->GetVertex(0)->GetPoint();
			const Point<3>& p2 =i->GetVertex(1)->GetPoint();
			const Point<3>& p3 =i->GetVertex(2)->GetPoint();

			Vec<3> normal = Cross(p2-p1,p3-p1);
			if (normal.Length() != 0)
			{
				normal /= (normal.Length());
			}

			outfile << normal[0]<< " " << normal[1] << " " << normal[2] << "\n";
			outfile << "outer loop\n";

			outfile << "vertex " << p1[0] << " " << p1[1] << " " << p1[2] << "\n";
			outfile << "vertex " << p2[0] << " " << p2[1] << " " << p2[2] << "\n";
			outfile << "vertex " << p3[0] << " " << p3[1] << " " << p3[2] << "\n";

			outfile << "endloop\n";
			outfile << "endfacet\n";
		}
		outfile << "endsolid" << endl;
	}

}
#include"externalcommon.h"
#include"internalcommon.h"
namespace meshwork
{
	int MW_STLSurfaceMeshing (STLGeometry & geom,class Mesh & mesh)
	{
		double h = mparam.maxh;
		double gminh = mparam.minh;
		if(mparam.minh==0)
			gminh = mparam.maxh/100.;
		vector<vector<int>> meshsurfaces;
		vector<STLLine> meshlines;
		vector<Point<3>> meshpoints;
		STLSurfaceMesh *pSFM,sfm;
		pSFM=&sfm;
		//限制剖分尺寸
		geom.RestrictLocalH(mesh, h, gminh);

		//填充线段
		for (int i = 0; i < geom.GetNL(); i++)
		{
			meshlines.push_back(geom.GetLine(i).Mesh(geom.GetPoints(), meshpoints, h, mesh)); 
		}

		//填充面
		for (int i = 0; i < geom.GetNF(); i++)
		{
			meshsurfaces.push_back(geom.GetFace(i).Mesh(geom.GetPoints(), geom.GetTriangles(),meshpoints,meshlines, h, mesh)); 
		}
		//将点插入mesh类
		POINTTYPE ptype;
		for(int i = 0;i < meshpoints.size();i++)
		{
			ptype=EDGEPOINT;
			MeshPoint mp(meshpoints[i],ptype);
			mesh.AddPoint(mp);
		}
		//形成面单元
		for (int i = 0; i< meshsurfaces.size();i++)
		{
			//cout<<"ddd"<<endl;
			//surfremesh *pSFM=new surfremesh();
			if(geom.GetFace(i).GetFaceTriangleNum()==0)
				continue;
			TransFromTripsToCGAL(geom.GetPoints(), geom.GetTriangles(),geom.GetFace(i).Triangles(),pSFM);
			TriangulationInKdtree *fvkdt = new TriangulationInKdtree(pSFM);
			SetBoundPropOnGEdge(fvkdt,pSFM,geom.GetPoints(),geom.GetLines(),geom.GetFace(i).Lines());
			//cout<<"mmmm"<<endl;
			MeshingWithNewNodes(pSFM, fvkdt,meshpoints,meshsurfaces[i],mesh);
			delete fvkdt;
			//cout<<"dfdf"<<endl;
			//提取面单元

			{
				FaceDescriptor fd(i);
				for(Triangulation::AllFacesIterator fi=pSFM->sft.AllFacesBegin();fi!=pSFM->sft.AllFacesEnd();fi++)
				{
					if(skipBeInfinite(fi,pSFM->sft)) continue;
					Element2d ele;
					ele.SetIndex(i);
					for(int k=0; k<3; k++)
					{

						if(fi->Info().eproperty[k]==pSFM->ptonbd) 
						{
							ele.SetEdgeType(k,1);
						}
						else
						{
							ele.SetEdgeType(k,0);
						}
						/*if(fi->info().eproperty[(k+1)%3]==pSFM->ptonbd||fi->info().eproperty[(k+2)%3]==pSFM->ptonbd) 
						{
						ptype=EDGEPOINT;
						}
						else
						{
						ptype=INNERPOINT;
						}
						MeshPoint mp(fi->GetVertex(k)->GetPoint(),ptype);*/
						ele.PNum(k)=fi->GetVertex(k)->Info().status;//mesh.AddPoint(mp);
					}

					fd.AddElement(mesh.AddElement2d(ele));
				}	
				mesh.AddFaceDescriptor(fd);
			}
			//delete pSFM;
		}

		return 1;
	}

	void TransFromTripsToCGAL(const vector<Point<3> >& ap,const vector<STLTriangle> &tris,const vector<int> trips,STLSurfaceMesh *pSFM)
	{
		vector <int> pois;
		int i , j , k;
		for( i = 0; i < trips.size(); i++)
		{
			for( j = 0; j < 3; j++)
			{
				int temp = pois.size();
				for( k = 0; k < pois.size(); k++)
				{
					if(tris[trips[i]][j]==pois[k])
						break;
				}
				if(k==temp)
					pois.push_back(tris[trips[i]][j]);
			}
		}

		vector <Point<3>> fpois;
		for(int i = 0; i < pois.size(); i++)
			fpois.push_back(ap[pois[i]]);

		vector<STLTriangle> ftrips;
		ftrips.resize(trips.size());
		int m;
		for(int i = 0; i < trips.size(); i++)
		{
			for(int j = 0; j <3; j++)
			{
				for(int k = 0; k < pois.size(); k++)
				{
					if(tris[trips[i]][j]==pois[k])
					{
						m=k;
						break;
					}
				}
				ftrips[i][j]=m;
			}
		}

		for(int i = 0; i< ftrips.size(); i++)
			for(int j = 0; j < 3; j++)
				ftrips[i].NBTrigNum(j)=-1;

		// Creation of the vertices
		if(pSFM->sft.NumberOfVertices() != 0) pSFM->sft.Clear();
		vector<Triangulation::VertexHandle > V(ap.size());
		vector<Triangulation::FaceHandle> F(tris.size());
		for(int i=0 ; i < pois.size(); ++i)
		{
			V[pois[i]] = pSFM->sft.CreateVertex();
			V[pois[i]]->SetPoint(Point<3>(ap[pois[i]][0],ap[pois[i]][1],ap[pois[i]][2]));
		}
		// Creation of the faces
		int index;
		for(int i = 0; i < trips.size(); ++i) 
		{
			F[trips[i]] = pSFM->sft.CreateFace() ;
			for(int j = 0; j < 3 ; ++j)
			{
				index=tris[trips[i]][j];
				F[trips[i]]->SetVertex( V[index],j);
				V[index]->SetFace(F[trips[i]]);
			}
		}
		//cout<<pSFM->sft.NumberOfVertices();
		//pSFM->sft.OutputSTL("test.stl");
		//create infinite vertex and face 1
		int (*tneighb)[3]=(int (*)[3]) new int[3*trips.size()];
		for(i=0; i<trips.size(); i++)
		{
			for(int j=0; j<3; j++)
				tneighb[i][j]=-1;
		}
		for(int i = 0; i < trips.size(); ++i)
		{
			for(int j = 0; j < 3; ++j)
			{
				index=tris[trips[i]].NBTrigNum(j);
				for(k= 0; k<trips.size(); k++)
				{
					if(index==trips[k])
					{
						F[trips[i]]->SetNeighbor(F[index],j);
						tneighb[i][j]=k;
					}
				}
			}
		}
		//create infinite vertex and face 2
		int nfinf=0;
		for(int i=0; i<trips.size(); i++)
			for(int j=0; j<3; j++)
				if(tneighb[i][j]==-1) nfinf++;
		if(nfinf==0) 
		{
			delete [] tneighb;
			return;
		}
		VertexHandle v0=pSFM->sft.CreateVertex();	
		std::vector<Triangulation::FaceHandle> finf(nfinf);
		int nf=0;
		FaceHandle fa;
		for(int i=0; i<trips.size(); i++)
		{
			for(int j=0; j<3; j++)
			{
				if(tneighb[i][j]!=-1) continue;
				finf[nf]=fa=pSFM->sft.CreateFace() ;
				fa->SetVertex(v0,0);
				fa->SetVertex(V[tris[trips[i]][(j+2)%3]],1);
				fa->SetVertex(V[tris[trips[i]][(j+1)%3]],2);
				fa->SetNeighbor(F[trips[i]],0);
				F[trips[i]]->SetNeighbor(fa,j);
				nf++;
			}
		}
		v0->SetFace(finf[0]);
		pSFM->sft.SetInfiniteVertex(v0);	
		for(int i=0; i<nfinf; i++)
		{
			FaceHandle fnei;
			FaceHandle fa = finf[i];
			fnei=fa->GetNeighbor(0);
			FaceHandle f0=fa;
			for(;;)
			{
				if(pSFM->sft.IsInfinite(fnei)) break;
				int ind=fnei->Index_f(f0);
				f0=fnei;
				fnei=fnei->GetNeighbor((ind+2)%3);
			}
			finf[i]->SetNeighbor(fnei,1);
			fnei->SetNeighbor(finf[i],2);
		}
		delete [] tneighb;
		//cout<<pSFM->sft.number_of_vertices()<<endl;
		//cout<<pSFM->sft.number_of_faces()<<endl;
		//pSFM->sft.OutputSTL("test.stl");

		//pSFM->sft.output_DXF("test.dxf");
	}

	bool SetBoundPropOnGEdge(TriangulationInKdtree *fvkdta,STLSurfaceMesh *psfma,const vector<Point<3> >& ap,const vector<STLLine> &lines,const vector<int> &facelines)
	{
		Triangulation::VertexHandle vh,vpre;
		Triangulation::AllFacesIterator fi;
		for(fi=psfma->sft.AllFacesBegin();fi!=psfma->sft.AllFacesEnd();fi++)
		{
			for(int i=0; i<3; i++)
			{
				fi->Info().eproperty[i]=0;
			}
		}
		for(int i = 0; i < facelines.size(); i++)
		{
			if(lines[facelines[i]].NP()==0)
				continue;
			fvkdta->GetVertexhandleFromPoint(ap[lines[facelines[i]].PNum(0)],vh);
			for(int j = 1; j < lines[facelines[i]].NP(); j++)
			{
				vpre=vh;
				fvkdta->GetVertexhandleFromPoint(ap[lines[facelines[i]].PNum(j)],vh);
				psfma->SetEdgeProperty(Edge2v(vh,vpre),psfma->ptonbd);
			}

		}

		return true;
	}

	void MeshingWithNewNodes(STLSurfaceMesh *pSFM,  TriangulationInKdtree *fvt,vector<Point<3>> &meshpoints,vector<int>& meshsurfaces,class Mesh & mesh)
	{
		//pSFM->sft.OutputSTL("test.stl");
		extern bool skipBeInfinite(VertexHandle &vt,Triangulation &sft);


		pSFM->InitializeAllVertexFlag(); //set all vertices as old
		for(int i = 0; i < meshsurfaces.size(); i++)
		{
			Point<3> pt;
			for(int j = 0; j < 3; j++)
			{
				pt[j]=meshpoints[meshsurfaces[i]][j];
			}
			double dist;
			bool flag=true;
			Triangulation::FaceHandle fh, *pFH=&fh;
			if(!fvt->GetFacehandleFromPoint(pt,*pFH,&dist))
			{
				cout<<"err here"<<endl;
			}
			Edge edge;
			double sqdist;
			double r=mesh.GetH(pt);
			bool tooClose=false;
			Point<3> pa;//,pb[3];
			double sqd,sqd0;
			sqd0=std::numeric_limits<double>::max(),sqd=0;
			int k0=-1;
			for(int k=0; k<3; k++)
			{
				pSFM->Get1VertexPtOfTdsTri(*pFH,k,pa);
				if((sqd=(pt-pa).Length2())<sqd0)
				{
					k0=k;
					sqd0=sqd;
				}
			}
			if(sqd0<0.00001*r*r)
			{//1.e-5){
				//pSFM->setFlagOfTdsVertex((*pFH)->vertex(k),surfremesh:: accept);
				if((*pFH)->GetVertex(k0)->Info().status>=0)
				{
					cout<<"err1 "<<(*pFH)->GetVertex(k0)->Info().status<<" "<<i<<endl;
					cout<<pt[0]<<" "<<pt[1]<<" "<<pt[2]<<endl;
					cout<<(*pFH)->GetVertex(k0)->GetPoint()[0]<<" "
						<<(*pFH)->GetVertex(k0)->GetPoint()[1]<<" "<<(*pFH)->GetVertex(k0)->GetPoint()[2]<<endl;
				}
				(*pFH)->GetVertex(k0)->Info().status=meshsurfaces[i];
				tooClose=true;
				//		break;
			}
			if(tooClose) continue;
			VertexHandle vh;
			pSFM->FindTheNearestEdgeOfTdsTriToPtWithSqdist(pt,*pFH,edge,sqdist);
			if(sqdist>0.001*r*r)
			{
				vh=Insert1PointInInteriorOfTriangle(pt,pFH,*pSFM,fvt,flag);
				//localImpAfterIns(v);
			}else vh=Insert1PointOnEdge(pt,pFH,edge,*pSFM,fvt,flag);

			LocalImproveAfterInsert(vh,*pSFM,fvt);
			vh->Info().bkvaddress=fvt->InsertVertexHandle(vh);
			vh->Info().radius=r;
			vh->Info().haveR=true;
			if(vh->Info().status>=0) cout<<"err6 "<<vh->Info().status<<endl;

			vh->Info().status=meshsurfaces[i];
		}
		//pSFM->sft.OutputSTL("test.stl");
		//cout<<"fsdgdfgdg"<<endl;
		int countDelete=0;
		int undelhd=0;
		int undelfour=0;
		Triangulation::AllVerticesIterator vit,vitnext;
		std::queue<VertexHandle> qvtaround,*pqvt=&qvtaround;
		vit=pSFM->sft.AllVerticesBegin();
		while(vit!=pSFM->sft.AllVerticesEnd())
		{
			vitnext=vit;
			vitnext++;
			if(!skipBeInfinite(vit,pSFM->sft)&&!pSFM->IsVertexRemain(vit) && pSFM->CanVertexDel(vit))
			{
				pSFM->DeleteOneOldNode(vit,pqvt);
				countDelete++;

			}
			pSFM->LocalImpAfter1VtDel(pSFM->qedgeforimp,qvtaround);
			vit=vitnext;
			//pSFM->sft.OutputSTL("test.stl");
		}
		//pSFM->sft.OutputSTL("test.stl");
		//countDelete=0;
		//for(Triangulation::AllVerticesIterator vit=pSFM->sft.AllVerticesBegin();vit!=pSFM->sft.AllVerticesEnd();vit++)
		//{
		//	if(skipBeInfinite(vit,pSFM->sft)) continue;
		//	if(!pSFM->IsVertexRemain(vit)) countDelete++;
		//}
		//cout<<"before IMP remain "<<countDelete;//<<endl;
		pSFM->ImpSurfaceMesh(pSFM->qedgeforimp);
		//pSFM->sft.OutputSTL("test.stl");
		vit=pSFM->sft.AllVerticesBegin();
		while(vit!=pSFM->sft.AllVerticesEnd())
		{
			vitnext=vit;
			vitnext++;
			if(!skipBeInfinite(vit,pSFM->sft)&&!pSFM->IsVertexRemain(vit) && pSFM->CanVertexDel(vit))
			{
				pSFM->DeleteOneOldNode(vit,pqvt);
				countDelete++;

			}
			pSFM->LocalImpAfter1VtDel(pSFM->qedgeforimp,qvtaround);
			vit=vitnext;
			//pSFM->sft.OutputSTL("test.stl");
		}
		countDelete=0;
		for(Triangulation::AllVerticesIterator vit=pSFM->sft.AllVerticesBegin();vit!=pSFM->sft.AllVerticesEnd();vit++)
		{
			if(skipBeInfinite(vit,pSFM->sft)) continue;
			if(!pSFM->IsVertexRemain(vit)) countDelete++;
		}
		cout<<" remain "<<countDelete;//<<endl;
	}

	VertexHandle Insert1PointInInteriorOfTriangle(Point<3> &pt,FaceHandle *pFH,STLSurfaceMesh &sfm,TriangulationInKdtree *fvt,bool flag)
	{
		VertexHandle v;
		FaceHandle tridiv=*pFH;
		//poly.polytree->deleteExinfo(tridiv->info().bkaddress,1);
		//delete tridiv->info().bkaddress;
		fvt->RemoveFaceHandle(tridiv);
		//	tridiv->info().bkaddress=0;
		//	if(!flag) delete pFH;
		sfm.Insert1PtInInteriorOfTdsTri(pt,tridiv,v);
		sfm.SetFlagOfTdsVertex(v,STLSurfaceMesh:: newone);
		FaceCirculator fc= sfm.sft.Incident_faces(v);
		FaceCirculator end=fc;
		//	do{
		//		if(skipBeInfinite(Face_handle(fc),sfm.sft)) jf_error("err ins1ptinerioroftdstri");//continue;	
		//		Triangulation::Face_handle fh=(Triangulation::Face_handle)fc;
		//		fh->info().bkaddress=fvt->insertFaceHandle(fh);
		//fc;
		//	}while(++fc!=end); will insert in fvt in the final stage.
		return v;
	}

	VertexHandle Insert1PointOnEdge(Point<3> &pt,FaceHandle *pFH,Edge edge,STLSurfaceMesh &sfm,TriangulationInKdtree *fvt,bool flag)
	{
		extern void PointProjectTo2PLine(double p[3] ,double p0[3] ,double p1[3] ,double pb[3]);
		extern bool skipBeInfinite(FaceHandle &fa,Triangulation &sft);

		VertexHandle v;
		FaceHandle tri;
		int i;
		double proj[3]={pt[0],pt[1],pt[2]};
		Point<3> pa,pb;
		tri=edge.first;
		i=edge.second;
		sfm.Get2VertexPtOfTdsEdge(edge,pa,pb);
		double ppt[3],ppa[3],ppb[3];
		for(int i = 0; i < 3; i++)
		{
			ppt[i]=pt[i];
			ppa[i]=pa[i];
			ppb[i]=pb[i];
		}
		PointProjectTo2PLine(ppt,ppa,ppb,proj);
		//	double r=radius_node(pt);
		//	if(SqDistance3D(proj,pa)<0.001*r*r){jf_error("once "); return tri->vertex((i+1)%3);}
		//	if(SqDistance3D(proj,pb)<0.001*r*r){jf_error("once "); return tri->vertex((i+2)%3);}
		//sfm.sft.clear();
		if(!skipBeInfinite(tri,sfm.sft))
			fvt->RemoveFaceHandle(tri);
		//sfm.sft.clear();
		if(!skipBeInfinite(tri->GetNeighbor(i),sfm.sft))
			fvt->RemoveFaceHandle(tri->GetNeighbor(i));

		//deleteExinfo(tri->info().bkaddress,1);
		//fvt->deleteExinfo(tri->neighbor(i)->info().bkaddress,1);
		//delete tri->info().bkaddress;
		//delete tri->neighbor(i)->info().bkaddress;
		sfm.Insert1PtOnTdsEdge(proj,edge,v);
		sfm.SetFlagOfTdsVertex(v,STLSurfaceMesh:: newone);
		if(!flag) delete pFH;
		FaceCirculator fc= sfm.sft.Incident_faces(v);
		FaceCirculator end=fc;
		//	do{
		//		if(skipBeInfinite(Face_handle(fc),sfm.sft)) continue;	
		//		Triangulation::Face_handle fh=(Triangulation::Face_handle)fc;
		//		fh->info().bkaddress=fvt->insertFaceHandle(fh);

		//	}while(++fc!=end); //2015-7-20, will insert in fvt in the final stage.

		return v;
	}

	void LocalImproveAfterInsert( VertexHandle &vt,STLSurfaceMesh &sfm, TriangulationInKdtree *fvt)
	{ //not used right now

		//extern bool skipBeInfinite(VertexHandle &vt,Triangulation &sft);
		extern bool skipBeInfinite(FaceHandle &fa,Triangulation &sft);
		int i=0;
		//	int a;
		FaceCirculator vttri=sfm.sft.Incident_faces(vt),done(vttri),nexttri;
		//	if(vttri!=0)
		//		do{
		//			tri[i]=vttri;
		//			i++;   //已修改，加入i++
		//		}while(++vttri!=done);
		//	if(fvt==0)
		do{
			int a=vttri.fi->Index(vt);
			FaceHandle neighbtri=vttri.fi->GetNeighbor(a);
			nexttri=vttri;
			nexttri++;
			if( sfm.CanTdsEdgeSwap(Edge(vttri.fi,a),true)&&sfm.IsImprovedAfterEdgeSwap(Edge(vttri.fi,a)) )
			{
				if(fvt) fvt->RemoveFaceHandle(neighbtri);
				sfm.SwapEdge(vttri.fi,a);
				done=nexttri;
				done--;
			}
		}while((vttri=nexttri)!=done);

		do
		{
			if(skipBeInfinite(vttri.fi,sfm.sft)) continue;	
			FaceHandle fh=vttri.fi;
			if(fvt) fh->Info().bkaddress=fvt->InsertFaceHandle(fh);
		}while(++vttri!=done); 

	}


}
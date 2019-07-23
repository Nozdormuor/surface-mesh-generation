#include"externalcommon.h"
#include"internalcommon.h"
namespace meshwork
{
	extern 	int cw(const int i) ;
	extern 	int ccw(const int i) ;
	extern void copyPtFromVersionStructToArray(Point<3> p, Point<3> &pt);
	extern double DistPointToSegm( double p[3] , double ps[3] , double pe[3] ,double pnear[3],int *ip);
	extern void normal_3p( double *p1 , double *p2 ,double *p3 ,double *normal );
	extern int  vec_unit(double *vector);
	extern double vec_dotproduct(double *vector1,double *vector2);
	extern double Volume_4p(Point<3> p0,Point<3> p1,Point<3> p2,Point<3> p3);
	extern double area4p(double p0[3] , double p1[3] ,double p2[3] ,double p3[3] );
	extern void vec_2p(double * , double * , double *) ;
	extern int vec_unit(double *) ;
	extern double compMinCos2Faceangle1pTo4p(double p[3],double p4t[4][3]);
	extern void rotatePoints(double (*pts)[3],int numpt);
	extern bool skipBeInfinite(Edge &ed,Triangulation &sft);


	void  STLSurfaceMesh::Get3VertexPtOfTdsTriFromaIndex(const FaceHandle &tri,int i,Point<3> &pa,Point<3> &pb,Point<3> &pc)
	{

		Get1VertexPtOfTdsTri(tri,i,pa);
		Get1VertexPtOfTdsTri(tri,ccw(i),pb);
		Get1VertexPtOfTdsTri(tri,cw(i),pc);

	}
	void  STLSurfaceMesh::Get1VertexPtOfTdsTri(const FaceHandle &tri,int i,Point<3> &p)
	{

		Point<3> pt=tri->GetVertex(i)->GetPoint();
		copyPtFromVersionStructToArray(pt,p);
	}

	bool STLSurfaceMesh::SetEdgeProperty(Edge2v &e2v,void * ptprop)
	{

		Edge ed;
		if(!ChangeEdgeFormFrom2vToTn(e2v,ed)) return false;
		SetEdgeProperty(ed.first,ed.second,ptprop);
		return true;
	}
	void STLSurfaceMesh::SetEdgeProperty(FaceHandle &tri, int index,void * ptprop)
	{
		tri->Info().eproperty[index]=ptprop;
		FaceHandle trin=tri->GetNeighbor(index);
		int ind=trin->Index(tri);
		trin->Info().eproperty[ind]=ptprop;
	}

	bool STLSurfaceMesh::ChangeEdgeFormFrom2vToTn(Edge2v e2v, Edge &e)
	{

		if(!sft.IsEdge(e2v.va,e2v.vb,e.first,e.second)){
			cout<<e2v.va->GetPoint()<<"and"<<e2v.vb->GetPoint()<<endl;
			cout<<("asdfasd change edgeform");
			if(! sft.IsEdge(e2v.vb,e2v.va,e.first,e.second)){
				cout<<("change edgeformerrag");
				return false;
			}
		}
		return true;
	}

	void STLSurfaceMesh::InitializeAllVertexFlag(void)
	{

		for(Triangulation::AllVerticesIterator vit=sft.AllVerticesBegin(); vit!=sft.AllVerticesEnd(); vit++)
			SetFlagOfTdsVertex(vit,STLSurfaceMesh::initial);
	}

	void  STLSurfaceMesh::SetFlagOfTdsVertex(const VertexHandle &v, int flag)
	{

		(v->Info()).status=flag;
	}

	void STLSurfaceMesh::FindTheNearestEdgeOfTdsTriToPtWithSqdist(Point<3> pt,const FaceHandle &tri,Edge &edge,double &dist)
	{

		double dist0=std::numeric_limits<double>::max();
		for(int i=0; i<3; i++)
			if((dist=SqDistPtToTdsEdge(pt,tri,i))<dist0)
			{
				dist0=dist;
				edge=Edge(tri,i);
			}
			dist=dist0;
	}
	double  STLSurfaceMesh::SqDistPtToTdsEdge(Point<3> pt,const FaceHandle &tri,int i)
	{

		Point<3> pa,pb;
		//	Face_handle tri;
		double ppt[3],ppa[3],ppb[3];
		Get2VertexPtOfTdsEdge(Edge(tri,i),pa,pb);
		for(int i = 0;i < 3; i++)
		{
			ppt[i]=pt[i];
			ppa[i]=pa[i];
			ppb[i]=pb[i];
		}
		double d=DistPointToSegm(ppt,ppa,ppb,0,0);
		return d*d;
	}

	void  STLSurfaceMesh::Get2VertexPtOfTdsEdge(Edge edge,Point<3> &pa,Point<3> &pb)
	{

		FaceHandle tri;
		int i;

		tri=edge.first;
		i=edge.second;
		Get1VertexPtOfTdsTri(tri,(i+1)%3,pa);
		Get1VertexPtOfTdsTri(tri,(i+2)%3,pb);

	}

	void  STLSurfaceMesh::Insert1PtInInteriorOfTdsTri(Point<3> &pt,const FaceHandle &tri,VertexHandle &v)
	{

		v=sft.InsertInFace(tri);
		Point<3> p(pt[0],pt[1],pt[2]);
		v->SetPoint(p);
		MaintainTriPropertyAroundaVertex(v);
	}

	void STLSurfaceMesh::MaintainTriPropertyAroundaVertex(VertexHandle &v)
	{

		FaceCirculator tri=sft.Incident_faces(v),done(tri);
		do
		{
			for(int i=0; i<3; i++)
				tri.fi->Info().eproperty[i]=0;
			int inda=tri.fi->Index(v);
			FaceHandle trin=tri.fi->GetNeighbor(inda);
			int ind=trin->Index(tri.fi);
			tri.fi->Info().eproperty[inda]=trin->Info().eproperty[ind];	
		}while(++tri!=done);

	}

	void  STLSurfaceMesh::Insert1PtOnTdsEdge(double pt[3],Edge edge,VertexHandle &v)
	{

		bool epropflag=false;
		if(IsTdsEdgeOnBoundary(edge)) epropflag=true;
		Edge2v e2v;
		ChangeEdgeFormFromTnTo2v(edge,e2v);
		v=sft.InsertInEdge(edge.first,edge.second);
		v->SetPoint(Point<3>(pt[0],pt[1],pt[2]));
		v->Degree();
		MaintainTriPropertyAroundaVertex(v);
		if(epropflag){
			SetEdgeProperty(Edge2v(v,e2v.va),ptonbd);	
			SetEdgeProperty(Edge2v(v,e2v.vb),ptonbd);
		}
	}

	bool  STLSurfaceMesh::IsTdsEdgeOnBoundary(Edge edge)
	{

		//	Face_handle ta=edge.first->neighbor(edge.second);
		//	if(edge.first->has_vertex(sft.infinite_vertex())||
		//		ta->has_vertex(sft.infinite_vertex()))
		//	Edge2v e2v;
		//	changeEdgeFormFromTnTo2v(edge,e2v);

		if(edge.first->Info().eproperty[edge.second]==ptonbd)
			return true;
		else 
			return false;

	}

	void STLSurfaceMesh::ChangeEdgeFormFromTnTo2v(Edge e, Edge2v &e2v)
	{

		Get2TdsVertexOfEdge(e,e2v.va,e2v.vb);
	}

	void  STLSurfaceMesh::Get2TdsVertexOfEdge(Edge edge,VertexHandle &va,VertexHandle &vb)
	{

		FaceHandle tri=edge.first;
		int i=edge.second;
		va=tri->GetVertex(ccw(i));
		vb=tri->GetVertex(cw(i));
	}

	bool STLSurfaceMesh::CanTdsEdgeSwap(Edge edge,bool geomapp)
	{
		FaceHandle tri=edge.first;
		int i=edge.second;
		int indn=tri->GetNeighbor(i)->Index(tri);
		//cout<<tri->GetVertex((i+0)%3)->GetPoint()<<endl;
		//cout<<tri->GetVertex((i+1)%3)->GetPoint()<<endl;
		//cout<<tri->GetVertex((i+2)%3)->GetPoint()<<endl;
		//int nnn=	tri->GetVertex((i+1)%3)->Degree();
		//int mmm=tri->GetVertex((i+2)%3)->Degree();
		//sft.is_edge(tri->GetVertex(i),tri->GetNeighbor(i)->GetVertex(indn));
		if(tri->GetVertex((i+1)%3)->Degree()<=3||
			tri->GetVertex((i+2)%3)->Degree()<=3||
			sft.IsEdge(tri->GetVertex(i),tri->GetNeighbor(i)->GetVertex(indn)))
			return false;
		//if (sft.is_infinite(tri) || sft.is_infinite(tri->neighbor(i))) //后加入，判断要swap的边是否在边界上
		//	return false;
		if(skipBeInfinite(tri,sft)||skipBeInfinite(tri->GetNeighbor(i),sft)) 
			return false;//jf_error("cantdsedgeswap");
		if(IsTdsEdgeOnBoundary(edge)) return false;

		double pa[3],pb[3],pc[3],pd[3],nabc[3],ndcb[3],nabd[3],ncad[3];
		Point<3>ppa,ppb,ppc,ppd;
		Get3VertexPtOfTdsTriFromaIndex(tri,i,ppa,ppb,ppc);
		GetSpecifiedNeighbVertexPtOfTdsTri(tri,i,ppd);
		for(int i = 0;i < 3;i++)
		{
			pa[i]=ppa[i];
			pb[i]=ppb[i];
			pc[i]=ppc[i];
			pd[i]=ppd[i];
		}
		normal_3p(pa,pb,pc,nabc); //后加入，判断swap前的状态
		vec_unit(nabc);
		normal_3p(pd,pc,pb,ndcb);
		vec_unit(ndcb);

		normal_3p(pa,pb,pd,nabd);
		if(vec_unit(nabd)==0) i=i;
		//jf_error("maybe 1wrong",1);
		normal_3p(pc,pa,pd,ncad);
		if(vec_unit(ncad)==0) i=i;
		//jf_error("maybe 2wrong",1);
		double tol=0.02866;
		//	if(geomapp) tol=0.999;
		//0.00015*mparam.maxh的选取有什么道理
		if(geomapp) return fabs(Volume_4p(ppa,ppb,ppc,ppd)/area4p(pa,pb,pc,pd))<0.00015*mparam.minh&&vec_dotproduct(nabd,ncad)>tol;
		//vec_dotp(nabc,ndcb)>tol;

		if(vec_dotproduct(nabd,ncad)>tol)//||vec_dotp(ndcb,nabc)<vec_dotp(nabd,ncad))
			return true; //false;
		else
			return false;//true;
	}

	void  STLSurfaceMesh::GetSpecifiedNeighbVertexPtOfTdsTri(FaceHandle tri,int i,Point<3> &p)
	{

		FaceHandle ta=tri->GetNeighbor(i);
		int index=ta->Index(tri);
		Get1VertexPtOfTdsTri(ta,index,p);
	}

	bool  STLSurfaceMesh::IsImprovedAfterEdgeSwap(Edge edge)
	{

		FaceHandle tri=edge.first;
		int i=edge.second;
		double cs0,cs1;
		Point<3> p4t[4];
		int k;
		Get4VertexPtOf2NeighbTri(tri,i,p4t);
		SelectTheDiagonalOf4Pt(p4t,k,cs0,cs1);
		if(k==0) return false;
		else
			return true;
	}
	void  STLSurfaceMesh::Get4VertexPtOf2NeighbTri(const FaceHandle &tri,int i,Point<3> p4t[4])
	{

		Get3VertexPtOfTdsTriFromaIndex(tri,i,p4t[0],p4t[1],p4t[2]);
		GetSpecifiedNeighbVertexPtOfTdsTri(tri,i,p4t[3]);
	}

	void  STLSurfaceMesh::SelectTheDiagonalOf4Pt(Point<3> p4t[4],int &k,double &cos0,double &cos1)
	{

		double cosa,cosb;//,cos0,cos1;
		//	int k;
		cosa=ComputeCosValueOfAngleInaTriangle(p4t[0],p4t[1],p4t[2]);
		cosb=ComputeCosValueOfAngleInaTriangle(p4t[2],p4t[3],p4t[1]);
		cos0=max(cosa,cosb);
		cosa=ComputeCosValueOfAngleInaTriangle(p4t[0],p4t[1],p4t[3]);
		cosb=ComputeCosValueOfAngleInaTriangle(p4t[3],p4t[2],p4t[0]);
		cos1=max(cosa,cosb);
		if(cos0>cos1)
			k=1;
		else 
			k=0;
	}

	double  STLSurfaceMesh::ComputeCosValueOfAngleInaTriangle(Point<3> pa,Point<3> pb,Point<3> pc)
	{

		double pab[3],pac[3],pbc[3],dabac,dbabc,dcacb;
		double ppa[3],ppb[3],ppc[3];
		for(int i = 0; i < 3 ;i++)
		{
			ppa[i]=pa[i];
			ppb[i]=pb[i];
			ppc[i]=pc[i];
		}
		vec_2p(ppa,ppb,pab);
		vec_unit(pab);
		vec_2p(ppa,ppc,pac);
		vec_unit(pac);
		vec_2p(ppb,ppc,pbc);
		vec_unit(pbc);
		dabac=vec_dotproduct(pab,pac);
		dbabc=-vec_dotproduct(pab,pbc); //
		dcacb=vec_dotproduct(pac,pbc);
		return max(dabac,max(dbabc,dcacb));
	}

	void STLSurfaceMesh::SwapEdge(FaceHandle tri1,int index)
	{

		FaceHandle tri2=tri1->GetNeighbor(index);
		sft.Flip(tri1,index);
		SetEdgeProperty(tri1,tri1->Index(tri2),0);
		MaintainTriProperty(tri1);
		MaintainTriProperty(tri2);	

	}

	void STLSurfaceMesh::MaintainTriProperty(FaceHandle &tri)
	{

		for(int i=0; i<3; i++)
		{
			FaceHandle trin=tri->GetNeighbor(i);
			int ind=trin->Index(tri);
			tri->Info().eproperty[i]=trin->Info().eproperty[ind];
		}
	}

	bool STLSurfaceMesh::IsVertexRemain(const VertexHandle &v)
	{

		return (v->Info()).status!=STLSurfaceMesh::initial;
	}

	bool STLSurfaceMesh::CanVertexDel(VertexHandle v)
	{  //后加入，判断该点是不是边界上能删的点

		return true;

	}

	bool STLSurfaceMesh::DeleteOneOldNode( const VertexHandle &v, std::queue<VertexHandle> *pqvt)
	{

		void *bkaddress= v->Info().bkvaddress; //not maintain qvt after deletion? since needn't to search?

		if(pqvt)
		{
			VertexCirculator vs=sft.Incident_vertices(v),ve=vs;
			do
			{
				if(skipBeInfinite(VertexHandle(ve.v),sft)) continue;	
				pqvt->push(ve.v);
			}while(++ve!=vs);
		}
		if(v->Degree()>4) ReduceDegreeBySwap(v);
		//sft.OutputSTL("test.stl");
		//cout<<v->Degree()<<endl;
		if(v->Degree()>4)
		{  
			return false;
		}
		int intsegcount;
		Edge2v e2v;
		CompPropertyOfVertex(v,intsegcount,e2v);
		int tdeg=v->Degree();
		if(v->Degree()<=2||(intsegcount!=0&&intsegcount!=2))
		{
			cout<<v->GetPoint()[0]<<" "<<v->GetPoint()[1]<<" "<<v->GetPoint()[2]<<" "<<intsegcount<<endl;
			return false;
		}
		//	jf_error("delete1oldnode",1);
		FaceHandle tria,trib;
		if(v->Degree()==3)
		{
			tria=v->GetFace(); //should be check carefully
			sft.RemoveDegree_3(v,tria);
			//sft.OutputSTL("test.stl");
			MaintainTriProperty(tria);
		}else
		{
			//暂时注释掉，看看效果
		
			VertexHandle vm=e2v.va;
			
			if((intsegcount==0||sft.IsEdge(e2v.va,e2v.vb))&& !ChooseAMergeTdsVertexFrom4Incidents(v,vm))
			{  
				return false;
			}
			
			if(!Can2TdsVertsJoin(v,vm))
			{   
				return false;
			}
			Merge2TdsVertexOf1TdsEdge(v,vm,tria,trib);
			SetEdgeProperty(tria,tria->Index(trib),0);
			MaintainTriProperty(tria);
			MaintainTriProperty(trib);
		}
		if(intsegcount==2) SetEdgeProperty(e2v,ptonbd);
		//	delete v->info().bkvaddress; v->info().bkvaddress=0;
		return true;
	}

	void  STLSurfaceMesh::ReduceDegreeBySwap(const VertexHandle &v)
	{

		FaceCirculator tri1=sft.Incident_faces(v),tri2(tri1),trin;
		tri2++;
		int swflag=v->Degree();
		do{
			trin=tri2, trin++; swflag--;	
			if( !skipBeInfinite(FaceHandle(tri1.fi),sft)&&!skipBeInfinite(FaceHandle(tri2.fi),sft)&&IsTriPairCorrect(tri1.fi,tri2.fi)&&CanTdsEdgeSwap(Edge(tri1.fi,tri1.fi->Index(tri2.fi))) )
			{
				SwapEdge(tri1.fi,tri1.fi->Index(tri2.fi));
				swflag=v->Degree();
				if(swflag==3) return;
			}
			tri2=tri1=trin;
			tri1--;
		}while(swflag>0);
	}

	bool  STLSurfaceMesh::IsTriPairCorrect(const FaceHandle &tri1,const FaceHandle &tri2)
	{

		if(!tri1->IsNeighbor(tri2))
			return false;
		else 
			return true;
	}

	void STLSurfaceMesh::CompPropertyOfVertex(const VertexHandle &v,int &intsegcount,Edge2v &e2v)
	{

		VertexCirculator vs=sft.Incident_vertices(v),ve=vs;
		intsegcount=0;
		VertexHandle edbd[2];
		do
		{
			Edge ed;
			ChangeEdgeFormFrom2vToTn(Edge2v(v,ve.v),ed);
			if(IsTdsEdgeOnBoundary(ed))
			{
				if(intsegcount<2)
					edbd[intsegcount++]=ve.v;
			}
		}while(++ve!=vs);
		e2v.va=edbd[0]; e2v.vb=edbd[1];
	}

	bool STLSurfaceMesh::ChooseAMergeTdsVertexFrom4Incidents(const VertexHandle &v,VertexHandle &vm)
	{

		FaceHandle tri;
		//	int i;
		//	if(sft.is_edge(v,sft.infinite_vertex(),tri,i)){
		//		vm=tri->vertex(i);
		//		return true;
		//	}
		//	int dg4v[4],id;
		Point<3> p4t[4];
		VertexCirculator vc=sft.Incident_vertices(v);
		for(int i=0; i<4; i++)
		{
			//		if(vc->degree()<=3) return false;
			//if(skipBeInfinite(vc,sft)) jf_error("err chooseamergetds");
			copyPtFromVersionStructToArray((vc++).v->GetPoint(),p4t[i]);
		}
		//	if(!(dg4v[0]>3&&dg4v[2]>3||dg4v[1]>3&&dg4v[3]>3))
		//		return false;
		Point<3> p;
		copyPtFromVersionStructToArray(v->GetPoint(),p);
		double pp[3],pp4t[4][3];
		for(int i = 0; i < 3; i++)
		{
			pp[i]=p[i];
			pp4t[0][i]=p4t[0][i];
			pp4t[1][i]=p4t[1][i];
			pp4t[2][i]=p4t[2][i];
			pp4t[3][i]=p4t[3][i];
		}
		double r1=compMinCos2Faceangle1pTo4p(pp,pp4t);
		rotatePoints(pp4t,4);
		double r2=compMinCos2Faceangle1pTo4p(pp,pp4t);
		if(r1==-1.&&r2==-1.)
			return false;
		//jf_error("err chooseamergetdsvert",1);
		else if(r1>r2)
			vm=vc.v;
		else
			vm=(++vc).v;
		return true;
	}

	bool  STLSurfaceMesh::Can2TdsVertsJoin(const VertexHandle &va,const VertexHandle &vb)
	{

		if(skipBeInfinite(VertexHandle(va),sft)||skipBeInfinite(VertexHandle(vb),sft))
		{
			cout<<("errcand2tdsvertsjoin");return false;
		}
		if(!sft.IsEdge(va,vb))
			return false;
		VertexCirculator vs=sft.Incident_vertices(va),ve=vs;
		int count=0;
		do
		{
			if(sft.IsEdge(ve.v,vb))
				count++;
			if(count>=3) 
				return false;
		}while(++ve!=vs);
		return true;
	}

	void  STLSurfaceMesh:: Merge2TdsVertexOf1TdsEdge(const VertexHandle &v,const VertexHandle &vm,FaceHandle &ta,FaceHandle &tb)
	{

		FaceHandle fr;
		int i;
		Point<3> p=vm->GetPoint();   //错误，v应为vm,已修改
		if(!sft.IsEdge(v,vm,fr,i))
			cout<<("err merge2tdsvertex")<<endl;
		int n=fr->Index(v);
		FaceHandle fn=fr->GetNeighbor(n);
		n=fn->Index(fr);
		VertexHandle vn=sft.JoinVertices(fr,i,vm);
		vn->SetPoint(p);
		ta=fn->GetNeighbor(n);
		int ind=3-ta->Index(fn)-ta->Index(vn);
		tb=ta->GetNeighbor(ind);
	}

	void STLSurfaceMesh:: LocalImpAfter1VtDel(std::list<Edge2v> &queueofedge, std::queue<VertexHandle> &quevt)
	{
		Edge edgen;
		while(!quevt.empty())
		{
			VertexHandle v=quevt.front();
			quevt.pop();
			VertexCirculator vs=sft.Incident_vertices(v);
			VertexCirculator ve=vs;
			do
			{
				ChangeEdgeFormFrom2vToTn(Edge2v(v,ve.v),edgen);
				if(skipBeInfinite(edgen,sft)) continue;	
				PushFreeEdgeInque(edgen.first,edgen.second,queueofedge);
			}while(++ve!=vs);
		}
		ImpSurfMeshFromEdgesInQue(queueofedge);
	}

	void STLSurfaceMesh:: PushFreeEdgeInque(FaceHandle tri,int ind,std::list<Edge2v> &queueofedge)
	{

		Edge2v e2v;

		if(tri->Info().eproperty[ind]!=0)
			return;
		ChangeEdgeFormFromTnTo2v(Edge(tri,ind),e2v);
		queueofedge.push_back(e2v); 
		SetEdgeProperty(tri,ind,ptinque);
	}

	void	STLSurfaceMesh:: ImpSurfMeshFromEdgesInQue(std::list<Edge2v> &queueofedge)
	{

		Edge2v e2v;
		Edge edgen;

		while(!queueofedge.empty())
		{
			e2v=*queueofedge.begin();		
			queueofedge.pop_front();
			ChangeEdgeFormFrom2vToTn(e2v,edgen);
			SetEdgeProperty(edgen.first,edgen.second,0);
			if(CanTdsEdgeSwap(edgen)&&
				IsImprovedAfterEdgeSwap(edgen))
			{
				FaceHandle fsecond=edgen.first->GetNeighbor(edgen.second);
				SwapEdge(edgen.first,edgen.second);
				UpdateQueueAndSwapStatusOfEdgesInQueue(queueofedge,edgen.first,fsecond);
			}

			// it is supposed that 2 Face_handles are still valid after flip.
			// otherwise more works should be done to achieve the 2 triangles.
		}
	}

	void STLSurfaceMesh::UpdateQueueAndSwapStatusOfEdgesInQueue(list<Edge2v> &queueofedge,FaceHandle &tria, FaceHandle &trib)
	{

		VertexHandle va,vb;

		for(int i=0; i<3; i++)
		{
			if(tria->GetNeighbor(i)!=trib)
				PushFreeEdgeInque(tria,i,queueofedge);
		}
		for(int i=0; i<3; i++)
		{
			if(trib->GetNeighbor(i)!=tria)
				PushFreeEdgeInque(trib,i,queueofedge);
		}			
	}

	void STLSurfaceMesh:: ImpSurfaceMesh(	std::list<Edge2v> &queueofedge)
	{

		/*Vertex_handle vc;
		int dcount[10];
		for(int i=0; i<10; i++) dcount[i]=0;
		for(Triangulation::All_vertices_iterator vc=sft.all_vertices_begin(); vc!=sft.all_vertices_end(); vc++)
		dcount[vc->degree()-3]++;
		for(int i=0; i<10; i++)
		cout<<dcount[i]<<" ";	
		cout<<endl;*/
		//	for(;;){
		Edge edgen;
		for(Triangulation::AllFacesIterator fit=sft.AllFacesBegin(); fit!=sft.AllFacesEnd(); fit++)
		{
			for(int i = 0; i < 3; i++)
				fit->SetEdgeStatus(0,i);
		}
		for(Triangulation::AllFacesIterator fit=sft.AllFacesBegin(); fit!=sft.AllFacesEnd(); fit++)
		{
			for(int i = 0; i < 3 ;i++)
			{
				if(fit->EdgeStatus(i)==1) continue;
				fit->SetEdgeStatus(1,i);
				fit->GetNeighbor(i)->SetEdgeStatus(1,fit->GetNeighbor(i)->Index(fit));
				Edge eit(fit,i);
				if(skipBeInfinite(eit,sft)) continue;	
				edgen=eit;
				Edge2v e2v;
				ChangeEdgeFormFromTnTo2v(edgen,e2v);
				//		if(e2v.va->degree()<8&&e2v.vb->degree()<8) continue;
				PushFreeEdgeInque(edgen.first,edgen.second,queueofedge);
			}

		}
		ImpSurfMeshFromEdgesInQue(queueofedge);

		/*for(int i=0; i<10; i++) dcount[i]=0;
		for(Triangulation::All_vertices_iterator vc=sft.all_vertices_begin(); vc!=sft.all_vertices_end(); vc++)
		dcount[vc->degree()-3]++;
		for(int i=0; i<10; i++)
		cout<<dcount[i]<<" ";	*/
		//	}
	}


	bool skipBeInfinite(VertexHandle &vt,Triangulation &sft)
	{

		return  sft.HaveInfinite()&&sft.IsInfinite(vt);
	}
	bool skipBeInfinite(FaceHandle &fa,Triangulation &sft)
	{

		return sft.HaveInfinite()&&sft.IsInfinite(fa);
	}
	bool skipBeInfinite(Edge &ed,Triangulation &sft)
	{

		return  sft.HaveInfinite()&&(sft.IsInfinite(ed.first)||sft.IsInfinite(ed.first->GetNeighbor(ed.second)));
	}

	void copyPtFromVersionStructToArray(Point<3> p, Point<3> &pt)
	{
		pt[0]=p[0];
		pt[1]=p[1];
		pt[2]=p[2];
	}

}
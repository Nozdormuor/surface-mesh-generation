#include"externalcommon.h"
#include"internalcommon.h"
namespace meshwork
{
	extern  void jf_error(char *);
	extern void copyPtFromVersionStructToArray(Point<3> p, Point<3> &pt);
	extern bool skipBeInfinite(VertexHandle &vt,Triangulation &sft);
	extern bool skipBeInfinite(FaceHandle &fa,Triangulation &sft);
	extern bool IsTriangleAndBoxOverlap(Point<3> p1 ,Point<3> p2 ,Point<3> p3 ,const double bd[6],double eps );
	extern double sqDistPointToTri(double p[3],double p0[3],double p1[3],double p2[3]);
	extern int is_point_same( Point<3>p1 , Point<3> p2 );
	using namespace std;
	const double TriangulationInKdtree::epsoverlap=0.000001;

	TriangulationInKdtree::TriangulationInKdtree(STLSurfaceMesh *pM)
	{
		Box<3> box;
		Triangulation::AllVerticesIterator vi, vn;
		box[0]=box[1]=box[2]=std::numeric_limits<double>::max();
		box[3]=box[4]=box[5]=-box[0];
		int count=0;

		for(vi=pM->sft.AllVerticesBegin();vi!=pM->sft.AllVerticesEnd();vi++)
		{
			if(skipBeInfinite(vi,pM->sft)) continue;
			Triangulation::VertexHandle vh=(Triangulation::VertexHandle)vi;
			box[0]=min(box[0],vh->GetPoint()[0]);
			box[1]=min(box[1],vh->GetPoint()[1]);
			box[2]=min(box[2],vh->GetPoint()[2]);
			box[3]=max(box[3],vh->GetPoint()[0]);
			box[4]=max(box[4],vh->GetPoint()[1]);
			box[5]=max(box[5],vh->GetPoint()[2]);
			count++;
		}
		double dinc=max(box[3]-box[0],box[4]-box[1]);
		dinc=max(dinc,box[5]-box[2])*0.0001;
		for(int i=0; i<3; i++)
		{
			box[i]-=dinc;
			box[i+3]+=dinc;
		}
		fvtree=new Kodtree(box,pofvforcoordnodes3,2,epsoverlap);
		fvtree->SetFuncExinfoShouldbeInCell(ifexinfoshouldbeincell);
		fvtree->SetFuncExinfoOverlapBox(ifexinfooverlapbox);
		pMesh=pM;
		InsertSurfremesh();
	}

	TriangulationInKdtree::~TriangulationInKdtree()
	{
		for(Triangulation::AllVerticesIterator vi=pMesh->sft.AllVerticesBegin();vi!=pMesh->sft.AllVerticesEnd();vi++)
		{
			//		delete vi->info().bkvaddress;
			//		vi->info().bkvaddress=0;
		}
		for(Triangulation::AllFacesIterator fi=pMesh->sft.AllFacesBegin();fi!=pMesh->sft.AllFacesEnd();fi++)
		{
			//		delete fi->info().bkaddress;
			//		fi->info().bkaddress=0;
		}
		delete fvtree;
	}

	void TriangulationInKdtree::InsertSurfremesh()
	{
		Triangulation *tri=&(pMesh->sft);
		for(Triangulation::AllVerticesIterator vi=tri->AllVerticesBegin();vi!=tri->AllVerticesEnd();vi++)
		{
			if(skipBeInfinite(vi,pMesh->sft)) continue;
			Triangulation::VertexHandle vh=vi;
			//		delete vh->info().bkvaddress;
			vh->Info().bkvaddress=InsertVertexHandle(vh);
			//		vh->info().haveR=false;
		}

		for(Triangulation::AllFacesIterator fi=tri->AllFacesBegin();fi!=tri->AllFacesEnd();fi++)
		{
			if(skipBeInfinite(fi,pMesh->sft)) continue;
			FaceHandle fh=fi;
			//		delete fh->info().bkaddress;
			fh->Info().bkaddress=InsertFaceHandle(fh);
		}
	}

	bool   TriangulationInKdtree::ifexinfoshouldbeincell(void *info,int infotype,CellNode *cnode)
	{
		if(infotype==1)
		{
			return true; // 2015-7-19 otherwise some cells will lose parts of exinfos after split. 
			Triangulation::FaceHandle pFH=*((Triangulation::FaceHandle *)info);
			//Triangulation::Vertex_handle v1=pFH->vertex(0);
			//Triangulation::Vertex_handle v2=pFH->vertex(1);
			//Triangulation::Vertex_handle v3=pFH->vertex(2);

			//double pt[3][3],p[3];

			//pofvforcoordnodes3(pt[0],&v1);
			//pofvforcoordnodes3(pt[1],&v2);
			//pofvforcoordnodes3(pt[2],&v3);

			for(int i=0;i<cnode->numvert;i++)
			{
				Triangulation::VertexHandle vh=*((Triangulation::VertexHandle*)(cnode->vert[i]->vt));
				//pofvforcoordnodes3(p,&vh);
				for(int j=0;j<3;j++)//{
					//	double temp=pointToPointDist(p,pt[j]);
					//	if(temp<0.00001){
					if(vh==pFH->GetVertex(j))
						return false;
				//	}
				//}
			}
		}
		return true;
	}

	bool   TriangulationInKdtree::ifexinfooverlapbox(void *info,int infotype,const Box<3> &bd,double eps)
	{
		double bbd[6];
		bbd[0]=bd[0];bbd[1]=bd[1];bbd[2]=bd[2];bbd[3]=bd[3];bbd[4]=bd[4];bbd[5]=bd[5];
		Point<3> pt1,pt2,pt3;
		Triangulation::VertexHandle vh1,vh2,vh3;

		if(infotype==1)
		{
			Triangulation::FaceHandle fh=*((Triangulation::FaceHandle *)info);

			vh1=fh->GetVertex(0);
			vh2=fh->GetVertex(1);
			vh3=fh->GetVertex(2);

			pofvforcoordnodes3(pt1,&vh1);
			pofvforcoordnodes3(pt2,&vh2);
			pofvforcoordnodes3(pt3,&vh3);

			return IsTriangleAndBoxOverlap(pt1,pt2,pt3,bbd,eps);
		}
		return false;
	}

	void  TriangulationInKdtree::pofvforcoordnodes3(Point<3> &p,void *pv)
	{

		Triangulation::VertexHandle vh=*((Triangulation::VertexHandle *)pv);
		Point<3> pt=vh->GetPoint();
		p[0]=pt[0];
		p[1]=pt[1];
		p[2]=pt[2];
	}
	bool TriangulationInKdtree::GetFacehandleFromPoint(Point<3> pp,Triangulation::FaceHandle &fh,double *pdist)
	{

		Point<3> p0,p1,p2;
		double distemp,dist=std::numeric_limits<double>::max();
		CellNode *pcell;
		if(!( pcell=fvtree->FindaLeafCellContainingPoint(fvtree->GetRoot(),pp)))
		{
			fvtree->FindaLeafCellContainingPoint(fvtree->GetRoot(),pp);
			jf_error("getfacehandle");
		}
		list<void *> linfo;
		fvtree->CollectExinfoWithCell(pcell,1,linfo);
		if(linfo.size()==0)
		{
			Box<3> bd;double dcent[3];
			for(int i=0; i<6; i++) bd[i]=pcell->bound[i];
			for(int i=0; i<3; i++) dcent[i]=(bd[i+3]+bd[i])/2.-pp[i];
			for(int i=0; i<6; i++) bd[i]=bd[i]-dcent[i%3];
			//    double dlp=0.1;
			//    for(int i=0; i<3; i++){
			//	   bd[i]=p[i]-dlp;
			//	   bd[3+i]=p[i]+dlp;
			//    }
			fvtree->CollectExinfoWithBox(bd,1,linfo);
			cout<<"bigcell "<<linfo.size()<<" ";
			if(linfo.size()==0) jf_error("getfacenhandleerrorbound");
		}
		for( list<void *>::iterator pos=linfo.begin();pos!=linfo.end();pos++)
		{
			FaceHandle fc=*((FaceHandle *)*pos);
			pMesh->Get3VertexPtOfTdsTriFromaIndex(fc,0,p0,p1,p2);
			double ppp[3],pp0[3],pp1[3],pp2[3];
			for(int i = 0; i< 3;i++)
			{
				ppp[i]=pp[i];
				pp0[i]=p0[i];
				pp1[i]=p1[i];
				pp2[i]=p2[i];
			}
			if((distemp=sqDistPointToTri(ppp,pp0,pp1,pp2))<dist)
			{
				dist=distemp;
				fh=fc;
			}
		}
		/*
		vector<void *> vecvert;
		fvtree->collectVertsWithCell(pcell,vecvert);
		for( vector<void *>::iterator pos=vecvert.begin();pos!=vecvert.end();pos++){
		Vertex_handle vh=*((Vertex_handle *)*pos);
		Triangulation::Face_circulator fc=pMesh->sft.incident_faces(vh),end=fc;
		do{
		if(skipBeInfinite(FaceHandle(fc),pMesh->sft)) continue;	
		pMesh->get3VertexPtOfTdsTriFromaIndex(fc,0,p0,p1,p2);
		if((distemp=sqDistPointToTri(p,p0,p1,p2))<dist){
		dist=distemp;
		fh=fc;
		}
		} while(++fc!=end);
		}
		*/
		return true;
	}
	bool TriangulationInKdtree:: GetVertexhandleFromPoint(Point<3> p,Triangulation::VertexHandle &vh)
	{

		Box<3> bd;
		double dlp=0.0000001;
		for(int i=0; i<3; i++)
		{
			bd[i]=p[i]-dlp;
			bd[3+i]=p[i]+dlp;
		}
		list<void *> lvert;
		fvtree->CollectVertsWithBox(bd,lvert);
		for( list<void *>::iterator pos=lvert.begin();pos!=lvert.end();pos++)
		{
			VertexHandle vc=*((VertexHandle *)*pos);
			Point<3> pvc;
			copyPtFromVersionStructToArray(vc->GetPoint(),pvc);
			//if(Distance3D(pvc,p)<0.000005)
			if(is_point_same(pvc,p))
			{
				vh=vc;
				return true;
			}
		}
		return false;
	}

}
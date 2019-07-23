#include"externalcommon.h"
#include"internalcommon.h"
namespace meshwork
{
	//找到包含该节点的第一个父节点
	CellNode * Kodtree:: FindTheNearestAncestorContainingPoint(CellNode *pcell0,Point<3> pcha)
	{
		CellNode *pcell=pcell0;
		for(;;)
		{
			if(pcell==0) 
				return 0;
			if(IsBoxContainPoint(pcha,pcell->bound,root->bound)) 
				return pcell;//need or not to follow the convention?
			else pcell=pcell->parent;
		}
	}
	//找到包含该点的叶子节点
	CellNode *  Kodtree::FindaLeafCellContainingPoint(CellNode *pcell,Point<3> p)
	{
		CellNode *rtpcell;

		if(!pcell||!IsBoxContainPoint(p,pcell->bound,root->bound))//need or not to follow the convention?
			return 0;
		if(pcell->IsLeaf())
			return pcell;
		for( int i=0; i<2; i++)
			if((rtpcell=FindaLeafCellContainingPoint(pcell->child[i],p))!=0)
				return rtpcell;
		jf_error("err findaleafcellcontainp");
	}
	//判断两个单元是否相交
	bool  Kodtree::If2CellNeighb(CellNode *pcell0, CellNode *pcell1)
	{

		if(!pcell0||!pcell1) 
			jf_error("err is2cellneigh");
		if(IsTwoBoxNeighber(pcell0->bound,pcell1->bound))
			return true;
		else
			return false;
	}


	Kodtree::Kodtree(const Box<3> &bd,Funcpointofvert pofvin,int capacity,double epsi)
	{

		double lcube=max(bd[3]-bd[0],max(bd[4]-bd[1],bd[5]-bd[2]));
		epscell=numeric_limits<double>::epsilon()*(1+lcube*10.);
		root=new CellNode(bd);
		pofv=pofvin;
		cellcapacity=capacity;
		epsoverlap=epsi;
	}

	Kodtree::Kodtree(void **vert, int numvert,Funcpointofvert pofvin,int capacity,double epsi)
	{

		Box<3> bd;
		BoxOfVertices(vert,numvert,bd,pofvin);
		double lcube=max(bd[3]-bd[0],max(bd[4]-bd[1],bd[5]-bd[2]));
		epscell=numeric_limits<double>::epsilon()*(1+lcube*10.);
		root=new CellNode(bd);
		pofv=pofvin;
		cellcapacity=capacity;
		epsoverlap=epsi;
		for( int i=0; i<numvert; i++)
			InsertVert(vert[i]);
	}

	Kodtree::Kodtree(void **vert, int numvert,const Box<3> &bd,Funcpointofvert pofvin,int capacity,double epsi)
	{

		//	Box bd;
		//	boxOfPoints(vert,numvert,bd);
		double lcube=max(bd[3]-bd[0],max(bd[4]-bd[1],bd[5]-bd[2]));
		epscell=numeric_limits<double>::epsilon()*(1+lcube*10.);
		root=new CellNode(bd);
		pofv=pofvin;
		cellcapacity=capacity;
		epsoverlap=epsi;
		for( int i=0; i<numvert; i++)
			InsertVert(vert[i]);
	}



	Kodtree::~Kodtree(){

		FreeSubTree(root);
	}


	//插入顶点
	void Kodtree::InsertWpVertInSubTree(const Point<3> &p, WpVert *v, CellNode *cnode)
	{

		if(!cnode)
			jf_error("err insvinst");
		if(!IsPointOverlapWithBox(p,cnode->bound,root->bound ,epsoverlap ))
			return ;
		if(!cnode->IsLeaf())
		{
			for(int i=0; i<2; i++)
				InsertWpVertInSubTree(p,v,cnode->child[i]);
			cnode->numvert++;
			return ;
		}
		if(!cnode->vert)
		{
			cnode->vert =(PtWpVert *) new PtWpVert[Kodtree::cellcapacity];
			cnode->nodecapacity=Kodtree::cellcapacity;
		}else if(cnode->numvert==cnode->nodecapacity&&!CanNodeSplit(cnode))
		{
			cnode->vert=renew(cnode->vert,cnode->nodecapacity,2*cnode->nodecapacity);
			cnode->nodecapacity*=2;
		}
		if(cnode->numvert<cnode->nodecapacity)
		{
			cnode->vert[cnode->numvert++]=v;
			v->rcount ++;
		}else
		{
			SplitNode(cnode);
			for(int i=0; i<2; i++)
				InsertWpVertInSubTree(p,v,cnode->child[i]);
			cnode->numvert++;
			return;
		}
	}

	//判断顶点是否在树中有记录
	bool Kodtree::IsVertRecordInSubTree(const Point<3> &p,void *v, CellNode *cnode)
	{

		if(!cnode)
			jf_error("err insvinst");
		if(cnode->numvert<=0||!IsPointOverlapWithBox(p,cnode->bound,root->bound ,epsoverlap ))
			return false;
		if(!cnode->IsLeaf())
		{
			for(int i=0; i<2; i++)
				if(IsVertRecordInSubTree(p,v,cnode->child[i])) return true;
			return false;
		}
		if(!cnode->vert)
			jf_error("err insvinst");
		for(int i=0; i<cnode->numvert; i++)
			if(cnode->vert[i]->vt==v) return true;
		return false;
	}
	//计算顶点的个数
	int comWpVertNum(CellNode *cnode, CellNode *cnsib)
	{

		int num=0;
		for(int i=0; i<cnsib->numvert; i++)
		{
			for(int j=0; j<cnode->numvert; j++)
				if(cnsib->vert[i]==cnode->vert[j])
				{
					num++; 
					break;
				}
		}
		return num;
	}
	//从树中删除顶点
	void Kodtree::DeleteVertInSubTree(const Point<3> &p,void *v, CellNode *cnode)
	{

		if(!cnode)
			jf_error("err insvinst");
		if(!IsPointOverlapWithBox(p,cnode->bound,root->bound ,epsoverlap ))
			return ;
		cnode->numvert--;
		if(!cnode->IsLeaf())
		{
			for(int i=0; i<2; i++)
				DeleteVertInSubTree(p,v,cnode->child[i]);
			return ;
		}
		if(!cnode->vert)
			jf_error("err deletevertinsubtree");
		int i;
		for(i=0; i<cnode->numvert; i++)
			if(cnode->vert[i]->vt==v)
				break;
		if(--(cnode->vert[i]->rcount)<=0) delete cnode->vert[i];
		if(i!=cnode->numvert)
			cnode->vert[i]=cnode->vert[cnode->numvert];
		if(cnode->numvert==0)
		{
			delete cnode->vert;
			cnode->vert=0;
		}

	}
	//在树中插入额外的信息
	void Kodtree::InsertWpInfoInSubTree(WpInfo *pwinfo, CellNode *cnode)
	{

		if(!cnode)
			jf_error("err insvinst");
		if(!ifExinfoOverlapBox(pwinfo->info,pwinfo->infotype ,cnode->bound,epsoverlap ))
			return ;
		if(!cnode->IsLeaf())
		{
			for(int i=0; i<2; i++)
				InsertWpInfoInSubTree(pwinfo,cnode->child[i]);
			return ;
		}
		if(!ifExinfoShouldbeInCell(pwinfo->info,pwinfo->infotype ,cnode ))
			return ;
		if(!cnode->lpwpinfo)
		{
			cnode->lpwpinfo=new std::list<WpInfo *>;
		}
		cnode->lpwpinfo->push_back(pwinfo);
		pwinfo->rcount++;
	}
	//删除树中额外的信息
	void Kodtree::DeleteExinfoInSubTree(void *info,int infotype, CellNode *cnode)
	{
		if(!cnode)
			jf_error("err insvinst");
		if(!ifExinfoOverlapBox(info,infotype ,cnode->bound,epsoverlap ))
			return ;
		if(!cnode->IsLeaf())
		{
			for(int i=0; i<2; i++)
				DeleteExinfoInSubTree(info,infotype,cnode->child[i]);
			return ;
		}
		if(!ifExinfoShouldbeInCell(info,infotype ,cnode ))
			return ;
		if(!cnode->lpwpinfo) return;
		std::list<WpInfo *>::iterator ite,iten;
		for( ite=cnode->lpwpinfo->begin();ite!=cnode->lpwpinfo->end(); ite=iten)
		{
			iten=ite; iten++;
			if((*ite)->info==info&&(*ite)->infotype==infotype)
			{
				if(--(*ite)->rcount<=0) delete *ite;
				cnode->lpwpinfo->erase(ite);
			}
		}
		if((cnode->lpwpinfo)->empty())
		{
			delete cnode->lpwpinfo;
			cnode->lpwpinfo=0;
		}
	}
	//删除点后检查合并子树
	void Kodtree::CheckAndMergeSubTreeAfterDelete(const Point<3> &p,CellNode *cnode)
	{

		if(!cnode||cnode->IsLeaf()|| ! IsPointOverlapWithBox(p,cnode->bound,root->bound ,epsoverlap ))
			return;
		else if(cnode->numvert<=Kodtree::cellcapacity)
		{
			MergeSubTree(cnode);
			CheckAndRemoveSurplusWpInfoAfterMerge(cnode);
		}else
			for(int i=0; i<2; i++)
				CheckAndMergeSubTreeAfterDelete(p,cnode->child[i]);
	}
	//合并后检查删除额外的信息
	void Kodtree::CheckAndRemoveSurplusWpInfoAfterMerge(CellNode *cnode)
	{

		if(!cnode->lpwpinfo) return;
		std::list<WpInfo *>::iterator ite,iten;
		for( ite=cnode->lpwpinfo->begin();ite!=cnode->lpwpinfo->end(); ite=iten)
		{
			iten=ite; iten++;
			if(!ifExinfoShouldbeInCell((*ite)->info,(*ite)->infotype,cnode))
			{
				if(--(*ite)->rcount<=0) delete *ite;
				cnode->lpwpinfo->erase(ite);
			}
		}
		if((cnode->lpwpinfo)->empty())
		{
			delete cnode->lpwpinfo;
			cnode->lpwpinfo=0;
		}
	}

	//合并子树
	void Kodtree::MergeSubTree(CellNode *cnode)
	{

		if(cnode==0) jf_error("err mergecellup");
		if(cnode->IsLeaf()) return;
		for(int i=0; i<2; i++)
			MergeSubTree(cnode->child[i]);
		Merge2SubCellWpVert(cnode);
		Merge2SubCellWpInfo(cnode);
		for(int i=0; i<2; i++)
		{
			delete cnode->child[i];
			cnode->child[i]=0;
		}
	}
	//合并两个单元中的节点
	void Kodtree::Merge2SubCellWpVert(CellNode *cnode)
	{

		cnode->vert =(PtWpVert *) new PtWpVert[Kodtree::cellcapacity];
		if(cnode->IsLeaf()) jf_error("err merge2subcellvert");
		for(int i=0; i<cnode->child[0]->numvert; i++)
		{
			cnode->vert[i]=cnode->child[0]->vert[i];
			cnode->vert[i]->vget=true;
			cnode->vert[i]->rcount++;
		}
		int count=cnode->child[0]->numvert;
		for(int i=0; i<cnode->child[1]->numvert; i++)
		{
			WpVert *v=cnode->child[1]->vert[i];
			if(v->vget ==false)
			{cnode->vert[count++]=v; 
			v->rcount++;
			}
			//		else v->rcount--;
		}
		for(int i=0; i<count; i++)
			cnode->vert[i]->vget=false;
		if(cnode->numvert!=count) jf_error("err merge2subcellvert1");
	}
	//合并两个单元中的额外信息
	void Kodtree::Merge2SubCellWpInfo(CellNode *cnode)
	{
		if(cnode->IsLeaf()) jf_error("err merge2subcellwpinfo");
		CellNode *left=cnode->child[0],*right=cnode->child[1];
		if(left->lpwpinfo==0&&right->lpwpinfo==0)
		{
			cnode->lpwpinfo=0;
			return;
		}
		if(left->lpwpinfo!=0)
		{
			if(right->lpwpinfo!=0)
			{
				for(std::list<WpInfo *>::iterator ite=left->lpwpinfo->begin();ite!=left->lpwpinfo->end(); ite++)
				{
					(*ite)->get=true;
					//(*ite)->rcount++;
				}
				for(std::list<WpInfo *>::iterator iten, ite=right->lpwpinfo->begin();ite!=right->lpwpinfo->end(); ite=iten)
				{
					iten=ite, iten++;
					if(!(*ite)->get)
					{
						//(*ite)->rcount++;
						left->lpwpinfo->splice(left->lpwpinfo->end(),*(right->lpwpinfo),ite);
					}
				}
				for(std::list<WpInfo *>::iterator ite=left->lpwpinfo->begin();ite!=left->lpwpinfo->end(); ite++)
					(*ite)->get=false;
			}
			cnode->lpwpinfo=left->lpwpinfo;
			left->lpwpinfo=0;
			//delete right->lpwpinfo;
			//right->lpwpinfo=0;
		}else
		{
			cnode->lpwpinfo=right->lpwpinfo;
			right->lpwpinfo=0;
		}
	}

	//搜集盒子中的点
	void Kodtree::CollectVertsWithBox(const Box<3> &bd, list<void *> &lvert)
	{

		std::list<WpVert *> lwpvert;
		CollectWpVertsWithBoxInSubTree(root,bd,lwpvert); // may be not unique.
		for(std::list<WpVert *>::iterator ite=lwpvert.begin();ite!=lwpvert.end(); ite++)
		{
			lvert.push_back((*ite)->vt);
			(*ite)->vget=false;
		}
	}
	//搜集盒子中的点
	void Kodtree::CollectWpVertsWithBoxInSubTree(CellNode *cnode,const Box<3> &bd,list<WpVert *> &lvert)
	{
		if(!cnode) return;
		if(!IsTwoBoxOverlap(bd,cnode->bound)) return;
		if(!cnode->IsLeaf())
		{
			CollectWpVertsWithBoxInSubTree(cnode->child[0],bd,lvert);
			CollectWpVertsWithBoxInSubTree(cnode->child[1],bd,lvert);
		}else
		{
			for(int i=0; i<cnode->numvert; i++)
			{
				Point<3> p;
				if(cnode->vert[i]->vget==true) continue;
				pofv(p,cnode->vert[i]->vt);
				if(IsBoxContainPoint(p,bd,bd))
				{
					lvert.push_back(cnode->vert[i]); //maybe not unique.
					cnode->vert[i]->vget =true;
				}
			}
		}
	}
	//搜集单元中的点
	void Kodtree::CollectVertsWithCell(CellNode *cnode, std::vector<void *> &vecvert)
	{

		for(int i=0; i<cnode->numvert; i++)
			vecvert.push_back(cnode->vert[i]->vt);
	}

	//搜集盒子中的额外信息
	void Kodtree::CollectExinfoWithBox(const Box<3> &bd, int infotype,list<void *> &lexinfo)
	{

		std::list<WpInfo *> lwpinfo;
		CollectWpinfoWithBoxInSubTree(root,bd,infotype,lwpinfo);
		for(std::list<WpInfo *>::iterator ite=lwpinfo.begin();ite!=lwpinfo.end(); ite++)
		{
			lexinfo.push_back((*ite)->info);
			(*ite)->get=false;
		}
	}
	//搜集盒子中的额外信息
	void Kodtree::CollectWpinfoWithBoxInSubTree(CellNode *cnode,const Box<3> &bd,int infotype,list<WpInfo *> &lwpinfo)
	{
		if(!cnode) return;
		if(!IsTwoBoxOverlap(bd,cnode->bound)) return;
		if(!cnode->IsLeaf())
		{
			CollectWpinfoWithBoxInSubTree(cnode->child[0],bd,infotype,lwpinfo);
			CollectWpinfoWithBoxInSubTree(cnode->child[1],bd,infotype,lwpinfo);
		}else
		{
			if(cnode->lpwpinfo==0) return;
			for(std::list<WpInfo *>::iterator ite=cnode->lpwpinfo->begin();ite!=cnode->lpwpinfo->end(); ite++)
			{
				if((*ite)->infotype!=infotype||(*ite)->get==true)
					continue;
				if(ifExinfoOverlapBox((*ite)->info,infotype,bd,epsoverlap))
				{
					lwpinfo.push_back((*ite));
					(*ite)->get=true;
				}

			}
		}
	}
	//搜集单元中的额外信息
	void Kodtree::CollectExinfoWithCell(CellNode *cnode, int infotype,list<void *> &lexinfo)
	{

		if(cnode->lpwpinfo==0) return;
		for(std::list<WpInfo *>::iterator ite=cnode->lpwpinfo->begin();ite!=cnode->lpwpinfo->end(); ite++)
			if((*ite)->infotype==infotype)
				lexinfo.push_back((*ite)->info);
	}
	//判断节点是否可以分裂
	bool Kodtree::CanNodeSplit(CellNode *cnode)
	{
		return (cnode->bound[3]-cnode->bound[0])*100000.> root->bound[3]-root->bound[0];
	}
	//分裂节点
	void Kodtree::SplitNode(CellNode *cnode)
	{

		for(int i=0; i<2; i++){
			cnode->child[i]=new CellNode(cnode->bound);
			cnode->child[i]->parent=cnode;
		}
		int di;
		GetTheLongestDistanceOfBox(cnode->bound,di);
		cnode->child[1]->bound[di]=cnode->child[0]->bound[di+3]=(cnode->bound[di]+cnode->bound[di+3])/2.;
		//	if(cnode->vert==0)
		//		return;
		for(int i=0; i<cnode->numvert; i++)
		{
			Point<3> p;
			pofv(p,cnode->vert[i]->vt);
			for(int j=0; j<2; j++)
				InsertWpVertInSubTree(p,cnode->vert[i],cnode->child[j]);
		}
		for(int i=0; i<cnode->numvert ; i++)
			cnode->vert[i]->rcount--;
		delete [] cnode->vert;
		cnode->vert=0;
		if(cnode->lpwpinfo==0) return;
		for(std::list<WpInfo *>::iterator ite=cnode->lpwpinfo->begin();ite!=cnode->lpwpinfo->end(); ite++)
		{
			(*ite)->rcount--;
			for(int i=0; i<2; i++)
				InsertWpInfoInSubTree(*ite,cnode->child[i]);
		}
		delete cnode->lpwpinfo;
		cnode->lpwpinfo=0;
		//	cnode->numvert=0;
	}

	//释放子树
	void Kodtree::FreeSubTree(CellNode *pcell)
	{

		if(pcell ==0) return;
		for(int i=0; i<2; i++)
			FreeSubTree(pcell->child[i]);
		delete pcell;
	}


	CellNode ::CellNode(const Box<3> &bd)
	{

		//	psegar=0;
		vert=0;
		numvert=0;
		nodecapacity=0;
		lpwpinfo=0;
		inoutattrib=-3;
		for(int i=0; i<6; i++)
			bound[i]=bd[i];
		child[0]=child[1]=0;
		parent=0;
	}

	CellNode ::~CellNode ()
	{
		//	delete psegar;
		if(vert!=0)
			for(int i=0; i<numvert; i++)
				if(--(vert[i]->rcount)<=0) delete vert[i];
		if(lpwpinfo!=0)
			for(std::list<WpInfo *>::iterator ite=lpwpinfo->begin();ite!=lpwpinfo->end(); ite++)
				if(--(*ite)->rcount<=0)	delete *ite;
		delete [] vert;
		delete lpwpinfo;
	}

}
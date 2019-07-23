#ifndef _KODTREE_
#define _KODTREE_

/* *************************************************************************/
/* File:   kodtree.hpp                                                 */
/* Function:   accelerate the search speed of vetices and triangles            */ 
/* Date:   01. march. 2016                                                   */
/* Class:  WpInfo      WpVert      CellNode     Kodtree  */
/* *************************************************************************/

namespace meshwork
{
	//额外的信息
	struct WpInfo
	{
		void *info;
		int infotype;
		bool get;
		int rcount;
		WpInfo( void *inf,int infot):info(inf),infotype(infot),get(false),rcount(0){}
	};
	//顶点的信息
	struct WpVert
	{
		void *vt;
		bool vget;
		int rcount;
		WpVert(void *vin):vt(vin),vget(false),rcount(0){}
		WpVert():vt(0),vget(false),rcount(0){}
	};
	typedef WpVert * PtWpVert;
	struct CellNode;
	template<typename   T>   
	T* renew(T* &p,size_t oldsize,size_t newsize)
	{   
		T* tmp=new T[newsize];   
		memcpy(tmp,p,oldsize*sizeof(T));   
		delete [] p;   
		return p=tmp;   
	}   

	class Kodtree{
	public:
		typedef void (*Funcpointofvert)(Point<3> &p,void *v);
		typedef bool (*Funcexinfoshouldbeincell)(void *info, int infotype, CellNode *cnode);
		typedef bool (*Funcexinfooverlapbox) (void *info, int infotype, const Box<3> &bd,double epsi);
		Kodtree(const Box<3> &bd,Funcpointofvert pofv, const int capacity=10,double epsi=0);
		Kodtree(void **vti, int numvi,Funcpointofvert pofv,const int capacity=10,double epsi=0);
		Kodtree(void **vti, int numvi,const Box<3> &bd,Funcpointofvert pofv, const int capacity=10,double epsi=0);
		~Kodtree();
		CellNode *GetRoot(void)
		{ 
			return root;
		}
		double GetEpsCell(void)
		{
			return epscell;
		}
		WpVert *InsertVert(void *v)
		{
			WpVert *nvert= new WpVert(v);
			Point<3> p;
			pofv(p,v);
			InsertWpVertInSubTree(p,nvert,root);
			if(nvert->rcount==0)
			{
				delete nvert;
				return 0;
			}
			return nvert;
		}
		bool DeleteVert(void *v)
		{
			Point<3> p;
			pofv(p,v);
			if(!IsVertRecordInSubTree(p,v,root)) return false;
			DeleteVertInSubTree(p,v,root);
			CheckAndMergeSubTreeAfterDelete(p,root);
			return true;
		}
		WpInfo * InsertExinfo(void *info,int infotype)
		{
			WpInfo *nwinf=new WpInfo(info,infotype);
			InsertWpInfoInSubTree(nwinf,root);
			if(nwinf->rcount==0)
			{
				delete nwinf;
				return 0;
			}
			return nwinf;
		}
		void DeleteExinfo(void *info,int infotype)
		{
			DeleteExinfoInSubTree(info,infotype,root);
		}
		void CollectVertsWithBox(const Box<3> &bd, std::list<void *> &lvert);
		void CollectVertsWithCell(CellNode *cnode, std::vector<void *> &vecvert);
		void CollectExinfoWithBox(const Box<3> &bd, int infotype,std::list<void *> &linfo);
		void CollectExinfoWithCell(CellNode *cnode, int infotype,std::list<void *> &lexinfo);
		void SetFuncExinfoShouldbeInCell(Funcexinfoshouldbeincell infunc) 
		{
			ifExinfoShouldbeInCell=infunc;
		}
		void SetFuncExinfoOverlapBox(Funcexinfooverlapbox infunc)
		{
			ifExinfoOverlapBox=infunc;
		}
		CellNode *FindaLeafCellContainingPoint(CellNode *pcell,Point<3> p);
		CellNode *FindTheNearestAncestorContainingPoint(CellNode *pcell,Point<3> pcha);
		void FreeSubTree(CellNode *pcell);
		Kodtree();
	private:
		bool IsVertRecordInSubTree( const Point<3> &p, void *v,CellNode *cnode);
		void InsertWpVertInSubTree( const Point<3> &p, WpVert *nv,CellNode *cnode);
		void InsertWpInfoInSubTree(WpInfo *pwinfo, CellNode *cnode);
		void CollectWpVertsWithBoxInSubTree(CellNode *cnode,const Box<3> &bd,std::list<WpVert *> &lvert);
		void CollectWpinfoWithBoxInSubTree(CellNode *cnode,const Box<3> &bd,int infotype,std::list<WpInfo *> &lwpinfo);
		void DeleteVertInSubTree(const Point<3> &p,void *v,CellNode *cnode);
		void DeleteExinfoInSubTree(void *info,int infotype, CellNode *cnode);
		void CheckAndRemoveSurplusWpInfoAfterMerge(CellNode *cnode);
		void CheckAndMergeSubTreeAfterDelete(const Point<3> &p,CellNode *cnode);
		void MergeSubTree(CellNode *cnode);
		void Merge2SubCellWpVert(CellNode *cnode);
		void Merge2SubCellWpInfo(CellNode *cnode);
		bool If2CellNeighb(CellNode *pcell0, CellNode *pcell1);
		bool CanNodeSplit(CellNode *cnode);
		void SplitNode(CellNode *cnode);

	private:
		//	static const double epsilonon;
		//	static const double epscoplanar;
		double epsoverlap;
		int cellcapacity;
		Funcpointofvert pofv;
		Funcexinfoshouldbeincell ifExinfoShouldbeInCell;
		Funcexinfooverlapbox ifExinfoOverlapBox;
		double epscell;
		CellNode *root;
	};

	struct CellNode
	{
		WpVert **vert;
		int numvert;
		int nodecapacity;
		std::list<WpInfo *> *lpwpinfo;
		Box<3> bound;
		CellNode *child[2];
		CellNode *parent;
		int inoutattrib;
		CellNode(const Box<3> &bd);
		~CellNode();
		CellNode *AnotherChild(CellNode *pcell)
		{
			if(pcell==0) return this;
			if(child[0]==pcell) return child[1];
			else return child[0];
		}
		bool IsEmpty()
		{ 
			return vert==0;
		}
		bool IsLeaf()
		{ 
			return !child[0];
		}
	};

}
#endif //_kodtree_
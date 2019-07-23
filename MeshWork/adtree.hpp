#ifndef FILE_ADTREE
#define FILE_ADTREE

/* *************************************************************************/
/* File:   adtree.hh                                                            */
/* Date:   12. jan. 2016                                                     */
/* *************************************************************************/

namespace meshwork
{

	/**
	Alternating Digital Tree
	*/

	/* ******************************* ADTreeNode3 ******************************* */
	/* 存储点的二叉树，和刘老师的比较还是有很大差别的，因为只是用来存储点 */
	class ADTreeNode3
	{
	public:
		ADTreeNode3 *left, *right, *father;
		double sep;
		double data[3];
		int pi;
		int nchilds;

		ADTreeNode3 ();
		void DeleteChilds ();
		friend class ADTree3;
	};


	class ADTree3
	{
		ADTreeNode3 * root;
		double cmin[3], cmax[3];
		vector<ADTreeNode3*> ela;

	public:
		ADTree3 (const double acmin[3], const double  acmax[3]);
		~ADTree3 ();

		void Insert (const double  p[3], int pi);
		void GetIntersecting (const double  bmin[3], const double  bmax[3],vector<int> & pis) const;

		void DeleteElement (int pi);

		void Print (ostream & ost) const
		{ 
			PrintRec (ost, root); 
		}

		void PrintRec (ostream & ost, const ADTreeNode3 * node) const;
	};

	/* ******************************* Point3dTree ******************************* */

	class Point3dTree 
	{
		ADTree3 * tree;

	public:
		Point3dTree (const Point<3> & pmin, const Point<3> & pmax);
		~Point3dTree ();
		void Insert (const Point<3> & p, int pi);
		void DeleteElement (int pi) 
		{ 
			tree->DeleteElement(pi);
		}
		void GetIntersecting (const Point<3> & pmin, const Point<3> & pmax, vector<int> & pis) const;
		const ADTree3 & Tree() const
		{ 
			return *tree; 
		}

	};

	/* ******************************* ADTree6 ******************************* */
	/* 这个结构将box作为节点来存储，表达的是box的问题，例如可以用来检查两条边挨得是否够近，可以先检查包含着两条边的box是否挨得近 */
	class ADTreeNode6
	{
	public:
		ADTreeNode6 *left, *right, *father;
		double sep;
		double data[6];
		int pi;
		int nchilds;

		ADTreeNode6 ();
		void DeleteChilds ();
		friend class ADTree6;

	};

	class ADTree6
	{
		ADTreeNode6 * root;
		double cmin[6], cmax[6];
		vector<ADTreeNode6*> ela;

	public:
		ADTree6 (const double * acmin, const double * acmax);
		~ADTree6 ();

		void Insert (const double * p, int pi);
		void GetIntersecting (const double * bmin, const double * bmax,vector<int> & pis) const;

		void DeleteElement (int pi);

		void Print (ostream & ost) const
		{ 
			PrintRec (ost, root); 
		}
		int Depth () const
		{ 
			return DepthRec (root); 
		}
		int Elements () const
		{ 
			return ElementsRec (root); 
		}

		void PrintRec (ostream & ost, const ADTreeNode6 * node) const;
		int DepthRec (const ADTreeNode6 * node) const;
		int ElementsRec (const ADTreeNode6 * node) const;

		void PrintMemInfo (ostream & ost) const;
	};

	/* ******************************* Box3dTree ******************************* */
	class Box3dTree
	{
		ADTree6 * tree;
		Point<3> boxpmin, boxpmax;
	public:
		Box3dTree (const Box<3> & abox);
		Box3dTree (const Point<3> & apmin, const Point<3> & apmax);
		~Box3dTree ();
		void Insert (const Point<3> & bmin, const Point<3> & bmax, int pi);
		void Insert (const Box<3> & box, int pi)
		{
			Insert (box.PMin(), box.PMax(), pi);
		}
		void DeleteElement (int pi) 
		{ 
			tree->DeleteElement(pi); 
		}
		void GetIntersecting (const Point<3> & pmin, const Point<3> & pmax, vector<int> & pis) const;

		const ADTree6 & Tree() const 
		{ 
			return *tree;
		}
	};

}

#endif

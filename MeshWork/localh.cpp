#include"externalcommon.h"
#include"internalcommon.h"


namespace meshwork
{

	GradingBox :: GradingBox (const double * ax1, const double * ax2)
	{
		h2 = 0.5 * (ax2[0] - ax1[0]);
		for (int i = 0; i < 3; i++)
			xmid[i] = 0.5 * (ax1[i] + ax2[i]);

		for (int i = 0; i < 8; i++)
			childs[i] = NULL;
		father = NULL;

		hopt = 2 * h2;
	}

	void GradingBox :: DeleteChilds()
	{
		for (int i = 0; i < 8; i++)
			if (childs[i])
			{
				childs[i]->DeleteChilds();
				delete childs[i];
				childs[i] = NULL;
			}
	}


	LocalH :: LocalH (const Point<3> & pmin, const Point<3> & pmax, double agrading)
	{
		double x1[3], x2[3];
		double hmax;

		boundingbox = Box<3> (pmin, pmax);
		grading = agrading;

		// a small enlargement, non-regular points 
		double val = 0.0879;
		for (int i = 1; i <= 3; i++)
		{
			x1[i-1] = (1 + val * i) * pmin(i-1) - val * i * pmax(i-1);
			x2[i-1] = 1.1 * pmax(i-1) - 0.1 * pmin(i-1);
		}

		hmax = x2[0] - x1[0];
		for (int i = 1; i <= 2; i++)
			if (x2[i] - x1[i] > hmax)
				hmax = x2[i] - x1[i];

		for (int i = 0; i <= 2; i++)
			x2[i] = x1[i] + hmax;

		root = new GradingBox (x1, x2);
		boxes.push_back (root);
	}


	LocalH :: LocalH (const Box<3> & box, double agrading)
	{
		Point<3> pmin = box.PMin();
		Point<3> pmax = box.PMax();

		double x1[3], x2[3];
		double hmax;

		boundingbox = Box<3> (pmin, pmax);
		grading = agrading;

		// a small enlargement, non-regular points 
		double val = 0.0879;
		for (int i = 1; i <= 3; i++)
		{
			x1[i-1] = (1 + val * i) * pmin(i-1) - val * i * pmax(i-1);
			x2[i-1] = 1.1 * pmax(i-1) - 0.1 * pmin(i-1);
		}

		hmax = x2[0] - x1[0];
		for (int i = 1; i <= 2; i++)
			if (x2[i] - x1[i] > hmax)
				hmax = x2[i] - x1[i];

		for (int i = 0; i <= 2; i++)
			x2[i] = x1[i] + hmax;

		root = new GradingBox (x1, x2);
		boxes.push_back (root);
	}




	LocalH :: ~LocalH ()
	{
		root->DeleteChilds();
		delete root;
	}

	void LocalH :: Delete ()
	{
		root->DeleteChilds();
	}

	//在p点设置局部尺寸
	void LocalH :: SetH (const Point<3> & p, double h)
	{
		/*
		(*testout) << "Set h at " << p << " to " << h << endl;
		if (h < 1e-8)
		{
		cout << "do not set h to " << h << endl;
		return;
		}
		*/

		//点落在了包容盒的外面，所以无法继续计算
		if (fabs (p[0] - root->xmid[0]) > root->h2 ||
			fabs (p[1] - root->xmid[1]) > root->h2 ||
			fabs (p[2] - root->xmid[2]) > root->h2)
			return;

		/*      
		if (p.X() < root->x1[0] || p.X() > root->x2[0] ||
		p.Y() < root->x1[1] || p.Y() > root->x2[1] ||
		p.Z() < root->x1[2] || p.Z() > root->x2[2])
		return;
		*/

		//盒子的localH小于1.2倍的h，即不操作
		if (GetH(p) <= 1.2*h) return;
	//	if (GetH(p) <= 1.2*h) return;

		GradingBox * box = root;
		GradingBox * nbox = root;
		GradingBox * ngb;
		int childnr;
		double x1[3], x2[3];

		while (nbox)
		{
			box = nbox;
			childnr = 0;
			if (p[0] > box->xmid[0]) childnr += 1;
			if (p[1] > box->xmid[1]) childnr += 2;
			if (p[2] > box->xmid[2]) childnr += 4;
			nbox = box->childs[childnr];
		};
		//while (1.6 * box->h2 > h)
		while (2 * box->h2 > h)
		{
			childnr = 0;
			if (p[0] > box->xmid[0]) childnr += 1;
			if (p[1] > box->xmid[1]) childnr += 2;
			if (p[2] > box->xmid[2]) childnr += 4;

			double h2 = box->h2;
			if (childnr & 1)
			{
				x1[0] = box->xmid[0];
				x2[0] = x1[0]+h2;   // box->x2[0];
			}
			else
			{
				x2[0] = box->xmid[0];
				x1[0] = x2[0]-h2;   // box->x1[0];
			}

			if (childnr & 2)
			{
				x1[1] = box->xmid[1];
				x2[1] = x1[1]+h2;   // box->x2[1];
			}
			else
			{
				x2[1] = box->xmid[1];
				x1[1] = x2[1]-h2;   // box->x1[1];
			}

			if (childnr & 4)
			{
				x1[2] = box->xmid[2];
				x2[2] = x1[2]+h2;  // box->x2[2];
			}
			else
			{
				x2[2] = box->xmid[2];
				x1[2] = x2[2]-h2;  // box->x1[2];
			}

			ngb = new GradingBox (x1, x2);
			box->childs[childnr] = ngb;
			ngb->father = box;

			boxes.push_back (ngb);
			box = box->childs[childnr];
		}

		box->hopt = h;


		double hbox = 2 * box->h2;  // box->x2[0] - box->x1[0];
		double hnp = h + grading * hbox;

		Point<3> np;
		for (int i = 0; i < 3; i++)
		{
			np = p;
			np[i] = p[i] + hbox;
			SetH (np, hnp);

			np[i] = p[i] - hbox;
			SetH (np, hnp);
		}
	}


	//获取点的局部尺寸，通过八叉树二分查找
	double LocalH :: GetH (const Point<3> & x) const
	{
		const GradingBox * box = root;

		while (1)
		{
			int childnr = 0;
			if (x[0] > box->xmid[0]) childnr += 1;
			if (x[1] > box->xmid[1]) childnr += 2;
			if (x[2] > box->xmid[2]) childnr += 4;

			if (box->childs[childnr])
				box = box->childs[childnr];
			else
				return box->hopt;
		}
	}


	/// minimal h in box (pmin, pmax)
	double LocalH :: GetMinH (const Point<3> & pmin, const Point<3> & pmax) const
	{ 
		Point<3> pmin2, pmax2;
		for (int j = 0; j < 3; j++)
			if (pmin[j] < pmax[j])
			{ pmin2[j] = pmin[j]; pmax2[j] = pmax[j]; }
			else
			{ pmin2[j] = pmax[j]; pmax2[j] = pmin[j]; }

			return GetMinHRec (pmin2, pmax2, root); 
	}


	double LocalH :: GetMinHRec (const Point<3> & pmin, const Point<3> & pmax,
		const GradingBox * box) const
	{
		double h2 = box->h2;
		if (pmax[0] < box->xmid[0]-h2 || pmin[0] > box->xmid[0]+h2 ||
			pmax[1] < box->xmid[1]-h2 || pmin[1] > box->xmid[1]+h2 ||
			pmax[2] < box->xmid[2]-h2 || pmin[2] > box->xmid[2]+h2)
			return 1e8;

		double hmin = 2 * box->h2; // box->x2[0] - box->x1[0];

		for (int i = 0; i < 8; i++)
			if (box->childs[i])
				hmin = min (hmin, GetMinHRec (pmin, pmax, box->childs[i]));

		return hmin;
	}

	void LocalH :: PrintMemInfo (ostream & ost) const
	{
		ost << "LocalH: " << boxes.size() << " boxes of " << sizeof(GradingBox)
			<< " bytes = " << boxes.size()*sizeof(GradingBox) << " bytes" << endl;
	}
}

#include"externalcommon.h"
#include"internalcommon.h"

namespace meshwork
{

  /* ******************************* ADTree3 ******************************* */

  ADTreeNode3 :: ADTreeNode3()
  {
    pi = -1;
    left = NULL;
    right = NULL;
    father = NULL;
    nchilds = 0;
  }

  void ADTreeNode3 :: DeleteChilds ()
  {
	  if (left)
	  {
		  left->DeleteChilds();
		  delete left;
		  left = NULL;
	  }
	  if (right)
	  {
		  right->DeleteChilds();
		  delete right;
		  right = NULL;
	  }
  }

  ADTree3 ::  ADTree3 (const double acmin[3], const double  acmax[3]): ela(0)
  {
	  for ( int i = 0; i < 3; i++)
	  {
		  cmin[i] =acmin[i];
		  cmax[i]=acmax[i];
	  }

    root = new ADTreeNode3;
    root->sep = (cmin[0] + cmax[0]) / 2;
  }

  ADTree3 :: ~ADTree3 ()
  {
    root->DeleteChilds();
    delete root;
  }

  void ADTree3 :: Insert (const double  p[3], int pi)
  {
	  ADTreeNode3 *node(NULL);
	  ADTreeNode3 *next;
	  int dir;
	  int lr(0);

	  double bmin[3];
	  double bmax[3];

	  for ( int i = 0; i < 3; i++)
	  {
		  bmin[i] =cmin[i];
		  bmax[i]=cmax[i];
	  }

	  next = root;
	  dir = 0;
	  while (next)
	  {
		  node = next;
		  if (node->pi == -1)
		  {    
			  for ( int i = 0; i < 3; i++)
			  {
				  node->data[i] =p[i];
			  }
			  node->pi = pi;
			  if (ela.size() < pi+1)
				  ela.resize (pi+1);
			  ela[pi] = node;
			  return;
		  }

		  if (node->sep > p[dir])
		  {
			  next = node->left;
			  bmax[dir] = node->sep;
			  lr = 0;
		  }
		  else
		  {
			  next = node->right;
			  bmin[dir] = node->sep;
			  lr = 1;
		  }
		  dir++;
		  if (dir == 3)
			  dir = 0;
	  }

	  next = new ADTreeNode3;
	  for ( int i = 0; i < 3; i++)
	  {
		  next->data[i] =p[i];
	  }
	  next->pi = pi;
	  next->sep = (bmin[dir] + bmax[dir]) / 2;

	  if (ela.size() < pi+1)
		  ela.resize (pi+1);
	  ela[pi] = next;		

	  if (lr)
		  node->right = next;
	  else
		  node->left = next;

	  next -> father = node;
	  while (node)
	  {
		  node->nchilds++;
		  node = node->father;
	  }
  }

  void ADTree3 :: DeleteElement (int pi)
  {
	  ADTreeNode3 * node = ela[pi];

	  node->pi = -1;

	  node = node->father;
	  while (node)
	  {
		  node->nchilds--;
		  node = node->father;
	  }
  }

  void ADTree3::GetIntersecting (const double  bmin[3], const double  bmax[3], vector<int> & pis) const
  {
	  vector<ADTreeNode3*> stack;
	  vector<int> stackdir;
	  pis.resize(0);
	  ADTreeNode3 * node;
	  int dir, stacks;

	  stack.push_back (root);
	  stackdir.push_back(0);
	  stacks = 1;

	  while (stacks)
	  {
		  node = stack[stacks-1];
		  dir = stackdir[stacks-1]; 
		  stacks--;
		  if (node->pi != -1)
		  {
			  if (node->data[0] >= bmin[0] && node->data[0] <= bmax[0] &&
				  node->data[1] >= bmin[1] && node->data[1] <= bmax[1] &&
				  node->data[2] >= bmin[2] && node->data[2] <= bmax[2])

				  pis.push_back (node->pi);
		  }

		  int ndir = dir+1;
		  if (ndir == 3)
			  ndir = 0;

		  if (node->left && bmin[dir] <= node->sep)
		  {
			  stacks++;
			  if(stack.size()<stacks)
			  {
				  stack.resize(stacks);
				  stackdir.resize(stacks);
			  }
			  stack[stacks-1]=(node->left);
			  stackdir[stacks-1]=(ndir);
		  }
		  if (node->right && bmax[dir] >= node->sep)
		  {
			  stacks++;
			  if(stack.size()<stacks)
			  {
				  stack.resize(stacks);
				  stackdir.resize(stacks);
			  }
			  stack[stacks-1]=(node->right);
			  stackdir[stacks-1]=(ndir);
		  }
	  }
  }

  void ADTree3 :: PrintRec (ostream & ost, const ADTreeNode3 * node) const
  {
	  if (node->data)
	  {
		  ost << node->pi << ": ";
		  ost << node->nchilds << " childs, ";
		  for (int i = 0; i < 3; i++)
			  ost << node->data[i] << " ";
		  ost << endl;
	  }
	  if (node->left)
		  PrintRec (ost, node->left);
	  if (node->right)
		  PrintRec (ost, node->right);
  }


  /* ************************************* Point3dTree ********************** */

  Point3dTree :: Point3dTree (const Point<3> & pmin, const Point<3> & pmax)
  {
	  double pmi[3], pma[3];
	  for (int i = 0; i < 3; i++)
	  {
		  pmi[i] = pmin(i);
		  pma[i] = pmax(i);
	  }
	  tree = new ADTree3 (pmi, pma);
  }

  Point3dTree :: ~Point3dTree ()
  {
    delete tree;
  }

  void Point3dTree :: Insert (const Point<3> & p, int pi)
  {
    double pd[3];
    pd[0] = p(0);
    pd[1] = p(1);
    pd[2] = p(2);
    tree->Insert (pd, pi);
  }

  void Point3dTree :: GetIntersecting (const Point<3> & pmin, const Point<3> & pmax, vector<int> & pis) const
  {
	  double pmi[3], pma[3];
	  for (int i = 0; i < 3; i++)
	  {
		  pmi[i] = pmin(i);
		  pma[i] = pmax(i);
	  }
	  tree->GetIntersecting (pmi, pma, pis);
  }

 /* ******************************* ADTree6 ******************************* */

   ADTreeNode6 :: ADTreeNode6()
  {
    pi = -1;
    left = NULL;
    right = NULL;
    father = NULL;
    nchilds = 0;
  }

  void ADTreeNode6 :: DeleteChilds ()
  {
	  if (left)
	  {
		  left->DeleteChilds();
		  delete left;
		  left = NULL;
	  }
	  if (right)
	  {
		  right->DeleteChilds();
		  delete right;
		  right = NULL;
	  }
  }

  ADTree6 :: ADTree6 (const double * acmin, const double * acmax) : ela(0)
  {
	  	for ( int i = 0; i < 6; i++)
	{
		cmin[i] =acmin[i];
		cmax[i]=acmax[i];
	}
    root = new ADTreeNode6;
    root->sep = (cmin[0] + cmax[0]) / 2;
  }

  ADTree6 :: ~ADTree6 ()
  {
    root->DeleteChilds();
    delete root;
  }

  void ADTree6 :: Insert (const double * p, int pi)
  {
	  ADTreeNode6 *node(NULL);
	  ADTreeNode6 *next;
	  int dir;
	  int lr(0);

	  double bmin[6];
	  double bmax[6];
	  for ( int i = 0; i < 6; i++)
	  {
		  bmin[i] =cmin[i];
		  bmax[i]=cmax[i];
	  }
	  next = root;
	  dir = 0;
	  while (next)
	  {
		  node = next;
		  if (node->pi == -1)
		  {    
			  for ( int i = 0; i < 6; i++)
			  {
				  node->data[i] =p[i];
			  }
			  node->pi = pi;
			  if (ela.size() < pi+1)
				  ela.resize (pi+1);
			  ela[pi] = node;
			  return;
		  }
		  if (node->sep > p[dir])
		  {
			  next = node->left;
			  bmax[dir] = node->sep;
			  lr = 0;
		  }
		  else
		  {
			  next = node->right;
			  bmin[dir] = node->sep;
			  lr = 1;
		  }
		  dir++;
		  if (dir == 6) dir = 0;
	  }

	  next = new ADTreeNode6;
	  for ( int i = 0; i < 6; i++)
	  {
		  next->data[i] =p[i];
	  }
	  next->pi = pi;
	  next->sep = (bmin[dir] + bmax[dir]) / 2;

	  if (ela.size() < pi+1)
		  ela.resize (pi+1);
	  ela[pi] = next;

	  if (lr)
		  node->right = next;
	  else
		  node->left = next;
	  next -> father = node;
	  while (node)
	  {
		  node->nchilds++;
		  node = node->father;
	  }
  }

  void ADTree6 :: DeleteElement (int pi)
  {
	  ADTreeNode6 * node = ela[pi];
	  node->pi = -1;
	  node = node->father;
	  while (node)
	  {
		  node->nchilds--;
		  node = node->father;
	  }
  }

  void ADTree6 :: PrintMemInfo (ostream & ost) const
  {
    ost << Elements() << " elements a " << sizeof(ADTreeNode6) 
	<< " Bytes = "
	<< Elements() * sizeof(ADTreeNode6) << endl;
    ost << "maxind = " << ela.size() << " = " << sizeof(ADTreeNode6*) * ela.size() << " Bytes" << endl;
  }

  void ADTree6 :: GetIntersecting (const double * bmin, const double * bmax,vector<int> & pis) const
  {
	  // static Array<inttn6> stack(10000);
	  // stack.SetSize (10000);
	  //ArrayMem<inttn6,10000> stack(10000);
	  vector<ADTreeNode6*> stack;
	  vector<int> stackdir;
	  pis.resize(0);
	  stack.push_back( root);
	  stackdir.push_back( 0);
	  int stacks = 1;
	  while (stacks)
	  {
		  ADTreeNode6 * node = stack[stacks-1];
		  int dir = stackdir[stacks-1]; 
		  stacks--;
		  if (node->pi != -1)
		  {
			  if (node->data[0] > bmax[0] || 
				  node->data[1] > bmax[1] || 
				  node->data[2] > bmax[2] || 
				  node->data[3] < bmin[3] || 
				  node->data[4] < bmin[4] || 
				  node->data[5] < bmin[5])
				  ;
			  else
			  {
				  pis.push_back (node->pi);
			  }
		  }

		  int ndir = (dir+1) % 6;
		  if (node->left && bmin[dir] <= node->sep)
		  {
			  stacks++;
			  if(stack.size()<stacks)
			  {
				  stack.resize(stacks);
				  stackdir.resize(stacks);
			  }
			  stack[stacks-1]=(node->left);
			  stackdir[stacks-1]=(ndir);
		  }
		  if (node->right && bmax[dir] >= node->sep)
		  {
			  stacks++;
			  if(stack.size()<stacks)
			  {
				  stack.resize(stacks);
				  stackdir.resize(stacks);
			  }
			  stack[stacks-1]=(node->right);
			  stackdir[stacks-1]=(ndir);
		  }
	  }
  }

  void ADTree6 :: PrintRec (ostream & ost, const ADTreeNode6 * node) const
  {
	  if (node->data)
	  {
		  ost << node->pi << ": ";
		  ost << node->nchilds << " childs, ";
		  for (int i = 0; i < 6; i++)
			  ost << node->data[i] << " ";
		  ost << endl;
	  }
	  if (node->left)
		  PrintRec (ost, node->left);
	  if (node->right)
		  PrintRec (ost, node->right);
  }

  int ADTree6 :: DepthRec (const ADTreeNode6 * node) const
  {
    int ldepth = 0;
    int rdepth = 0;

    if (node->left)
      ldepth = DepthRec(node->left);
    if (node->right)
      rdepth = DepthRec(node->right);
    return 1 + max (ldepth, rdepth);
  }

  int ADTree6 :: ElementsRec (const ADTreeNode6 * node) const
  {
    int els = 1;
    if (node->left)
      els += ElementsRec(node->left);
    if (node->right)
      els += ElementsRec(node->right);
    return els;
  }

 /* ******************************* Box3dTree ******************************* */
  Box3dTree :: Box3dTree (const Box<3> & abox)
  {
	  boxpmin = abox.PMin();
	  boxpmax = abox.PMax();
	  double tpmin[6], tpmax[6];
	  for (int i = 0; i < 3; i++)
	  {
		  tpmin[i] = tpmin[i+3] = boxpmin(i);
		  tpmax[i] = tpmax[i+3] = boxpmax(i);
	  }
	  tree = new ADTree6 (tpmin, tpmax);
  }

  Box3dTree :: Box3dTree (const Point<3> & apmin, const Point<3> & apmax)
  {
	  boxpmin = apmin;
	  boxpmax = apmax;
	  double tpmin[6], tpmax[6];
	  for (int i = 0; i < 3; i++)
	  {
		  tpmin[i] = tpmin[i+3] = boxpmin(i);
		  tpmax[i] = tpmax[i+3] = boxpmax(i);
	  }
	  tree = new ADTree6 (tpmin, tpmax);
  }

  Box3dTree :: ~Box3dTree ()
  {
    delete tree;
  }

  void Box3dTree :: Insert (const Point<3> & bmin, const Point<3> & bmax, int pi)
  {
	  double tp[6];
	  for (int i = 0; i < 3; i++)
	  {
		  tp[i] = bmin(i);
		  tp[i+3] = bmax(i);
	  }
	  tree->Insert (tp, pi);
  }

  void Box3dTree ::GetIntersecting (const Point<3> & pmin, const Point<3> & pmax, vector<int> & pis) const
  {
	  double tpmin[6];
	  double tpmax[6];

	  for (int i = 0; i < 3; i++)
	  {
		  tpmin[i] = boxpmin(i);
		  tpmax[i] = pmax(i);

		  tpmin[i+3] = pmin(i);
		  tpmax[i+3] = boxpmax(i);
	  }

	  tree->GetIntersecting (tpmin, tpmax, pis);
  }

}
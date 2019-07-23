#include"externalcommon.h"
#include"internalcommon.h"
namespace meshwork
{
   Mesh:: Mesh()
   {
		lochfunc=NULL;
		hglob = 1e10;
		hmin = 0;
   }

   void Mesh :: SetGlobalH (double h)
   {
	   hglob = h;
   }

   void Mesh :: SetLocalH (const Point<3> & pmin, const Point<3> & pmax, double grading)
   {
	   Point<3> c = Center (pmin, pmax);
	   double d = max (pmax[0]-pmin[0],
		   max(pmax[1]-pmin[1],
		   pmax[2]-pmin[2]));
	   d /= 2;
	   Point<3> pmin2 = c - Vec<3> (d, d, d);
	   Point<3> pmax2 = c + Vec<3> (d, d, d);

	   if(lochfunc!=NULL)
	   delete lochfunc;
	   lochfunc = new LocalH (pmin2, pmax2, grading);
   }

   void Mesh :: RestrictLocalH (const Point<3> & p, double hloc)
   {
	   if(hloc < hmin)
		   hloc = hmin;

	   //cout << "restrict h in " << p << " to " << hloc << endl;
	   if (!lochfunc)
	   {
		   cout<<("RestrictLocalH called, creating mesh-size tree")<<endl;
	   }
	   //在p点设置局部尺寸
	   lochfunc -> SetH (p, hloc);
   }

   void Mesh :: RestrictLocalHLine (const Point<3> & p1, const Point<3> & p2,double hloc)
   {
	   if(hloc < hmin)
		   hloc = hmin;

	   // cout << "restrict h along " << p1 << " - " << p2 << " to " << hloc << endl;
	   int i;
	   int steps = int (Dist (p1, p2) / hloc) + 2;
	   Vec<3> v=p2-p1;

	   for (i = 0; i <= steps; i++)
	   {
		   Point<3> p = p1 + (double(i)/double(steps) * v);
		   RestrictLocalH (p, hloc);
	   }
   }

   double Mesh :: GetH (const Point<3> & p) const
   {
	   double hmin = hglob;
	   if (lochfunc)
	   {
		   double hl = lochfunc->GetH (p);
		   if (hl < hglob)
			   hmin = hl;
	   }
	   return hmin;
   }

  int Mesh ::AddPoint(MeshPoint mp)
  {
	  //for(int i = 0; i<nodes.size(); i++)
	  //{
		 // if(nodes[i]==mp)
			//  return i;
	  //}
	  nodes.push_back(mp);
	  return nodes.size()-1;
  }

  int Mesh ::AddElement2d(Element2d ele)
  {
	  surfelements.push_back(ele);
	  return surfelements.size()-1;
  }

  int Mesh ::AddFaceDescriptor(FaceDescriptor fd)
  {
	  facedecoding.push_back(fd);
	  return facedecoding.size()-1;
  }


}
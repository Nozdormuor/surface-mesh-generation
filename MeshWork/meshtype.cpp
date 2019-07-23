#include"externalcommon.h"
#include"internalcommon.h"
namespace meshwork
{

	 Element2d :: Element2d ()
	 {
		 for(int i = 0 ; i < 3; i++)
		 {
			 pnum[i]=-1;
		 }
		 findex=-1;
	 }

	  Element2d ::  Element2d (int pi1, int pi2, int pi3)
	   {
		   pnum[0]=pi1;
		   pnum[1]=pi2;
		   pnum[3]=pi3;
		   findex=-1;
	   }

	  Element2d :: Element2d (int pi1, int pi2, int pi3, int face)
	  {
		  pnum[0]=pi1;
		  pnum[1]=pi2;
		  pnum[3]=pi3;
		  findex=face;
	  }

	  void Element2d :: GetBox (const vector<MeshPoint> & points, Box<3> & box) const
	  {
		  box.Set (points[pnum[0]]);
		  for (unsigned i = 1; i < 3; i++)
			  box.Add (points[pnum[i]]);
	  }

	  inline void Element2d :: Invert()
	  {
		  int temp =pnum[2];
		  pnum[2]=pnum[3];
		  pnum[3]=temp;

	  }

	  inline void Element2d :: NormalizeNumbering ()
	  {

			  if (PNum(0) < PNum(1) && PNum(0) < PNum(2))
				  return;
			  else
			  {
				  if (PNum(1) < PNum(2))
				  {
					  int pi1 = PNum(1);
					  PNum(1) = PNum(2);
					  PNum(2) = PNum(0);
					  PNum(0) = pi1;
				  }
				  else
				  {
					  int pi1 = PNum(2);
					  PNum(2) = PNum(1);
					  PNum(1) = PNum(0);
					  PNum(0) = pi1;
				  }
			  }

	  }

	  ostream & operator<<(ostream  & s, const Element2d & el)
	  {
		  for (int j = 0; j < 3; j++)
			  s << " " << el.PNum(j);
		  return s;
	  }

	  FaceDescriptor ::  FaceDescriptor()
	  { 
		  surfnr = domin = domout  = -1; 

	  }

	  FaceDescriptor ::FaceDescriptor(int surfnri)
	  {
		  surfnr = surfnri;
		  domin = domout  = -1; 
	  }

	  FaceDescriptor ::  FaceDescriptor(const FaceDescriptor& other)
		  : surfnr(other.surfnr), domin(other.domin), domout(other.domout)
	  { 
		  ;
	  }
	  ostream & operator<<(ostream  & s, const FaceDescriptor & fd)
	  {
		  s << "surfnr = " << fd.SurfNr() 
			  << ", domin = " << fd.DomainIn()
			  << ", domout = " << fd.DomainOut();
		  for(int i = 0; i < fd.element.size(); i++)
			  cout<<" "<<fd.element[i]<<" ";
		  return s;
	  }

}
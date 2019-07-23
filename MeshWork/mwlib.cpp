#include"externalcommon.h"
#include"internalcommon.h"
namespace meshwork
{
	 void MW_STL_GenerateSurfaceMesh(STLGeometry *pgeom, Mesh *pmesh)
	{
		int retval = MW_STLSurfaceMeshing (*pgeom, *pmesh);
	}

	 void MW_WriteSTLFormat (Mesh *mesh, const char *filname)
	{
		cout << "\nWrite STL Surface Mesh" << endl;

		ostream *outfile;
		string filename(filname);

		outfile = new ofstream(filename.c_str());

		int i;

		outfile->precision(10);

		*outfile << "solid" << endl;

		for (i = 0; i < mesh->GetNSE(); i++)
		{
			*outfile << "facet normal ";
			const Point<3>& p1 = mesh->GetPoint(mesh->GetSurfaceElement(i).PNum(0));
			const Point<3>& p2 = mesh->GetPoint(mesh->GetSurfaceElement(i).PNum(1));
			const Point<3>& p3 = mesh->GetPoint(mesh->GetSurfaceElement(i).PNum(2));

			Vec<3> normal = Cross(p2-p1,p3-p1);
			if (normal.Length() != 0)
			{
				normal /= (normal.Length());
			}

			*outfile << normal[0]<< " " << normal[1] << " " << normal[2] << "\n";
			*outfile << "outer loop\n";

			*outfile << "vertex " << p1[0] << " " << p1[1] << " " << p1[2] << "\n";
			*outfile << "vertex " << p2[0] << " " << p2[1] << " " << p2[2] << "\n";
			*outfile << "vertex " << p3[0] << " " << p3[1] << " " << p3[2] << "\n";

			*outfile << "endloop\n";
			*outfile << "endfacet\n";
		}
		*outfile << "endsolid" << endl;
	}



}
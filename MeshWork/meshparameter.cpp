#include"externalcommon.h"
#include"internalcommon.h"
namespace meshwork
{

	MeshingParameters :: MeshingParameters ()
	{
		optimize3d = "cmdmustm";
		//optimize3d = "cmdmstm";
		optsteps3d = 3;
		optimize2d = "smsmsmSmSmSm";
		optsteps2d = 3;
		opterrpow = 2;
		blockfill = 1;
		filldist = 0.1;
		safety = 5;
		relinnersafety = 3;
		uselocalh = 1;
		grading = 0.12;
		delaunay = 1;
		maxh = 10;
		minh = 0;
		meshsizefilename = NULL;
		startinsurface = 0;
		checkoverlap = 1;
		checkoverlappingboundary = 1;
		checkchartboundary = 1;
		curvaturesafety = 2;
		segmentsperedge = 1;
		parthread = 0;

		elsizeweight = 0.2;
		giveuptol2d = 200;
		giveuptol = 10;
		maxoutersteps = 10;
		starshapeclass = 5;
		baseelnp = 0;
		sloppy = 1;

		badellimit = 175;
		check_impossible = 0;
		secondorder = 0;
	}

	void MeshingParameters :: Print (ostream & ost) const
	{
		ost << "Meshing parameters: " << endl
			<< "optimize3d = " << optimize3d << endl
			<< "optsteps3d = " << optsteps3d << endl
			<< " optimize2d = " <<  optimize2d << endl
			<< " optsteps2d = " <<  optsteps2d << endl
			<< " opterrpow = " <<  opterrpow << endl
			<< " blockfill = " <<  blockfill << endl
			<< " filldist = " <<  filldist << endl
			<< " safety = " <<  safety << endl
			<< " relinnersafety = " <<  relinnersafety << endl
			<< " uselocalh = " <<  uselocalh << endl
			<< " grading = " <<  grading << endl
			<< " delaunay = " <<  delaunay << endl
			<< " maxh = " <<  maxh << endl;
		if(meshsizefilename)
			ost << " meshsizefilename = " <<  meshsizefilename << endl;
		else
			ost << " meshsizefilename = NULL" << endl;
		ost << " startinsurface = " <<  startinsurface << endl
			<< " checkoverlap = " <<  checkoverlap << endl
			<< " checkchartboundary = " <<  checkchartboundary << endl
			<< " curvaturesafety = " <<  curvaturesafety << endl
			<< " segmentsperedge = " <<  segmentsperedge << endl
			<< " parthread = " <<  parthread << endl
			<< " elsizeweight = " <<  elsizeweight << endl
			<< " giveuptol2d = " <<  giveuptol2d << endl
			<< " giveuptol = " <<  giveuptol << endl
			<< " maxoutersteps = " <<  maxoutersteps << endl
			<< " starshapeclass = " <<  starshapeclass << endl
			<< " baseelnp        = " <<  baseelnp        << endl
			<< " sloppy = " <<  sloppy << endl
			<< " badellimit = " <<  badellimit << endl
			<< " secondorder = " <<  secondorder << endl
			<< " elementorder = " <<  elementorder << endl
			<< " quad = " <<  quad << endl
			<< " inverttets = " <<  inverttets << endl
			<< " inverttrigs = " <<  inverttrigs << endl;
	}

	void MeshingParameters :: CopyFrom(const MeshingParameters & other)
	{
		//strcpy(optimize3d,other.optimize3d); 
		optimize3d = other.optimize3d;
		optsteps3d = other.optsteps3d;
		//strcpy(optimize2d,other.optimize2d); 
		optimize2d = other.optimize2d;
		optsteps2d = other.optsteps2d;
		opterrpow = other.opterrpow;
		blockfill = other.blockfill;
		filldist = other.filldist;
		safety = other.safety;
		relinnersafety = other.relinnersafety;
		uselocalh = other.uselocalh;
		grading = other.grading;
		delaunay = other.delaunay;
		maxh = other.maxh;
		//strcpy(const_cast<char*>(meshsizefilename), other.meshsizefilename);
		//const_cast<char*>(meshsizefilename) = other.meshsizefilename; //???
		startinsurface = other.startinsurface;
		checkoverlap = other.checkoverlap;
		checkoverlappingboundary = other.checkoverlappingboundary;
		checkchartboundary = other.checkchartboundary;
		curvaturesafety = other.curvaturesafety;
		segmentsperedge = other.segmentsperedge;
		parthread = other.parthread;
		elsizeweight = other.elsizeweight;
		giveuptol2d = other.giveuptol2d;
		giveuptol = other.giveuptol;
		maxoutersteps = other.maxoutersteps;
		starshapeclass = other.starshapeclass;
		baseelnp = other.baseelnp;       
		sloppy = other.sloppy;
		badellimit = other.badellimit;
		secondorder = other.secondorder;
		elementorder = other.elementorder;
		quad = other.quad;
		inverttets = other.inverttets;
		inverttrigs = other.inverttrigs;
	}

	STLParameters ::   STLParameters()
	{
		yangle = 30;
		contyangle = 20;
		edgecornerangle = 60;
		chartangle = 15;
		outerchartangle = 70;

		usesearchtree = 0;
		atlasminh = 1E-4;
		resthsurfcurvfac = 2.0;
		resthsurfcurvenable = 1;
		resthatlasfac = 2;
		resthatlasenable = 1;
		resthchartdistfac = 1.2;
		resthchartdistenable = 0;
		resthlinelengthfac = 2.0;
		resthlinelengthenable = 1;
		resthcloseedgefac = 1.0;
		resthcloseedgeenable = 0.5;
		resthedgeanglefac = 2.0;
		resthedgeangleenable = 1;
		resthsurfmeshcurvfac = 1.;
		resthsurfmeshcurvenable = 0;
		recalc_h_opt = 1;
	}

	void STLParameters :: Print (ostream & ost) const
	{
		ost << "STL parameters:" << endl
			<< "yellow angle = " << yangle << endl
			<< "continued yellow angle = " << contyangle << endl
			<< "edgecornerangle = " << edgecornerangle << endl
			<< "chartangle = " << chartangle << endl
			<< "outerchartangle = " << outerchartangle << endl
			<< "restrict h due to ..., enable and safety factor: " << endl
			<< "surface curvature: " << resthsurfcurvenable
			<< ", fac = " << resthsurfcurvfac << endl
			<< "atlas surface curvature: " << resthatlasenable
			<< ", fac = " << resthatlasfac << endl
			<< "chart distance: " << resthchartdistenable
			<< ", fac = " << resthchartdistfac << endl
			<< "line length: " << resthlinelengthenable
			<< ", fac = " << resthlinelengthfac << endl
			<< "close edges: " << resthcloseedgeenable
			<< ", fac = " << resthcloseedgefac << endl
			<< "edge angle: " << resthedgeangleenable
			<< ", fac = " << resthedgeanglefac << endl;
	}

}
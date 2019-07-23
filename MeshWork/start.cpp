#include"externalcommon.h"
#include"internalcommon.h"
#include<iostream>
#include<cstring>
#include<fstream>
#include <time.h>
using namespace std;
using namespace meshwork;
void main (){
	string strMeshFile = "finallamp.stl";
	ifstream ist(strMeshFile);
	if(!ist){
		cout<<"file error"<<endl;
		return;
	}
	mparam.maxh=2.5;
	mparam.minh=0.20;
	STLGeometry geom;
	geom.Load(ist);
	Mesh mh;
	//mparam.maxh=20;
	mh.SetGlobalH(mparam.maxh);
	clock_t ts=clock();
	MW_STLSurfaceMeshing(geom,mh);
	clock_t te=clock();
	cout<<0<<" tests "<<(te-ts)/1000.<<" sec."<<endl;
	MW_WriteSTLFormat(&mh, "remeshlamp.stl");
	cout<<"complete!!!!"<<endl;
}

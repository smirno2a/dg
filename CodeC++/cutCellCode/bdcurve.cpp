#include <stdio.h>
#include <math.h>
#include "mImportExport.h"
#include "mDGMesh.h"
#include "mEntity.h"
#include "Vertex.h"
#include "Face.h"
#include "Edge.h"
#include "mMesh.h"
#include <time.h>
#include <map>
//using namespace std;

int bdcurve(double x, double y){
  double R=1.0;
	double eps=0.0000000001;
	if(x*x+y*y>R*R+eps)
		return 1;
	else if (x*x+y*y<R*R-eps)
		return -1;
	else return 0;
}

mPoint mDGMesh::intersection(mPoint* A, Vertex* B){
	//find the intersection with boundary curve
	double R=1.0;
	 // mPoint aa=A->point();
	mPoint aa=*A;	
	mPoint bb=B->point();
		double k=(bb(1)-aa(1))/(bb(0)-aa(0));
		double b=aa(1)-k*aa(0);
		double c=sqrt(4.0*k*k*b*b-4*(k*k+1)*(b*b-R*R));
		double r1=(-2.0*k*b+c)/2.0/(k*k+1);
		double r2=(-2.0*k*b-c)/2.0/(k*k+1);
		mPoint newv;
		if ((r1>aa(0)&&r1<bb(0))||(r1<aa(0)&&r1>bb(0))){
		    double	ry1=sqrt(R*R-r1*r1);
			if (( ry1 > aa(1) && ry1 < bb(1) )||( ry1 < aa(1) && ry1 > bb(1) ))
				return newv = mPoint (r1,ry1,0);
			else
				return newv = mPoint (r1,-ry1,0);
		}
		else {
   		  double  ry1=sqrt(R*R-r2*r2);
		    if (( ry1 > aa(1) && ry1 < bb(1) )||( ry1 < aa(1) && ry1 > bb(1) ))
            return newv = mPoint (r2,ry1,0);
			  else
				    return newv = mPoint (r2,-ry1,0);
		}
}
#include "DGAnalysis.h"
#include "FieldEvaluator.h"
#include "mPoint.h"
#include "ScalarLaw.h"
#include "LinearSystem.h"
#include "DGLimiter.h"
#include "DGSensor.h"
#include "mDGMesh.h"
#include <stdio.h>
#include <math.h>
#include "mVector.h"
#include "Constants.h"

class ExactOne : public FieldEvaluator
{
public :
  void eval (const mPoint &space, double time, double *val) const
  {
	 //if (space(1)-time*.402<-.625*(space(0)-.666*time)-0.8&& space(0)<2.) val[0]=1.; else val[0]=0; 
	//  val[0] = pow(sin(mPi*(2.*space(0)+1.*space(1)-3.*time)),1);
    //val[0] = 1;
	  //val[0]=pow(space(0)-2.*space(1)+space(2),1);
	  //val[0]=space(0)-2.*space(1);
     val[0]=0.0;
    if (max(fabs(space(0) + 0. - time*1.0),fabs(space(1) + 0. - time*1.0) ) <= 0.25 ) //(sqrt(pow(space(0) + 0.0 - time*0.,2) + pow(space(1) + 0.0 - time*0.,2)) < 0.2)// 
      val[0] = 1.0;
    else
    {
       val[0] = 0.0;
    }
    /*double x = space(0), y = space(1);
    //Rotating shapes
     double r , x0 = -0.5, y0 = 0, R = 0.15;
     r = sqrt(pow(x-x0,2) + pow(y,2));
     val[0] = 0.0;
     if (max(fabs(space(0) - 0.35),fabs(space(1) ) ) <= 0.25) 
      val[0] = 1.0;
     else if (r <= 0.25)
     {
       val[0] = pow(cos(2.*M_PI*r),2);
     }*/

   // val[0] = 1.*x*x*x + 0.5*y*y*y + 1.5*x*x + 1.*y*y + 3.*x + 2.*y + 0.1;//space(0) - space(1);
  //val[0] = (x - 1.*time) + (y - 0.*time);
  // double r , x0 = -0.25, y0 = -0.25, R = 0.15;
  // r = sqrt(pow(x-x0 - 1.*time ,2) + pow(y-y0 - 1.*time,2));
  // val[0] = 2.5*exp(-0.5*r*r/(R*R));

  //val[0] = sin(2*M_PI*(space(0) + 1.*space(1) - time*1.));
   /* double r=sqrt(pow(space(0)+0.45 -time*1.0,2) + space(1)*space(1));
	 if (r<=0.35) val[0] = 1-r/0.35; 
	 else val[0]=0.0;*/
   // double r=sqrt(pow(space(0)  + 0.4 - time*1.,2) + space(1)*space(1));
  //  if (r<=0.5) val[0] = pow(cos(mPi*r),2); 
  //  else val[0]=0.0;
     //val[0] = pow(2.*space(0)-3.*space(1)+space(2),4);
     //val[0] = pow(2.*space(0)-3.*space(1)+space(2),2);
	// val[0] =sin(mPi*(-time+space(0)-2.*space(1)));  
  /*double r=sqrt(pow(space(0)+0.5 - time*1.,2) + space(1)*space(1));
    if (r<=0.4) val[0] = pow(cos((1./1.6)*2*mPi*r),2); 
    else val[0]=0.0;*/
  }
  int order() const
  {
    return 2;
  }
};

// Initial Conditiong

class CaseOne : public FieldEvaluator
{
public :
  void eval (const mPoint &space, double time, double *val) const
  {
     // if (space(0)>0.1 && space(0)<0.6 && space(1)>-0.25 && space(1) <0.25) val[0] = 1.0;
	  //else 
	  //val[0]=0.0;
    /*if (sqrt(pow(space(0)+0.0,2) + pow(space(1) + 0.0,2)) < 0.2)//(max(fabs(space(0) + 0.75),fabs(space(1) + 0.) ) <= 0.2 )//
      val[0] = 1.0;
    else
    {
       val[0] = 0.0;
    }*/
    //val[0] = space(0) + space(1);
   // val[0] = sin(2*M_PI*(space(0) - time*1.));
    //val[0] = cos(2*M_PI*(space(0) - time*1.));
    /*double r=sqrt(pow(space(0) + 0.4 - time*1.,2) + space(1)*space(1));
    if (r<=0.5) val[0] = pow(cos(mPi*r),2); 
    else val[0]=0.0;*/
    
/*	 else   
	 {
	 double r=sqrt(pow(space(0)+0.45,2) + space(1)*space(1));
	 if (r<=0.35) val[0] = 1-r/0.35; 
	 else val[0]=0.0;
	 }*/
  /* double r=sqrt(pow(space(0)+0.45,2) + space(1)*space(1));
	 if (r<=0.35) val[0] = 1-r/0.35; 
	 else val[0]=0.0;*/
    
  /*double r=sqrt(pow(space(0)+0.5,2) + space(1)*space(1));
    if (r<=0.4) val[0] = pow(cos((1./1.6)*2*mPi*r),2); 
    else val[0]=0.0;*/
    //val[0] = pow(space(0),2)*pow(space(1),3)+pow(space(0),3)*pow(space(1),2)*3. +space(1)*space(1)+space(0)*space(0)+space(1)+space(0);
  }
  int order() const
  {
    return 3;
  }
};

class VelocityOne : public FieldEvaluator
{
public :
  void eval (const mPoint &space, double time, double *val) const
  {
	 val[0]= 1.0;
	  val[1]= 1.0;
    //val[0] =2.;
    //val[1] =1.;
  //val[0] = -2.*mPi*space(1);
  // val[1] = 2.*mPi*space(0);
    val[2] = 0.0;
  }
  int order() const
  {
    return 1;
  }
};

void ScalarPb (mDGMesh *theMesh, int order)
{

  CaseOne i;
  VelocityOne v;
  ExactOne e;

  //DGLineSensor sensor1  ("t1.dat",mPoint(-1,0),mPoint(1,0),256,1.0);
  //DGLineSensor sensor2  ("t2.dat",mPoint(-1,-1),mPoint(1,1),256,.1);
  //DGPointSensor sensor2 ("t2.dat",mPoint(0.0,0.0),1.0);
  
  //DGVertexLimiter myLimiter;
  //DGQuadMomentLimiter myLimiter;
  //DGBarthLimiter myLimiter;
  //DGSuperBee myLimiter;
  // DGMomentLimiter myLimiter;
  //ScalarLaw theLaw (&i,&v,&i);
  ScalarLaw theLaw (&e,&v,&e);
  //DGAnalysis analysis(theMesh,&theLaw,&myLimiter,order);
   DGAnalysis analysis(theMesh,&theLaw,0,order);
   //analysis.addSensor(&sensor1);
  // analysis.addSensor(&sensor2);
  DGGmshSensor sensor0  ("square",0.05);
  analysis.addSensor(&sensor0);
   //analysis.addSensor(&sensor2);
  analysis.run();
}

/********************Linear System****************************/

class Matrix1 : public FieldEvaluator
{
public :
  void eval (const mPoint &space, double time, double *val) const
  {
    val[0] =8.;
    val[1] = -9.;
	val[2] =0.;
    val[3] =6.;
    val[4] = -13.;
	val[5] = 0.;
	val[6] = 0.;
	val[7] = 0.;
	val[8] = 0.;
	  double u=-3.;
	  double rho=1;
	  double K=1.;
/*	 val[0] =u;
    val[1] = K;
	val[2] = 0;
    val[3] = 1./rho;
    val[4] = u;
	val[5] = 0.;
	val[6] = 0.;
	val[7] = 0.;
	val[8] = u;	*/
  }
  int order() const { return 1;}
};

class Matrix2 : public FieldEvaluator
{
public :
  void eval (const mPoint &space, double time, double *val) const
  {
    val[0] =13.;
    val[1] = -9.;
	val[2] = 0;
    val[3] = 6.;
    val[4] = -8.;
	val[5] = 0.;
	val[6] = 0.;
	val[7] = 0.;
	val[8] = 0.;
	  double v=-2.;
	  double rho=1.;
	  double K=1.;
	/* val[0] =v;
    val[1] = 0;
	val[2] = K;
    val[3] = 0;
    val[4] = v;
	val[5] = 0.;
	val[6] = 1./rho;
	val[7] = 0.;
	val[8] = v;*/
  }
  int order() const { return 1;}
};

class ExactSolution : public FieldEvaluator
{
public :
  void eval (const mPoint &space, double time, double *val) const
  {
    val[0]=sin(0.5*(space(0)+space(1)));
    val[1]=cos(0.5*(space(0)+space(1)));
	val[2]=0.0; 
	/*val[0]=sin(0.5*(space(0)+space(1)));
    val[1]=space(0)*space(0);
	val[2]=space(1)*space(1);*/
	//val[0]=val[1]=val[2]=1.;
  }
  int order() const {return 2;}
};

void LinSystem (mDGMesh *theMesh, int order)
{
  Matrix1 a1;
  Matrix2 a2;
  ExactSolution e;

  DGLineSensor sensor1  ("t1.dat",mPoint(0,1),mPoint(1,0),3000,5.);
  DGPointSensor sensor2 ("t2.dat",mPoint(0.0,0.0),1.0);
  
  // DGVertexLimiter myLimiter;
  // DGBarthLimiter myLimiter;
  //DGSuperBee myLimiter;

  LinearSystem theLaw (&e,&e,&a1,&a2);
  
   //DGAnalysis analysis(theMesh,&theLaw,&myLimiter,order);
     DGAnalysis analysis(theMesh,&theLaw,0,order);
  DGGmshSensor sensor0  ("squar",5.0);
  analysis.addSensor(&sensor0);
  analysis.addSensor(&sensor1);
  
  analysis.run();
}


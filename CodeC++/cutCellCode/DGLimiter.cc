#include "mEntity.h"
#include "Vertex.h"
#include "Mapping.h"
#include "Integrator.h"
#include "ConservationLaw.h"
#include "FunctionSpace.h"
#include "DGLimiter.h"
#include "DGCell.h"
#include "FieldEvaluator.h"
#include "Constants.h"
#include "gEntity.h"
#include <stdio.h>
#include <math.h>
#include <list>
#define absmin(a,b) (fabs(a)<fabs(b) ? a:b)
double minmod(double,double,double,bool&);
double minmod(double,double,bool &);
double minmod(double a, double b) {
	if (a*b > 0) return absmin(a,b);
	else return 0;
}

using namespace std;

static void recurGetAllsubs(mEntity *ent, list<mEntity*> &lis)
{
  int n = ent->getLevel();
  if(!ent->isAdjacencyCreated(n))lis.push_back(ent);
  else
    {
      for(int i=0;i<ent->size(n);i++)recurGetAllsubs(ent->get(n,i),lis);
    }
}


/******************************************************************************/

void DGSuperBee::limit(DGCell *vcell, double time)
{
  double Max[MaxNbEqn],Min[MaxNbEqn];
  int n = vcell->theMeshEntity->getLevel();
  int kk = 0;
  int fSize = vcell->fSize;
  int fOrder = vcell->theFunctionSpace->order();
  
  double a = 1./3.464101615137754;
  double I[2][2],II[2][2],III[2][2];
  double u,v,w,weight,val[MaxNbEqn];
  double uvol,vvol,wvol,x,y,z;
  double mean[4][5];
  
  I[0][0] = 0.5; I[0][1] = 0.5; I[1][0] = -a; I[1][1] = a;
  II[0][0] = 0.0; II[0][1] = -0.5; II[1][0] = 2*a; II[1][1] = a;
  III[0][0] = -0.5; III[0][1] = 0.0; III[1][0] = -a; III[1][1] = -2*a;
  
  for(int i=0;i<vcell->theMeshEntity->size(n-1);i++)
    {
      mEntity *bound = vcell->theMeshEntity->get(n-1,i);
      if(bound->size(n) == 2)
	{
	  DGCell *other;
	  if(bound->get(n,0) != vcell->theMeshEntity)
	    other = (DGCell*)bound->get(n,0)->getCell();
	  else
	    other = (DGCell*)bound->get(n,1)->getCell();
	  if(!kk)
	    for(int j=0;j<other->cSize;j++)
	      {
		Max[j] = Min[j] = other->theMean[j];
		mean[1][j] = other->theMean[j]; 
	      }
	  else
	    {
	      for(int j=0;j<other->cSize;j++)
		{
		  Max[j]  = (other->theMean[j] > Max[j])?other->theMean[j]:Max[j];
		  Min[j]  = (other->theMean[j] < Min[j])?other->theMean[j]:Min[j];
		  mean[kk+1][j] = other->theMean[j]; 
		}
	      // if ( Min[0]<0 || Min[3] <0 ) printf("NEGATIVE MEAN  \n");
	    }
	  kk++;
	}
    }
  
  list<mEntity *> allSubs;
  for(int k=0;k<vcell->theMeshEntity->size(n-1);k++)
  recurGetAllsubs(vcell->theMeshEntity->get(n-1,k),allSubs);
    for(list<mEntity*>::const_iterator it = allSubs.begin();it!=allSubs.end();++it)
    {
      //MODIFY!!!
      DGBoundaryCell *cell = (DGBoundaryCell*)(*it)->getCell();      
            DGCell *right = 0;
      DGCell *left = (DGCell*)cell->mleft->getCell();
      if(cell->mright)right = (DGCell*)cell->mright->getCell();
      
      // printf(" left %e \n",left->theMean[0]);
      //if(cell->mright) printf(" right %e \n",right->theMean[0]);
      int order = cell->computeOrder();
      order =1;
      GaussIntegrator gauss;
      int Nb = gauss.nbIntegrationPoints(cell->theBoundaryEntity,order);
      
      for(int i=0;i<Nb;i++)
	{
	    gauss.iPoint(cell->theBoundaryEntity,i,order,u,v,w,weight);
	    cell->theMapping->eval(u,v,w,x,y,z);
	    vcell->theMapping->invert(x,y,z,uvol,vvol,wvol);
	    vcell->interpolate(uvol,vvol,wvol,val);
	    //  printf("u,v,w %f % f %f\n", u,v,w);
	    // printf("x,y,z %f % f %f\n",x,y,z );
	    // printf("uvol,vvol,wvol %f % f %f\n", uvol,vvol,wvol); 
	    // if values are not physical
	    
	    for(int j=0;j<vcell->cSize;j++)
	      {
		if(fabs(val[j] - vcell->theMean[j]) > 1.e-5 * fabs(val[j] + vcell->theMean[j]) && fabs( vcell->theMean[j]) > 1.e-13 )	    
		  {
		    
		    //    printf("MAX %f  MIn %f Mean %f val %f \n",Max[j], Min[j], vcell->theMean[j], val[j]);
		    
		    if(val[j] > Max[j] || val[j] < Min[j] )
		      {
			if (Max[j] > vcell->theMean[j] && Min[j] > vcell->theMean[j])
			  { 
			    for (int k=1;k<fSize;k++)
			      vcell->theFieldsCoefficients->get(j,k) = 0.0; return;
			  }
			if (Max[j]<vcell->theMean[j] && Min[j]<vcell->theMean[j]) 
			  {
			    for (int k=1;k<fSize;k++)
			      vcell->theFieldsCoefficients->get(j,k) = 0.0; return;
			  }
			
			double b[2];
			double diff[4];
			diff[1] = fabs(mean[1][j]- vcell->theMean[j]);
			diff[2] = fabs(mean[2][j]- vcell->theMean[j]);
			diff[3] = fabs(mean[3][j]- vcell->theMean[j]);
			double pos =0.0;
			double neg = 0.0;
			if (diff[1]>0.0) pos += diff[1];else neg += -diff[1];
			if (diff[2]>0.0) pos += diff[2];else neg += -diff[2];
			if (diff[3]>0.0) pos += diff[3];else neg += -diff[3];
	    		//printf(" means %e %e %e %e\n",mean[1][j],mean[2][j],mean[3][j],vcell->theMean[j]);
			//printf("uvol,vvol,wvol %f % f %f\n", uvol,vvol,wvol); 
			//printf("theMean= %e val=%e\n",vcell->theMean[j],val[j]);
			//printf("Original coeff %e %e %e \n",vcell->theFieldsCoefficients[0][0],vcell->theFieldsCoefficients[0][1],vcell->theFieldsCoefficients[0][2]);
			double phiplus, phiminus;
			phiplus = (neg>pos) ? 1.0 : neg/pos; 
			phiminus = (pos>neg) ? 1.0 : pos/neg; 

			if (diff[1]>0.0) diff[1] *= phiplus;
			else diff[1]*=phiminus;
			if (diff[2]>0.0) diff[2] *= phiplus;
			else diff[2]*=phiminus;
			if (diff[3]>0.0) diff[3] *= phiplus;
			else diff[3]*=phiminus;
			
			b[0]=diff[1];
			b[1]= diff[2];
			//printf("b %e %e \n",b[0],b[1]);
			//printf("CASE 1\n");
			vcell->theFieldsCoefficients->get(j,1) = I[0][0]*b[0]+I[0][1]*b[1];
			vcell->theFieldsCoefficients->get(j,2)=I[1][0]*b[0]+I[1][1]*b[1];
			vcell->interpolate(0.0,0.5,0.0,val);
			if (val[j]<Min[j] || val[j]>Max[j])
			  {
			    vcell->theFieldsCoefficients->get(j,1)=0.0;
			    vcell->theFieldsCoefficients->get(j,2)=0.0;
			    printf("set coef. to zero at 1!!!!!!!!\n");
			  }
			// printf("coeff %e %e %e \n",vcell->theFieldsCoefficients[0][0],vcell->theFieldsCoefficients[0][1],vcell->theFieldsCoefficients[0][2]);
			
			//	printf("j= %d\n",j);
			/* vcell->interpolate(0.5,0.0,-0.5,val);
		     printf("Updated val=%e",val[0]);
		     vcell->interpolate(0.5,0.5,-0.5,val);
		     printf("   val=%e",val[0]);
		     vcell->interpolate(0.0,0.5,-0.5,val);
		     printf("   val=%e",val[0]);
		     printf("\n \n");*/
		      }
		  }
	      }
	} 
    }  
}

/******************************************************************************/

void DGBarthLimiter::limit(DGCell *vcell, double time)
{
  int fSize = vcell->fSize;
  int fOrder = vcell->theFunctionSpace->order();
  int cSize = vcell->cSize;
  int n = vcell->theMeshEntity->getLevel();
  vcell->limit=1;
  int i,kk=0;
  double phin,val[MaxNbEqn], Max[MaxNbEqn],Min[MaxNbEqn];
//  double u,v,w,uvol,vvol,wvol,x,y,z,weight;
  
  list<mEntity*>::const_iterator it,allSubs_end;
  double dx = vcell->cellSize;
  if (fabs(vcell->fullJump/(pow(dx,2*(fOrder+1))*vcell->getPerimeter()))>3. )
    vcell->limit = 1;
 
  if (vcell->limit) {     
    int nbFaces = vcell->theMeshEntity->size(n-1);
    for(i=0;i<nbFaces;i++) {
      mEntity *bound = vcell->theMeshEntity->get(n-1,i);
      if(bound->size(n) == 2)
	{
	  DGCell *other;
	  if(bound->get(n,0) != vcell->theMeshEntity)
	    other = (DGCell*)bound->get(n,0)->getCell();
	  else
	    other = (DGCell*)bound->get(n,1)->getCell();
	  if(!kk)
	    for(int j=0;j<cSize;j++)
	      Max[j] = Min[j] = other->theMean[j];
	  else
	    {
	      for(int j=0;j<cSize;j++)
		{
		  Max[j]  = (other->theMean[j] > Max[j])?other->theMean[j]:Max[j];
		  Min[j]  = (other->theMean[j] < Min[j])?other->theMean[j]:Min[j];
		}
	    }
	  kk++;
	}
    }    
    
    for(int j=0;j<vcell->cSize;j++) vcell->limSlope[j] = 1.0;
    
    allSubs_end=vcell->allSubs.end();
    boundaryGaussPoint *pg;
    for(it = vcell->allSubs.begin();it!=allSubs_end;++it) {
      DGBoundaryCell *cell = (DGBoundaryCell*)(*it)->getCell();      
      DGCell *right = 0;
      DGCell *left = (DGCell*)cell->mleft->getCell();
      if(cell->mright)right = (DGCell*)cell->mright->getCell();  
      
      int Nb = cell->getNbPtGauss();
      
      for(int i=0;i<Nb;i++)
	{
	  pg = cell->pt(i);
	  vcell->interpolate(pg->fctleft, val);
	  for(int j=0;j<cSize;j++)
	    {
	      //if(fabs(val[j] - vcell->theMean[j]) > 1.e-5 * fabs(val[j] + vcell->theMean[j]) &&fabs( vcell->theMean[j]) > 1.e-13 )
	      //if(pg->x>0.5 )
		  if(true)
 {
		
		if(val[j] > Max[j]) {
		  phin=(Max[j]-vcell->theMean[j])/(val[j]-vcell->theMean[j]);
		  vcell->limit=2;
		}
		else if(val[j] < Min[j]) {
		  phin=(Min[j]-vcell->theMean[j])/(val[j]-vcell->theMean[j]);
		  vcell->limit=2;
		}
		else {phin = 1.0;vcell->limit=3;}
	      }
	      else {phin = 1.0; vcell->limit = 4;
	      }
	      if(phin<0.0) {phin = 0.0;vcell->limit=5;}
	      vcell->limSlope[j] = (vcell->limSlope[j]<phin) ? 
		vcell->limSlope[j]:phin;
	    }
	}
    }
    //printf("alha 0 %e \n",vcell->limSlope[0]);
    while(1)
      { 
	int fSize1 = vcell->theFunctionSpace->size(1);
	int fSizei = vcell->theFunctionSpace->size(fOrder);
	int fSizej = vcell->theFunctionSpace->size(fOrder-1);
	
	for(int j=0;j<cSize;j++)
	  {	
	    if(fOrder == 1)
	      for (int k=1;k<fSize1;k++)
		vcell->theFieldsCoefficients->get(j,k) *= vcell->limSlope[j];
	    else	
	      for (int k=fSizej;k<fSizei;k++)
		if(vcell->limSlope[j] < 0.9)vcell->theFieldsCoefficients->get(j,k)= 0.0;
	  }	 
	if(fOrder == 1) break;
	fOrder --;
      }
  }  
}  


//Limiting at the center point

/*
void DGBarthLimiter::limit(DGCell *vcell, double time)
{
  int fSize = vcell->fSize;
  int fOrder = vcell->theFunctionSpace->order();
  int n = vcell->theMeshEntity->getLevel();
  vcell->limit=1;
  int i,k, kk=0;
  double u,v,w,weight,val[MaxNbEqn], Max[MaxNbEqn],Min[MaxNbEqn];
  double uvol,vvol,wvol,x,y,z,phin;
  
  list<mEntity *> allSubs;
  list<mEntity*>::const_iterator it,allSubs_end;
  double dx = vcell->getMaxH();
  if (fabs(vcell->fullJump/(pow(dx,2*(fOrder+1))*vcell->getPerimeter()))>1. )
    vcell->limit = 1;
  //printf("Jump %f, JUmp^3 %f, maxVal %f, dx %f\n",vcell->fullJump,pow(vcell->fullJump,3),vcell->maxVal,dx);
  
  if (vcell->limit)    
    {  
      int Size = vcell->theMeshEntity->size(n-1);
      for(i=0;i<Size;i++)
	{
	  mEntity *bound = vcell->theMeshEntity->get(n-1,i);
	  if(bound->size(n) == 2)
	    {
	      DGCell *other;
	      if(bound->get(n,0) != vcell->theMeshEntity)
		other = (DGCell*)bound->get(n,0)->getCell();
	      else
		other = (DGCell*)bound->get(n,1)->getCell();
	      if(!kk)
		for(int j=0;j<other->cSize;j++)
		  Max[j] = Min[j] = other->theMean[j];
	      else
		{
		  for(int j=0;j<other->cSize;j++)
		    {
		      Max[j]  = (other->theMean[j] > Max[j])?other->theMean[j]:Max[j];
		      Min[j]  = (other->theMean[j] < Min[j])?other->theMean[j]:Min[j];
		    }
		}
	      kk++;
	    }
	}    
      
      for(k=0;k<vcell->theMeshEntity->size(n-1);k++)
	recurGetAllsubs(vcell->theMeshEntity->get(n-1,k),allSubs);
      
      for(int j=0;j<vcell->cSize;j++) vcell->limSlope[j] = 1.0;
      allSubs_end=allSubs.end();
      
      for(it = allSubs.begin();it!=allSubs_end;++it)
	{if ((*it)->getClassification()->getId()!=20000) {
	DGBoundaryCell *cell = (DGBoundaryCell*)(*it)->getCell();      
	int order =1;
	GaussIntegrator gauss;
	int Nb = gauss.nbIntegrationPoints(cell->theBoundaryEntity,order);
	for(int i=0;i<Nb;i++)
	    {
	      gauss.iPoint(cell->theBoundaryEntity,i,order,u,v,w,weight);
	      cell->theMapping->eval(u,v,w,x,y,z);
	      vcell->theMapping->invert(x,y,z,uvol,vvol,wvol);
	      if (x>0.5) {
		vcell->interpolate(uvol,vvol,wvol,val);
		for(int j=0;j<vcell->cSize;j++) {
		  //		  if(fabs(val[j] - vcell->theMean[j]) > 1.e-5 * fabs(val[j] + vcell->theMean[j]) &&fabs( vcell->theMean[j]) > 1.e-13 )
		  if(true)
		    {
		      if(val[j] > Max[j])
			{
			  //printf("j=%d x=%e y=%e val=%e max=%e\n",j,x,y,val[j],Max[j]);
			  phin=(Max[j]-vcell->theMean[j])/(val[j]-vcell->theMean[j]);
			  vcell->limit=2;
			}
		      else if(val[j] < Min[j])
			{
			  //printf("j=%d x=%e y=%e val=%e min = %e\n",j,x,y,val[j],Min[j]);
			  phin=(Min[j]-vcell->theMean[j])/(val[j]-vcell->theMean[j]);
			  vcell->limit=2;
			}
		      else {phin = 1.0;vcell->limit=3;}
		    }
		  else phin = 1.0;
		  if(phin<0.0) {phin = 0.0;vcell->limit=5;}
		  vcell->limSlope[j] = (vcell->limSlope[j]<phin) ? 
		    vcell->limSlope[j]:phin;
		}
	      }
	    }
	}
	}
      //printf("alha 0 %e \n",vcell->limSlope[0]);
      while(1)
	{ 
	  int fSize1 = vcell->theFunctionSpace->size(1);
	  int fSizei = vcell->theFunctionSpace->size(fOrder);
	  int fSizej = vcell->theFunctionSpace->size(fOrder-1);
	  
	  for(int j=0;j<vcell->cSize;j++)
	    {	
	      if(fOrder == 1)
		for (int k=1;k<fSize1;k++)
		  vcell->theFieldsCoefficients->get(j,k) *= vcell->limSlope[j];
	      else	
		for (int k=fSizej;k<fSizei;k++)
		  if(vcell->limSlope[j] < 0.9)vcell->theFieldsCoefficients->get(j,k)= 0.0;
	    }	 
	  if(fOrder == 1) break;
	  fOrder --;
	}
    }  
}
*/

/******************************************************************************/

void DGBarthLimiterEuler::limit(DGCell *vcell, double time)
{
  int fSize = vcell->fSize;
  int fOrder = vcell->fOrder;
  int cSize  = vcell->cSize;
  int n = vcell->theMeshEntity->getLevel();
  vcell->limit=5;
  int j,k,kk=0;
  double rhom,rhoM,um,uM,vm,vM,pm,pM;
//  double bExtrema[MaxNbEqn],Max[MaxNbEqn],Min[MaxNbEqn];
  double val[MaxNbEqn],alfa[MaxNbEqn];
  double const GAMMA=1.4;
  double const GAMMA_1=0.4;
list<mEntity *> allSubs;
  list<mEntity*>::const_iterator it;
    if (vcell->limit)    
    {  
      
      for(int i=0;i<vcell->theMeshEntity->size(n-1);i++)
	{
	  mEntity *bound = vcell->theMeshEntity->get(n-1,i);
	  if(bound->size(n) == 2)
	    {
	      DGCell *other;
	      if(bound->get(n,0) != vcell->theMeshEntity)
		other = (DGCell*)bound->get(n,0)->getCell();
	      else
		other = (DGCell*)bound->get(n,1)->getCell();
	      if(!kk)
		{
		  rhoM = rhom = other->theMean[0];
		  uM   = um   = other->theMean[1]/rhoM;
		  vM   = vm   = other->theMean[2]/rhoM;
		  pM   = pm   = (GAMMA_1)*(other->theMean[3] - (uM*uM + vM*vM) *0.5*rhoM);
		}
	      else
		{
		  double rhoO = other->theMean[0]; 
		  double invrhoO = 1./other->theMean[0];
		  double uO =other->theMean[1]*invrhoO;  
		  double vO =other->theMean[2]*invrhoO;
		  double pO   = (GAMMA_1)*(other->theMean[3] - (uO*uO + vO*vO) *0.5*rhoO);
		  
		  rhoM = ( rhoO > rhoM ) ? rhoO : rhoM ; 
		  rhom = ( rhoO < rhom ) ? rhoO : rhom ; 
		  uM   = ( uO   >  uM  ) ?   uO : uM ; 
		  um   = ( uO   <  um  ) ?   uO : um ; 
		  vM   = ( vO   >  vM  ) ?   vO : vM ; 		  
		  vm   = ( vO   <  vm  ) ?   vO : vm ; 		  
		  pM   = ( pO   >  pM  ) ?   pO : pM ; 		  
		  pm   = ( pO   <  pm  ) ?   pO : pm ; 		  
		}
	      kk++;
	    }
	}
       
      for(j=0;j<cSize;j++)  vcell->limSlope[j]= alfa[j] =1.0;
      
      for(k=0;k<vcell->theMeshEntity->size(n-1);k++)
	recurGetAllsubs(vcell->theMeshEntity->get(n-1,k),allSubs);      
     
      for(it = vcell->allSubs.begin();it!=vcell->allSubs.end();++it)
	{
	  DGBoundaryCell *cell = (DGBoundaryCell*)(*it)->getCell();      
	  
	  for(int i=0;i<cell->nbPtGauss;i++)
	    {
	      boundaryGaussPoint *pg = cell->pt(i);
	      if(vcell==cell->left)vcell->interpolate(pg->fctleft,val);
	      else vcell->interpolate(pg->fctrght,val);
	      double u,v,p;
	      u = val[1]/val[0]; v=val[2]/val[0]; p = GAMMA_1*(val[3]-0.5*(u*u +v*v)*val[0]);
	      //if(fabs(val[j] - vcell->theMean[j]) > 1.e-5 * fabs(val[j] + vcell->theMean[j]) &&fabs( vcell->theMean[j]) > 1.e-13 )     
		{
		  if(val[0]>rhoM) alfa[0]=(rhoM-vcell->theMean[0])/(val[0]-vcell->theMean[0]);
		  if(val[0] < rhom) alfa[0] = (rhom-vcell->theMean[0])/(val[0]-vcell->theMean[0]);

		  if(u>uM) alfa[1]=(uM*val[0]-vcell->theMean[1])/(val[1]-vcell->theMean[1]);
		  if (u<um) alfa[1]=(um*val[0]-vcell->theMean[1])/(val[1]-vcell->theMean[1]);
		  
		  if(v>vM) alfa[2]=(vM*val[0]-vcell->theMean[2])/(val[2]-vcell->theMean[2]);
		  if(v<vm) alfa[2]=(vm*val[0]-vcell->theMean[2])/(val[2]-vcell->theMean[2]);	 

		  if(p>pM) alfa[3]=(pM/(GAMMA-1)+0.5*(u*u+v*v)*val[0]-vcell->theMean[3])/(val[3]-vcell->theMean[3]);  
		  if(p<pm) alfa[3]=(pm/(GAMMA-1)+0.5*(u*u+v*v)*val[0]-vcell->theMean[3])/(val[3]-vcell->theMean[3]);  
		  /*
		  {
		      phin = (Max[j]-vcell->theMean[j])/(val[j]-vcell->theMean[j]);
		      if (j==0) vcell->limit=0;
		      if (j==1) vcell->limit=1;
		      if (j==2) vcell->limit=2;
		      if (j==3) vcell->limit=3;
		      //printf("limiting %d\n",j);
		    }
		  else if(val[j] < Min[j])
		    {
		      phin = (Min[j]-vcell->theMean[j])/(val[j]-vcell->theMean[j]);
		      if (j==0) vcell->limit=0;
		      if (j==1) vcell->limit=1;
		      if (j==2) vcell->limit=2;
		      if (j==3) vcell->limit=3;
		      //vcell->limit=3;
		      //printf("limiting %d\n",j);
		      //printf("val<Min val[%d]=%e Min=%e Mean=%e\n",j,val[j],Min[j],vcell->theMean[j]);
		    }
		    */
		}
	      for(int j=0;j<cSize;j++) 
		{
		  if (alfa[j]<0.0) alfa[j] = 0.0;
		  vcell->limSlope[j] = (vcell->limSlope[j]<alfa[j])?vcell->limSlope[j]:alfa[j];
		}
	    }
	}
	 // printf("alha 0 %e \n",vcell->limSlope[0]);
     // for(int j=0;j<cSize;j++) printf("vcell->limSlope[j] %e ",vcell->limSlope[j]);
	  //printf("\n");
      while(1)
	{ 
	  int fSize1 = vcell->theFunctionSpace->size(1);
	  int fSizei = vcell->theFunctionSpace->size(fOrder);
	  int fSizej = vcell->theFunctionSpace->size(fOrder-1);
	  
	  for(int j=0;j<cSize;j++)
	    {	
	      if(fOrder == 1)
		for (int k=1;k<fSize1;k++)
		  vcell->theFieldsCoefficients ->get(j,k)*= vcell->limSlope[j];
	      else	
		for (int k=fSizej;k<fSizei;k++)
		  if(vcell->limSlope[j] < 0.9)vcell->theFieldsCoefficients->get(j,k)= 0.0;
	    }	 
	  if(fOrder == 1) break;
	  fOrder --;
	}
    }
  
}  

/******************************************************************************/

void DGQuadMomentLimiter::limit(DGCell *vcell, double time)
{
  int fSize = vcell->fSize;
  int fOrder = vcell->theFunctionSpace->order();
  int n = vcell->theMeshEntity->getLevel();
  vcell->limit=5;
  int kk=0;
//s  double u,v,w,weight,val[MaxNbEqn], Max[MaxNbEqn],Min[MaxNbEqn];
//  double uvol,vvol,wvol,x,y,z;
  double phin;
  list<mEntity *> allSubs;
  list<mEntity*>::const_iterator it;
  mEntity *ent;
  //double dx = vcell->getMaxH();
  //  if (fabs(pow(vcell->fullJump,3)/vcell->maxVal/sqrt(pow(dx,3*(fOrder+1)))/dx)>1. )
  //vcell->limit = 1;
  //printf("Jump %f, JUmp^3 %f, maxVal %f, dx %f\n",vcell->fullJump,pow(vcell->fullJump,3),vcell->maxVal,dx);
  mPoint v[4], v_neigh[4], centroid(0.,0.0);
  int no_vertex = vcell->theMeshEntity->size(0);
  int no_edges = vcell->theMeshEntity->size(1);
  for(int k=0;k<no_vertex;k++)
      {
        ent = vcell->theMeshEntity->get(0,k); 
        v[k] = ((Vertex*)ent)->point();
        centroid += v[k];
      }
      centroid = centroid * 0.25;

  /* for(int k =0 ; k < no_edges; ++k){
		ent = vcell->theMeshEntity->get(1,k);
		if(ent->getClassification()->getId() == 2 ||ent->getClassification()->getId() == 30000 || ent->getClassification()->getId() == 610){
			return;
		}
	}*/
  
if (vcell->limit)    
  {  
	 // printf("DGQuadLimiter is being called\n");

	double l3 = 1.0*sqrt(1.), l2 = 1.*sqrt(1.), l1= 1.*sqrt(1.);
	double scalingfactor = 3./(3. + fabs(vcell->alpha[0]) + fabs(vcell->alpha[1]));//1.0;// 0.75;
    double Max[MaxNbEqn],Min[MaxNbEqn];
    int n = vcell->theMeshEntity->getLevel();
    int k;
    int kk = 0;
    int fSize = vcell->fSize;
    int fOrder = vcell->theFunctionSpace->order();
	int flag = 0;
	double tolerance  = 1e-3;
	double temp_c3 = 0.0, temp_c1 = 0.0, temp_c2 = 0;
	double temp_c11 = 0.0, temp_c21 = 0;
	double temp_l11 = 0.0, temp_l21 = 0.0, temp_l12 = 0.0, temp_l22 = 0.0;
	double c31, c32, c3_f, c3_b;
	double mean = 0.0;
    
    double val[MaxNbEqn];
	double C0_neigh[8], C1_neigh[8], C2_neigh[8], Mean_neigh[12];
	bool isLimited1 = 1, isLimited2 = 1, isLimited3 = 1, isLimited4 = 1;
	double Temp[4] = {0,0,0,0};
	double Temp2[4] = {0,0,0,0};
	int Neighbor_flag = 0;

	for(int j =0 ; j < vcell->cSize; ++j){
		mean = vcell->theMean[j];
		for(int i = 0; i< 8; ++i){
			C0_neigh[i] = (vcell->Neighbours[i] != NULL)? (vcell->Neighbours[i])->theFieldsCoefficients ->getValCopy(j,0) : 0.0;
			C1_neigh[i] = (vcell->Neighbours[i] != NULL)? (vcell->Neighbours[i])->theFieldsCoefficients ->getValCopy(j,1) : 0;
			C2_neigh[i] = (vcell->Neighbours[i] != NULL)? (vcell->Neighbours[i])->theFieldsCoefficients ->getValCopy(j,2) : 0;
		    Mean_neigh[i] = (vcell->Neighbours[i] != NULL)? (vcell->Neighbours[i])->theMean[j] : 0.0;

			if (vcell->Neighbours[i] != NULL) Neighbor_flag += 1;
		}
		for(int i = 8; i< 12; ++i){
			//C0_neigh[i] = (vcell->Neighbours[i] != NULL)? (vcell->Neighbours[i])->theFieldsCoefficients ->getValCopy(j,0) : 0.0;
			//C1_neigh[i] = (vcell->Neighbours[i] != NULL)? (vcell->Neighbours[i])->theFieldsCoefficients ->getValCopy(j,1) : 0;
			//C2_neigh[i] = (vcell->Neighbours[i] != NULL)? (vcell->Neighbours[i])->theFieldsCoefficients ->getValCopy(j,2) : 0;
		    Mean_neigh[i] = (vcell->Neighbours[i] != NULL)? (vcell->Neighbours[i])->theMean[j] : 0.0;

			if (vcell->Neighbours[i] != NULL) Neighbor_flag += 1;
		}
		// To not limit boundary adjacent elements
		if (Neighbor_flag < 12) {return ;}

		// Reconstructing c3
		temp_c3 = vcell->theFieldsCoefficients ->getValCopy(j,3);
		if(vcell->Neighbours[0] && vcell->Neighbours[1]){
			c31 = vcell->coeff_alpha[0] * C1_neigh[0] + vcell->coeff_beta[0] * C2_neigh[0];
			c32 = vcell->coeff_alpha[1] * C1_neigh[1] + vcell->coeff_beta[1] * C2_neigh[1];
			c3_f = vcell->neighbour_weights[0]* c31 + vcell->neighbour_weights[1]* c32; 
			//Temp[0] = c3_f;//l3 *(vcell->scaling_derivative[0]/(sqrt(3.) * vcell->dist_intersection[0]) * (c3_f - vcell->theFieldsCoefficients ->getCopy(j,2)) );
			//temp_c3 = minmod(temp_c3, l3 *(vcell->scaling_derivative[0]/(sqrt(3.) * vcell->dist_intersection[0]) * (c3_f - vcell->theFieldsCoefficients ->getCopy(j,2)) ), isLimited1) ;
		    temp_c3 = minmod(temp_c3, l3 * (1./(sqrt(3.)) * (c3_f - vcell->theFieldsCoefficients ->getCopy(j,2)) ), isLimited1);
			Temp2[0] = (vcell->theFieldsCoefficients ->getValCopy(j,2) + sqrtf(3.)*vcell->theFieldsCoefficients ->getValCopy(j,3)*vcell->intersection_compvar[0]);
		}
		if(vcell->Neighbours[2] && vcell->Neighbours[3]){
			c31 = vcell->coeff_alpha[2] * C1_neigh[2] + vcell->coeff_beta[2] * C2_neigh[2];
			c32 = vcell->coeff_alpha[3] * C1_neigh[3] + vcell->coeff_beta[3] * C2_neigh[3];
			c3_b = vcell->neighbour_weights[2]* c31 + vcell->neighbour_weights[3]* c32; 
			//Temp[1] = c3_b;//l3 *(vcell->scaling_derivative[0]/(sqrt(3.) * vcell->dist_intersection[1]) * (-c3_b  + vcell->theFieldsCoefficients ->getCopy(j,2)) );
			//temp_c3 = minmod(temp_c3, l3 *(vcell->scaling_derivative[0]/(sqrt(3.) * vcell->dist_intersection[1]) * (-c3_b  + vcell->theFieldsCoefficients ->getCopy(j,2)) ), isLimited2) ;
            temp_c3 = minmod(temp_c3, l3 *(1./(sqrt(3.)) * (-c3_b  + vcell->theFieldsCoefficients ->getCopy(j,2)) ), isLimited2) ;
          
		   Temp2[1] = (vcell->theFieldsCoefficients ->getValCopy(j,2) + sqrtf(3.)*vcell->theFieldsCoefficients ->getValCopy(j,3)*vcell->intersection_compvar[1]);

		}

		
		if(vcell->Neighbours[4] && vcell->Neighbours[5]){
			c31 = vcell->coeff_alpha[4] * C1_neigh[4] + vcell->coeff_beta[4] * C2_neigh[4];
			c32 = vcell->coeff_alpha[5] * C1_neigh[5] + vcell->coeff_beta[5] * C2_neigh[5];
			c3_f = vcell->neighbour_weights[4]* c31 + vcell->neighbour_weights[5]* c32; 
			Temp[2] = l3 * (vcell->scaling_derivative[1]/(sqrt(3.) * vcell->dist_intersection[2]) * (c3_f - vcell->theFieldsCoefficients ->getCopy(j,1)) );
			//temp_c3 = minmod(temp_c3, l3 * (vcell->scaling_derivative[1]/(sqrt(3.) * vcell->dist_intersection[2]) * (c3_f - vcell->theFieldsCoefficients ->getCopy(j,1)) ), isLimited3) ;
            temp_c3 = minmod(temp_c3, l3 * (1./(sqrt(3.)) * (c3_f - vcell->theFieldsCoefficients ->getCopy(j,1)) ), isLimited3) ;
		
		}
		if(vcell->Neighbours[6] && vcell->Neighbours[7]){
			c31 = vcell->coeff_alpha[6] * C1_neigh[6] + vcell->coeff_beta[6] * C2_neigh[6];
			c32 = vcell->coeff_alpha[7] * C1_neigh[7] + vcell->coeff_beta[7] * C2_neigh[7];
			c3_b = vcell->neighbour_weights[6]* c31 + vcell->neighbour_weights[7]* c32; 
			Temp[3] =  l3 * (vcell->scaling_derivative[1] /(sqrt(3.) * vcell->dist_intersection[3]) * (-c3_b  + vcell->theFieldsCoefficients ->getCopy(j,1)) );
			//temp_c3 = minmod(temp_c3, l3 * (vcell->scaling_derivative[1] /(sqrt(3.) * vcell->dist_intersection[3]) * (-c3_b  + vcell->theFieldsCoefficients ->getCopy(j,1)) ), isLimited4) ;
		    temp_c3 = minmod(temp_c3, l3 * (1. /(sqrt(3.)) * (-c3_b  + vcell->theFieldsCoefficients ->getCopy(j,1)) ), isLimited4) ;

		}

		

		//flag = ( fabs(temp_c3) < 1e-5*tolerance || fabs(temp_c3 - vcell->theFieldsCoefficients->get(j,3))  > tolerance * fabs(vcell->theFieldsCoefficients->get(j,3)))?1:0;
        //flag = ( fabs(temp_c3/vcell->theFieldsCoefficients->get(j,3) - 1.)  > tolerance) ?1:0;
        //if( -0.76 < centroid(0) && centroid(0) < -0.73 &&  0.03 < centroid(1) && centroid(1) < 0.05){
		//if( -0.9 < centroid(0) && centroid(0) < -0.8 &&  -0.9 < centroid(1) && centroid(1) < -0.8){
		//if( -0.156 < centroid(0) && centroid(0) < -0.148 &&  -0.177 < centroid(1) && centroid(1) < -0.167){
		//if( -0.14 < centroid(0) && centroid(0) < -0.13 &&  0.16 < centroid(1) && centroid(1) < 0.167){
		//if( -0.284 < centroid(0) && centroid(0) < -0.26 &&  -0.8 < centroid(1) && centroid(1) < -0.76){
		//if( -0.26 < centroid(0) && centroid(0) < -0.24 &&  0.05 < centroid(1) && centroid(1) < 0.0667){
				//printf("C3 interpolation values : %.12e, %.12e, %.12e, %.12e, %.12e\n",vcell->theFieldsCoefficients->get(j,3),Temp[0], Temp[1], Temp[2], Temp[3]);
				//printf("Centroid : (%.12e,%.12e)\n",centroid(0),centroid(1));
				//printf("C3 interpolation values : %.15e, %.15e, %.15e, %.15e\n", Temp[0], Temp[1], Temp2[0], Temp2[1]);
				//printf("scaling Jacs : %.15e, %.15e \n", vcell->alpha[0], vcell->alpha[1]);
				//printf("Intersection ref elem variables : %.12e, %.12e, %.12e, %.12e\n",vcell->intersection_compvar[0],vcell->intersection_compvar[1],vcell->intersection_compvar[2],vcell->intersection_compvar[3]);
				//}
 //vcell->theFieldsCoefficients->get(j,3) = temp_c3;
        if (isLimited1 || isLimited2 || isLimited3 || isLimited4)//(flag)
		{   
			vcell->limit = 3;
			//if( -0.15 < centroid(0) && centroid(0) < -0.14 &&  0.16 < centroid(1) && centroid(1) < 0.167){;}
			//else if( -0.167 < centroid(0) && centroid(0) < -0.155 &&  -0.175 < centroid(1) && centroid(1) < -0.167){;}
			//else {vcell->theFieldsCoefficients->get(j,3) = 1.*temp_c3;}    
			vcell->theFieldsCoefficients->get(j,3) = 1.*temp_c3;	// modify c3 coefficient
			//printf("C3 coeffs : %.12e, %.12e\n", temp_c3,vcell->theFieldsCoefficients->getValCopy(j,3));

			// Reconstructing C2
			temp_c2 = vcell->theFieldsCoefficients->getValCopy(j,2);
			temp_c21 = vcell->theFieldsCoefficients->getValCopy(j,2);
			if(vcell->Neighbours[4] && vcell->Neighbours[10]){
			c31 = Mean_neigh[4];
			c32 = Mean_neigh[10];
			c3_f = vcell->neighbour_weights_cm[4]* c31 + vcell->neighbour_weights_cm[5]* c32; 
			temp_l21 = l2 * (2.*vcell->scaling_derivative[1]/(sqrt(3.) *  vcell->dist_intersection[2])) * (c3_f - mean );//(2./sqrt(3.))* 0.75 * (c3_f - mean );
			//temp_c2 = minmod(temp_c2, l2 * (2.*vcell->scaling_derivative[1]/(sqrt(3.) *  vcell->dist_intersection[2])) * (c3_f - mean ), isLimited3) ;
			temp_c2 = minmod(temp_c2, (2./sqrt(3.))* scalingfactor * (c3_f - mean ), isLimited3) ;
			}
			if(vcell->Neighbours[6] && vcell->Neighbours[11]){
				c31 = Mean_neigh[6];
				c32 = Mean_neigh[11];
				c3_b = vcell->neighbour_weights_cm[6]* c31 + vcell->neighbour_weights_cm[7]* c32;  
				//temp_c2 = minmod(temp_c2, l2 * (2.*vcell->scaling_derivative[1] /(sqrt(3.) * vcell->dist_intersection[3])) * (-c3_b  + mean), isLimited4) ;
			temp_c2 = minmod(temp_c2, (2./sqrt(3.))* scalingfactor * (-c3_b  + mean), isLimited4) ;
			
			}

			/*if(vcell->Neighbours[4] && vcell->Neighbours[5]){
			c31 = C0_neigh[4];
			c32 = C0_neigh[5];
			c3_f = vcell->neighbour_weights[4]* c31 + vcell->neighbour_weights[5]* c32; 
			temp_l22 = l2 * (1.*vcell->scaling_derivative[1]/(sqrt(3.) *  vcell->dist_intersection[2])) * (c3_f - vcell->theFieldsCoefficients->get(j,0) );
			Temp[0] = max(max(fabs(C1_neigh[4]),fabs(C2_neigh[4])),max(fabs(C1_neigh[5]),fabs(C2_neigh[5])));
			temp_c21 = minmod(temp_c21, l2 * (1.*vcell->scaling_derivative[1]/(sqrt(3.) *  vcell->dist_intersection[2])) * (c3_f - vcell->theFieldsCoefficients->get(j,0) ), isLimited3) ;
			}
			if(vcell->Neighbours[6] && vcell->Neighbours[7]){
				c31 = C0_neigh[6];
				c32 = C0_neigh[7];
				c3_b = vcell->neighbour_weights[6]* c31 + vcell->neighbour_weights[7]* c32;  
				temp_c21 = minmod(temp_c21, l2 * (1.*vcell->scaling_derivative[1] /(sqrt(3.) * vcell->dist_intersection[3])) * (-c3_b  +vcell->theFieldsCoefficients->get(j,0)), isLimited4) ;
			}*/

			//flag = (fabs(temp_c2 - vcell->theFieldsCoefficients->get(j,2))  > tolerance * fabs(vcell->theFieldsCoefficients->get(j,2)))?1:0;
			//if (flag) 
			//isLimited3 = 0;
			if (isLimited3 || isLimited4) {vcell->theFieldsCoefficients->get(j,2) = temp_c2; vcell->limit = 2;}
            
			
			// Reconstructing C1
			temp_c1 = vcell->theFieldsCoefficients->getValCopy(j,1);
			temp_c11 = vcell->theFieldsCoefficients->getValCopy(j,1);
			if(vcell->Neighbours[0] && vcell->Neighbours[8]){
			c31 = Mean_neigh[0];
			c32 = Mean_neigh[8];
			c3_f = vcell->neighbour_weights_cm[0]* c31 + vcell->neighbour_weights_cm[1]* c32; 
			temp_l11 = l1 * (2.*vcell->scaling_derivative[0] /(sqrt(3.) * vcell->dist_intersection[0]) * (c3_f - mean) );//(2./sqrt(3.))* 0.75 * (c3_f - mean);
			//temp_c1 = minmod(temp_c1, l1 * (2.*vcell->scaling_derivative[0] /(sqrt(3.) * vcell->dist_intersection[0]) * (c3_f - mean) ), isLimited1) ;
			temp_c1 = minmod(temp_c1, (2./sqrt(3.))* scalingfactor * (c3_f - mean), isLimited1) ;
			
			}
			if(vcell->Neighbours[2] && vcell->Neighbours[9]){
				c31 = Mean_neigh[2];
				c32 = Mean_neigh[9];
				c3_b = vcell->neighbour_weights_cm[2]* c31 + vcell->neighbour_weights_cm[3]* c32;  
				//temp_c1 = minmod(temp_c1, l1 * (2.*vcell->scaling_derivative[0] /(sqrt(3.) * vcell->dist_intersection[1]) * (-c3_b  + mean) ), isLimited2) ;
			temp_c1 = minmod(temp_c1, (2./sqrt(3.))* scalingfactor * (-c3_b  + mean), isLimited2) ;
		
			}
            /*if(vcell->Neighbours[0] && vcell->Neighbours[1]){
			c31 = C0_neigh[0];
			c32 = C0_neigh[1];
			c3_f = vcell->neighbour_weights[0]* c31 + vcell->neighbour_weights[1]* c32; 
			temp_l12 = l1 * (1.*vcell->scaling_derivative[0] /(sqrt(3.) * vcell->dist_intersection[0]) * (c3_f - vcell->theFieldsCoefficients->get(j,0)) );
			Temp[1] = max(max(fabs(C1_neigh[0]),fabs(C2_neigh[0])),max(fabs(C1_neigh[1]),fabs(C2_neigh[1])));
			
			temp_c11 = minmod(temp_c11, l1 * (1.*vcell->scaling_derivative[0] /(sqrt(3.) * vcell->dist_intersection[0]) * (c3_f - vcell->theFieldsCoefficients->get(j,0)) ), isLimited1) ;
			}
			if(vcell->Neighbours[2] && vcell->Neighbours[3]){
				c31 = C0_neigh[2];
				c32 = C0_neigh[3];
				c3_b = vcell->neighbour_weights[2]* c31 + vcell->neighbour_weights[3]* c32;  
				temp_c11 = minmod(temp_c11, l1 * (1.*vcell->scaling_derivative[0] /(sqrt(3.) * vcell->dist_intersection[1]) * (-c3_b  + vcell->theFieldsCoefficients->get(j,0)) ), isLimited2) ;
			}*/
			//flag = (fabs(temp_c2 - vcell->theFieldsCoefficients->get(j,2))  > tolerance * fabs(vcell->theFieldsCoefficients->get(j,2)))?1:0;
			//if (flag) 
			//isLimited2 = 0;
			if (isLimited1 || isLimited2) {vcell->theFieldsCoefficients->get(j,1) = temp_c1;  vcell->limit = 1;}
            

			double max_lincoeff = max(fabs(vcell->theFieldsCoefficients->getValCopy(j,1)),fabs(vcell->theFieldsCoefficients->getValCopy(j,2)));
			Temp[0] = max(Temp[0],Temp[1]);
			max_lincoeff = max(max_lincoeff,Temp[0]);
			vcell->LinearCoefferror[0] = 0.0;
			vcell->LinearCoefferror[1] = 0.0;
			/*if(max_lincoeff){
				vcell->LinearCoefferror[0] = max(fabs(temp_l21-temp_l22),fabs(temp_l11-temp_l12))/max_lincoeff;
				vcell->LinearCoefferror[1] = fabs(temp_c21-temp_c2)/max_lincoeff;
			}*/
			// 
			//else
			//if (fabs(vcell->theFieldsCoefficients->getValCopy(j,1)) < fabs(vcell->theFieldsCoefficients->getValCopy(j,2)))  
			if (fabs(vcell->theFieldsCoefficients->getValCopy(j,1)) > 1e-14)   vcell->LinearCoefferror[0] = fabs(vcell->theFieldsCoefficients->getValCopy(j,1) - temp_c1);///fabs(vcell->theFieldsCoefficients->getValCopy(j,1));
			if (fabs(vcell->theFieldsCoefficients->getValCopy(j,2)) > 1e-14) vcell->LinearCoefferror[1] = fabs(vcell->theFieldsCoefficients->getValCopy(j,2) - temp_c2);///fabs(vcell->theFieldsCoefficients->getValCopy(j,2));
			//if (fabs(vcell->theFieldsCoefficients->getValCopy(j,1))/vcell->scaling_derivative[0] < fabs(vcell->theFieldsCoefficients->getValCopy(j,2))/vcell->scaling_derivative[1])
            //    vcell->LinearCoefferror[0] = vcell->LinearCoefferror[1];
			//vcell->LinearCoefferror[0] = max(vcell->LinearCoefferror[0],vcell->LinearCoefferror[1]);
			vcell->LinearCoefferror[0] = fabs(vcell->alpha[0]) + fabs(vcell->alpha[1]);
			/*if( -0.2 < centroid(0) && centroid(0) < -0.15 && -.2 < centroid(1) && centroid(1) < -0.15){
				printf("C2 interpolation values : %.12e, %.12e, %.12e, %.12e, %.12e\n",C0_neigh[2], C0_neigh[3], c3_b,c3_f,vcell->theFieldsCoefficients->getValCopy(j,1));
				printf("C1 approx : %.12e, %.12e\n",(vcell->scaling_derivative[0] /(sqrt(3.) * vcell->dist_intersection[1]) * (-c3_b  + vcell->theFieldsCoefficients ->getCopy(j,0)) ),(vcell->scaling_derivative[0] /(sqrt(3.) * vcell->dist_intersection[0]) * (c3_f - vcell->theFieldsCoefficients ->getCopy(j,0)) ));
				printf("C0 : %.12e\n",vcell->theFieldsCoefficients ->getCopy(j,0));
				printf("Neighbour centroid : (%.12e,%.12e), %.12e\n",(vcell->Neighbours[2])->cellCentroid(0),(vcell->Neighbours[2])->cellCentroid(1),(vcell->Neighbours[2])->theFieldsCoefficients->getValCopy(j,0));
				printf("Neighbour centroid : (%.12e,%.12e), %.12e\n",(vcell->Neighbours[3])->cellCentroid(0),(vcell->Neighbours[3])->cellCentroid(1),(vcell->Neighbours[3])->theFieldsCoefficients->getValCopy(j,0));
		        printf("Weight coefficients : %.12e,%.12e\n",vcell->neighbour_weights[2],vcell->neighbour_weights[3] );
			}*/

			// Reconstructing C0
			// Error maintaining conservation! Roundoff error accumulating??
			double temp = 0.0;
			if (isLimited1 || isLimited2) {temp += (sqrt(3.)/3.) *  vcell->alpha[0] * (temp_c1 - vcell->theFieldsCoefficients->getValCopy(j,1)) ;}
			if (isLimited3 || isLimited4) {temp += (sqrt(3.)/3.) * (vcell->alpha[1] * (temp_c2 - vcell->theFieldsCoefficients->getValCopy(j,2)));} 
			vcell->theFieldsCoefficients->get(j,0) -= temp;
			//vcell->theFieldsCoefficients->get(j,0) = 2.0 * vcell->theMean[j] - (sqrt(3.)/3.) * (vcell->alpha[0] * temp_c1 +  vcell->alpha[1] * temp_c2 );
			//if (  (0.5*vcell->theFieldsCoefficients->get(j,0) + (sqrt(3.)/2.) * ( vcell->alpha[0] *(vcell->theFieldsCoefficients->get(j,1)) + vcell->alpha[1] * (vcell->theFieldsCoefficients->get(j,2))) ) < -0.5){
			//	printf("Coordinates : (%.12e,%.12e)\n",centroid(0),centroid(1));
			//}

			//double temp_mean = 0.0;
			//temp_mean = vcell->theFieldsCoefficients->get(j,0) * 0.5 + 0.5 * (sqrt(3.)/3.) *  vcell->alpha[0] * vcell->theFieldsCoefficients->get(j,1) + 0.5 * (sqrt(3.)/3.) *  vcell->alpha[1] * vcell->theFieldsCoefficients->get(j,2); 
		     //if (abs(temp_mean - mean) > 1e-10)
			//  {printf("Error in mean : %.12e\n", abs(temp_mean - mean));}
		}  
	}
  }
  return;
}
/******************************************************************************/

void DGVertexLimiter::limit(DGCell *vcell, double time)
{
  int fSize = vcell->fSize;
  int fOrder = vcell->theFunctionSpace->order();
  int n = vcell->theMeshEntity->getLevel();
  vcell->limit=5;
  int kk=0;
//s  double u,v,w,weight,val[MaxNbEqn], Max[MaxNbEqn],Min[MaxNbEqn];
//  double uvol,vvol,wvol,x,y,z;
  double phin;
  list<mEntity *> allSubs;
  list<mEntity*>::const_iterator it;
  //double dx = vcell->getMaxH();
  //  if (fabs(pow(vcell->fullJump,3)/vcell->maxVal/sqrt(pow(dx,3*(fOrder+1)))/dx)>1. )
  //vcell->limit = 1;
  //printf("Jump %f, JUmp^3 %f, maxVal %f, dx %f\n",vcell->fullJump,pow(vcell->fullJump,3),vcell->maxVal,dx);
  
if (vcell->limit)    
  {  
    double Max[MaxNbEqn],Min[MaxNbEqn];
    int n = vcell->theMeshEntity->getLevel();
    int k;
    int kk = 0;
    int fSize = vcell->fSize;
    int fOrder = vcell->theFunctionSpace->order();
    
    double val[MaxNbEqn];
    list<mEntity*>::const_iterator it;
    
    list<mEntity *> allVert;
    mEntity *ent;

	// Quad BJ limiter variables
	mPoint v[4], v_neigh[4];
	double coeffx[4], coeff_neighx[4], coeffy[4], coeff_neighy[4];
	double Jaccoeff[3], Jaccoeff_neigh[3];
	double Max_uder[MaxNbEqn], Min_uder[MaxNbEqn], Coef[3];
	double temp_der = 0, det_Jac = 0, temp_val[4] = {0.0,0.0,0.0,0.};
	int flag_limit[MaxNbEqn];
	int no_vertex = vcell->theMeshEntity->size(0);
	int no_edges = vcell->theMeshEntity->size(1);
	double tolerance = 1e-10, elem_tolerance = 1e-10;
	double norm_coeff3 = 1.0;

	for(k =0 ; k < no_edges; ++k){
		ent = vcell->theMeshEntity->get(1,k);
		if(ent->getClassification()->getId() == 510 || ent->getClassification()->getId() == 610){
			return;
		}
	}
    //printf("%d no of vertices \n", temp);
	//ent = vcell->theMeshEntity->get(0,2);
    for(k=0;k<no_vertex;k++)
      {
	ent = vcell->theMeshEntity->get(0,k); 
	v[k] = ((Vertex*)ent)->point();
	for(int i=0;i<ent->size(2);i++)
	  {
	    allVert.push_back(ent->get(2,i));
	    //ent->get(n,i)->print();
	  }
      }
	
	// Setting up the coefficients of determinant of Jacobian
	coeffx[0] = 0;
	coeffx[1] = 0.25*(v[1](0) + v[2](0) - (v[0](0) + v[3](0)));
	coeffx[2] = 0.25*(v[3](0) + v[2](0) - (v[0](0) + v[1](0)));
	coeffx[3] = 0.25*(v[0](0) + v[2](0) - (v[3](0) + v[1](0)));
	coeffy[0] = 0;
	coeffy[1] = 0.25*(v[1](1) + v[2](1) - (v[0](1) + v[3](1)));
	coeffy[2] = 0.25*(v[3](1) + v[2](1) - (v[0](1) + v[1](1)));
	coeffy[3] = 0.25*(v[0](1) + v[2](1) - (v[3](1) + v[1](1)));
	Jaccoeff[0] = coeffx[1]*coeffy[2] - coeffx[2]*coeffy[1];
	Jaccoeff[1] = coeffx[1]*coeffy[3] - coeffx[3]*coeffy[1];
	Jaccoeff[2] = coeffx[3]*coeffy[2] - coeffx[2]*coeffy[3];
	if (sqrt(pow(coeffx[3],2)) > tolerance || sqrt(pow(coeffx[3],2)) > tolerance){
		norm_coeff3 = sqrt(pow(coeffx[3],2) + pow(coeffx[3],2));
	}
    

    for(it = allVert.begin();it!=allVert.end();++it)
      {
		DGCell *cell = (DGCell*)(*it)->getCell();
         // Setting up detector for Quads
		/*for(k=0;k<vcell->theMeshEntity->size(0);k++)
		{
			v_neigh[k] = ((Vertex*)cell->theMeshEntity->get(0,k))->point();
		}

		// Setting up the coefficients of determinant of Jacobian
		coeff_neighx[0] = 0;
		coeff_neighx[1] = 0.25*(v_neigh[1](0) + v_neigh[2](0) - (v_neigh[0](0) + v_neigh[3](0)));
		coeff_neighx[2] = 0.25*(v_neigh[3](0) + v_neigh[2](0) - (v_neigh[0](0) + v_neigh[1](0)));
		coeff_neighx[3] = 0.25*(v_neigh[0](0) + v_neigh[2](0) - (v_neigh[3](0) + v_neigh[1](0)));
		coeff_neighy[0] = 0;
		coeff_neighy[1] = 0.25*(v_neigh[1](1) + v_neigh[2](1) - (v_neigh[0](1) + v_neigh[3](1)));
		coeff_neighy[2] = 0.25*(v_neigh[3](1) + v_neigh[2](1) - (v_neigh[0](1) + v_neigh[1](1)));
		coeff_neighy[3] = 0.25*(v_neigh[0](1) + v_neigh[2](1) - (v_neigh[3](1) + v_neigh[1](1)));
		Jaccoeff_neigh[0] = coeff_neighx[1]*coeff_neighy[2] - coeff_neighx[2]*coeff_neighy[1];
		Jaccoeff_neigh[1] = coeff_neighx[1]*coeff_neighy[3] - coeff_neighx[3]*coeff_neighy[1];
		Jaccoeff_neigh[2] = coeff_neighx[3]*coeff_neighy[2] - coeff_neighx[2]*coeff_neighy[3];*/
         
		int fSize_neigh1 = cell->theFunctionSpace->size(1);
	
	// if(cell->theMeshEntity->size() == 2)
	if(!kk)
	  {
	    for(int j=0;j<cell->cSize;j++){
		  Max[j] = Min[j] = cell->theMean[j];
		  //Coef[0] = cell->theFieldsCoefficients->getCopy(j,1);
		  //Coef[1] = cell->theFieldsCoefficients->getCopy(j,2);
		  //Max_uder[j] = Min_uder[j] =  (1./norm_coeff3)*(sqrt(3.)*0.5/Jaccoeff_neigh[0])*(coeffx[3]*(Coef[0]*coeff_neighy[2] - Coef[1]*coeff_neighy[1]) + coeffy[3]*(Coef[0]*coeff_neighx[2] - Coef[1]*coeff_neighx[1]));

		 }
	    }
	else
	  {
	    for(int j=0;j<cell->cSize;j++)
	      {
		Max[j]  = (cell->theMean[j] > Max[j])?cell->theMean[j]:Max[j];
		Min[j]  = (cell->theMean[j] < Min[j])?cell->theMean[j]:Min[j];

		//Coef[0] = cell->theFieldsCoefficients ->getCopy(j,1);
		//Coef[1] = cell->theFieldsCoefficients ->getCopy(j,2);
		//temp_der = (1./norm_coeff3)*(sqrt(3.)*0.5/Jaccoeff_neigh[0])*(coeffx[3]*(Coef[0]*coeff_neighy[2] - Coef[1]*coeff_neighy[1]) + coeffy[3]*(Coef[0]*coeff_neighx[2] - Coef[1]*coeff_neighx[1]));
        //Max_uder[j] = (temp_der > Max_uder[j])?temp_der:Max_uder[j];
		//Min_uder[j] = (temp_der < Min_uder[j])?temp_der:Min_uder[j]; 

		// printf("Min= %e Max=%e Mean =%e\n",Min[j],Max[j],cell->theMean[j]);
		//(*it)->print();
	      }
	    //if ( Min[0]<0 || Min[3] <0 ) printf("NEGATIVE MEAN  \n");
	  }
	kk++;
      }
  //printf("Min= %e Max=%e\n",Min[0],Max[0]);
  
  for(int j=0;j<vcell->cSize;j++)
    vcell->limSlope[j] = 1.0;
  
  for(k=0;k<vcell->theMeshEntity->size(n-1);k++)
    recurGetAllsubs(vcell->theMeshEntity->get(n-1,k),allSubs);
  /*for(it = allSubs.begin();it!=allSubs.end();++it)
    {
      DGBoundaryCell *cell = (DGBoundaryCell*)(*it)->getCell();      
      
      DGCell *right = 0;
      DGCell *left = (DGCell*)cell->mleft->getCell();
      if(cell->mright)right = (DGCell*)cell->mright->getCell();*/
      
      double V[4][3];
	  // Shock detector?? - Need corner points
     /* V[0][0] = -1.0; V[0][1] = -1.0; V[0][2] = 0.0;
      V[1][0] = 1.0; V[1][1] = -1.0; V[1][2] = 0.0;
      V[2][0] = 1.0; V[2][1] = 1.0; V[2][2] = 0.0;
	  V[3][0] = -1.0; V[3][1] = 1.0; V[3][2] = 0.0;
	  for(int i=0; i<4; ++i){
		  
		  for(int j=0;j<vcell->cSize;j++){
			  if(!i) flag_limit[j] = 0;
			  Coef[0] = vcell->theFieldsCoefficients ->get(j,1);
			  Coef[1] = vcell->theFieldsCoefficients ->get(j,2);
			  Coef[2] = vcell->theFieldsCoefficients ->get(j,3);

              det_Jac = Jaccoeff[0] + Jaccoeff[1]*V[i][0] + Jaccoeff[2]*V[i][1];
			  temp_der = (sqrt(3.0)*0.5)*(Coef[0]*Jaccoeff[2] + Coef[1]*Jaccoeff[1] - sqrt(3.)*Coef[2]*Jaccoeff[0])/det_Jac + 1.5*Coef[2];
			  temp_der = (1./norm_coeff3)*temp_der;
			  // Check for linear solution
			  //if(sqrt(pow(Coef[0]*Jaccoeff[2] + Coef[1]*Jaccoeff[1] - sqrt(3.)*Coef[2]*Jaccoeff[0],2)) < tolerance){
				//  flag_limit[j] = 1;
				//  break;
			  //}
              // Check for near cartesian elements?
			  if(sqrt(pow(coeffx[3],2)) < elem_tolerance && sqrt(pow(coeffy[3],2)) < elem_tolerance){
				  flag_limit[j] = 1;
				  break;
			  }
			  if((temp_der - Max_uder[j]) > tolerance || (Min_uder[j] - temp_der) > tolerance){
				  flag_limit[j] = 1;
				  break;
			  }
		  }
	  }*/

	  // vertex neighbourhood limiter
      V[0][0] = 0.0; V[0][1] = -1.0; V[0][2] = 0.0;
      V[1][0] = 1.0; V[1][1] = 0.0; V[1][2] = 0.0;
      V[2][0] = 0.0; V[2][1] = 1.0; V[2][2] = 0.0;
	  V[3][0] = -1.0; V[3][1] = 0.0; V[3][2] = 0.0;
	  /*V[0][0] = -1.0; V[0][1] = -1.0; V[0][2] = 0.0;
      V[1][0] = 1.0; V[1][1] = -1.0; V[1][2] = 0.0;
      V[2][0] = 1.0; V[2][1] = 1.0; V[2][2] = 0.0;
	  V[3][0] = -1.0; V[3][1] = 1.0; V[3][2] = 0.0;*/
	  /*V[0][0] = 0.5; V[0][1] = 0.0; V[0][2] = 0.0;
      V[1][0] = 0.5; V[1][1] = 0.5; V[1][2] = 0.0;
      V[2][0] = 0.0; V[2][1] = 0.5; V[2][2] = 0.0;
	  V[3][0] = 0.0; V[3][1] = 0.5; V[3][2] = 0.0;*/
      for(int i=0;i<4;i++){
		  vcell->interpolate(V[i][0],V[i][1],V[i][2],val);
		for(int j=0;j<vcell->cSize;j++)
		{   
			//if(flag_limit[j])
			{
				
				double x;
				if(true)
				//	    if(fabs(val[j] - vcell->theMean[j]) > 1.e-8 * fabs(val[j] + vcell->theMean[j]) &&fabs( vcell->theMean[j]) > 1.e-13 )	    
				{
				if(val[j]  > Max[j])
				{
					x = Max[j];
					phin = (x-vcell->theMean[j])/(val[j]-vcell->theMean[j]);
					vcell->limit=3;
				}
				else if (val[j]  < Min[j])
				{
					x = Min[j];
					phin = (x-vcell->theMean[j])/(val[j]-vcell->theMean[j]);
					vcell->limit=3;
				}
				else  phin = 1.0;
				}
				else       phin = 1.0;

				if(phin<0.0) {phin = 0.0;vcell->limit=2;}
				vcell->limSlope[j] = (vcell->limSlope[j]<phin)?vcell->limSlope[j]:phin;
			}
		}
	  }
    //}   

  while(1)
    { 
      int fSize1 = vcell->theFunctionSpace->size(1);
      int fSizei = vcell->theFunctionSpace->size(fOrder);
      int fSizej = vcell->theFunctionSpace->size(fOrder-1);
      
      for(int j=0;j<vcell->cSize;j++)
	{	
	  if(fOrder == 1){
		  /*if(vcell->limit < 5){
			  temp_der = vcell->theFieldsCoefficients->get(j,1)*Jaccoeff[2]/Jaccoeff[0] + vcell->theFieldsCoefficients->get(j,2)*Jaccoeff[1]/Jaccoeff[0] - sqrt(3.0)*vcell->theFieldsCoefficients->get(j,3); 
		  	  temp_val[0] = 0.5*sqrt(3.)*temp_der/(1.0 + Jaccoeff[1]/Jaccoeff[0]) + 1.5*vcell->theFieldsCoefficients->get(j,3);
			  temp_val[1] = 0.5*sqrt(3.)*temp_der/(1.0 + Jaccoeff[2]/Jaccoeff[0]) + 1.5*vcell->theFieldsCoefficients->get(j,3);
			  temp_val[2] = 0.5*sqrt(3.)*temp_der/(1.0 - Jaccoeff[1]/Jaccoeff[0]) + 1.5*vcell->theFieldsCoefficients->get(j,3);
			  temp_val[3] = 0.5*sqrt(3.)*temp_der/(1.0 - Jaccoeff[2]/Jaccoeff[0]) + 1.5*vcell->theFieldsCoefficients->get(j,3);
			  //printf("Gradients at the midpoints at : %1.2e, %.12e, %.12e, %.12e \n ", temp_val[0], temp_val[1], temp_val[2], temp_val[3]);
			  //printf("Curvature : %.12e, alpha : %.12e,  beta : %.12e \n",temp_der,Jaccoeff[1]/Jaccoeff[0],Jaccoeff[2]/Jaccoeff[0]);
			  printf("Curvature : %.12e\n",temp_der);
		  }*/
		//temp_der = sqrt(pow(Jaccoeff[2]/Jaccoeff[0],2)) + sqrt(pow(Jaccoeff[1]/Jaccoeff[0],2));
		//if(temp_der > 0.5) printf("Curvature : %.12e, alpha : %.12e,  beta : %.12e \n",temp_der,Jaccoeff[1]/Jaccoeff[0],Jaccoeff[2]/Jaccoeff[0]);

		 // if(vcell->limSlope[j] < 1.0) printf("Limited slope : %f, %d\n", vcell->limSlope[j], fSize1);
		  vcell->theFieldsCoefficients ->get(j,0) = 2.0*(1.-vcell->limSlope[j])*vcell->theMean[j] + vcell->limSlope[j]*vcell->theFieldsCoefficients ->get(j,0);
			for (int k=1;k<fSize1;k++)
	     		vcell->theFieldsCoefficients->get(j,k) *= vcell->limSlope[j];
	  }
	    
	  else	
	    for (int k=fSizej;k<fSizei;k++)
	      if(vcell->limSlope[j] < 0.9)vcell->theFieldsCoefficients ->get(j,k)= 0.0;
	}	 
      if(fOrder == 1) break;
      fOrder --;
    } 
  /*
///MODIFY!!!
    int order = cell->computeOrder();
    order=1;
      GaussIntegrator gauss(cell->theBoundaryEntity);
      for(int i=0;i<gauss.nbIntegrationPoints(order);i++)
	{
	  gauss.iPoint(i,order,u,v,w,weight);
	  cell->theMapping->eval(u,v,w,x,y,z);
	  vcell->theMapping->invert(x,y,z,uvol,vvol,wvol);
	  vcell->interpolate(uvol,vvol,wvol,val);
	  // printf("u=%f v=%f, x=%f y=%f,uvol=%f vvol=%f\n",u,v,x,y,uvol,vvol);
	  for(int j=0;j<vcell->cSize;j++)
	    {
	      double x;
	      if(fabs(val[j] - vcell->theMean[j]) > 1.e-5 * fabs(val[j] + vcell->theMean[j]) &&fabs( vcell->theMean[j]) > 1.e-13 )	    
		{
		  if(val[j] > Max[j])
		    {
		      x = Max[j];
		      phin = (x-vcell->theMean[j])/(val[j]-vcell->theMean[j]);
		      vcell->limit=3;
		      //    printf("val>Max val[%d]=%e Max=%e Mean=%e\n",j,val[j],Max[j],vcell->theMean[j]);
		    }
		  else if(val[j] < Min[j])
		    {
		      x = Min[j];
		      phin = (x-vcell->theMean[j])/(val[j]-vcell->theMean[j]);
		      vcell->limit=3;
		      //printf("val<Min val[%d]=%e Min=%e Mean=%e\n",j,val[j],Min[j],vcell->theMean[j]);
		    }
		  else 
		    {
		      phin = 1.0;
		    }
		}
	      else
		phin = 1.0;
	      if(phin<0.0) {phin = 0.0;vcell->limit=2;}
	      vcell->limSlope[j] = (vcell->limSlope[j]<phin)?vcell->limSlope[j]:phin;
	      }
	      }*/
  }
} 

/************************************************/
void VertLimiter::limit(DGCell *vcell, double time)
{
  const int fSize = vcell->fSize;
  int fOrder = vcell->theFunctionSpace->order();
  const int cSize = vcell->cSize;
  const int n = vcell->theMeshEntity->getLevel();
  vcell->limit=1;
  int  kk=0;
  //double u,v,w,weight;
  //double uvol,vvol,wvol,x,y,z,phin;
  double phin,val[MaxNbEqn], Max[MaxNbEqn],Min[MaxNbEqn];
  list<mEntity*>::const_iterator it,itt;
  list<mEntity*>::const_iterator vertEnd=vcell->allVert.end(); 
  //double dx = vcell->getMaxH();
  //  if (fabs(pow(vcell->fullJump,3)/vcell->maxVal/sqrt(pow(dx,3*(fOrder+1)))/dx)>1. )
  //vcell->limit = 1;
  //printf("Jump %f, JUmp^3 %f, maxVal %f, dx %f\n",vcell->fullJump,pow(vcell->fullJump,3),vcell->maxVal,dx);
  
  if (vcell->limit) {  
    for(it = vcell->allVert.begin();it!=vertEnd;++it)
      {
	DGCell *cell = (DGCell*)(*it)->getCell();
	
	if(!kk)
	  {
	    for(int j=0;j<cSize;j++)
	      Max[j] = Min[j] = cell->theMean[j];
	  }
	else
	    for(int j=0;j<cSize;j++)
	      {
		Max[j]  = (cell->theMean[j] > Max[j])?cell->theMean[j]:Max[j];
		Min[j]  = (cell->theMean[j] < Min[j])?cell->theMean[j]:Min[j];
		//printf("Min= %e Max=%e Mean =%e\n",Min[j],Max[j],cell->theMean[j]);
		//(*it)->print();
	      }
	//if ( Min[0]<0 || Min[3] <0 ) printf("NEGATIVE MEAN  \n");
	
	kk++;
      }
    
    for(int j=0;j<cSize;j++)
      vcell->limSlope[j] = 1.0;
         boundaryGaussPoint *pg;
    list<mEntity*>::const_iterator endE = vcell->allSubs.end();
    for(it = vcell->allSubs.begin();it!=endE;++it)
      { if ((*it)->getClassification()->getId()!=20000) {
	DGBoundaryCell *cell = (DGBoundaryCell*)(*it)->getCell();  
	//one-two poit switch starts here    
	/*int order=1;
	GaussIntegrator gauss;
	int Nb = gauss.nbIntegrationPoints( cell->theBoundaryEntity,order);
	for(int i=0;i<Nb;i++)
	  {
	    gauss.iPoint(cell->theBoundaryEntity,i,order,u,v,w,weight);
	    cell->theMapping->eval(u,v,w,x,y,z);
	    vcell->theMapping->invert(x,y,z,uvol,vvol,wvol);*/
	// and ends here
	int Nb = cell->getNbPtGauss();
	for(int i=0;i<Nb;i++)
	  {
	    pg = cell->pt(i);
	    vcell->interpolate(pg->fctleft, val);
	    //vcell->interpolate(uvol,vvol,wvol,val);
	    
	    for(int j=0;j<cSize;j++)
	      {
		double x;
		//if(fabs(val[j] - vcell->theMean[j]) > 1.e-5 * fabs(val[j] + vcell->theMean[j]) &&fabs( vcell->theMean[j]) > 1.e-13 )	
		if (1)//(pg->x>0.5)
		  {
		    if(val[j] > Max[j])
		      {
		      x = Max[j];
		      phin = (x-vcell->theMean[j])/(val[j]-vcell->theMean[j]);
		      vcell->limit=2;
		       if (j==0)  printf("val>Max val[%d]=%e Max=%e Mean=%e\n",j,val[j],Max[j],vcell->theMean[j]);
		      }
		  else if(val[j] < Min[j])
		    {
		      x = Min[j];
		      phin = (x-vcell->theMean[j])/(val[j]-vcell->theMean[j]);
		      vcell->limit=2;
		     if (j==0) printf("val<Min val[%d]=%e Min=%e Mean=%e\n",j,val[j],Min[j],vcell->theMean[j]);
		    }
		    else 
		      {
			phin = 1.0;vcell->limit=3;
		      }
		  }
		else
		  phin = 1.0;
		if(phin<0.0) {phin = 0.0;vcell->limit=5;}
		vcell->limSlope[j] = (vcell->limSlope[j]<phin)?vcell->limSlope[j]:phin;
		if (vcell->limSlope[j]<1.) {
		  for(itt = vcell->allVert.begin();itt!=vertEnd;++itt){
		    DGCell *ccell = (DGCell*)(*itt)->getCell();
	if (j==0)	    printf(" %e ",ccell->theMean[j]);}
		 if (j==0) printf("val [%d] =%e MIN[j] =%e MAx[j] =%e vcell->limSlope[j]=%e \n",j,val[j],Min[j],Max[j],vcell->limSlope[j]);
		  if (j==0) printf(" %e %e\n",pg->x,pg->y);
		  printf("\n");
		  }
	      }
	  }
      }
      }
    while(1)
      { 
	int fSize1 = vcell->theFunctionSpace->size(1);
	int fSizei = vcell->theFunctionSpace->size(fOrder);
      int fSizej = vcell->theFunctionSpace->size(fOrder-1);
      
      for(int j=0;j<cSize;j++)
	{	
	  if(fOrder == 1)
	    for (int k=1;k<fSize1;k++)
	      vcell->theFieldsCoefficients ->get(j,k)*= vcell->limSlope[j];
	  else	
	    for (int k=fSizej;k<fSizei;k++)
	      if(vcell->limSlope[j] < 0.9)vcell->theFieldsCoefficients ->get(j,k)= 0.0;
	}	 
      if(fOrder == 1) break;
      fOrder --;
      } 
  }
} 


/**************************************************************/

void DGLimiter::computeMinMaxEdge(DGCell *vcell)
{
  double Max[MaxNbEqn],Min[MaxNbEqn];
  int n = vcell->theMeshEntity->getLevel();
  int k;
  int kk = 0;
  int fSize = vcell->fSize;
  int fOrder = vcell->theFunctionSpace->order();
  
  double u,v,w,weight,val[MaxNbEqn];
  double uvol,vvol,wvol,x,y,z,phin;
  
  

  for(int i=0;i<vcell->theMeshEntity->size(n-1);i++)
    {
      mEntity *bound = vcell->theMeshEntity->get(n-1,i);
      if(bound->size(n) == 2)
	{
	  DGCell *other;
	  if(bound->get(n,0) != vcell->theMeshEntity)
	    other = (DGCell*)bound->get(n,0)->getCell();
	  else
	    other = (DGCell*)bound->get(n,1)->getCell();
	  if(!kk)
	    {
	      for(int j=0;j<other->cSize;j++)
		Max[j] = Min[j] = other->theMean[j];
	    }
	  else
	    {
	      for(int j=0;j<other->cSize;j++)
		{
		  Max[j]  = (other->theMean[j] > Max[j])?other->theMean[j]:Max[j];
		  Min[j]  = (other->theMean[j] < Min[j])?other->theMean[j]:Min[j];
		  //  printf("Min= %e Max=%e\n",Min[j],Max[j]);
		}
	      // if ( Min[0]<0 || Min[3] <0 ) printf("NEGATIVE MEAN  \n");
	    }
	  kk++;
	}
    }
  
  //printf("Min= %e Max=%e\n",Min[0],Max[0]);

  list<mEntity *> allSubs;
  list<mEntity*>::const_iterator it;
  for(k=0;k<vcell->theMeshEntity->size(n-1);k++)
    recurGetAllsubs(vcell->theMeshEntity->get(n-1,k),allSubs);

  for(it = allSubs.begin();it!=allSubs.end();++it)
    {
      DGBoundaryCell *cell = (DGBoundaryCell*)(*it)->getCell();      
      ///DOn't need this!
      DGCell *right = 0;
      DGCell *left = (DGCell*)cell->mleft->getCell();
      if(cell->mright)right = (DGCell*)cell->mright->getCell();  
  for(int j=0;j<vcell->cSize;j++)
    vcell->limSlope[j] = 1.0;
  
  int order = cell->computeOrder();
  order=1;
  GaussIntegrator gauss;
  int Nb = gauss.nbIntegrationPoints(cell->theBoundaryEntity,order);
  for(int i=0;i<Nb;i++)
    {
      gauss.iPoint(cell->theBoundaryEntity,i,order,u,v,w,weight);
      cell->theMapping->eval(u,v,w,x,y,z);
      vcell->theMapping->invert(x,y,z,uvol,vvol,wvol);
      vcell->interpolate(uvol,vvol,wvol,val);
      // printf("u=%f v=%f, x=%f y=%f,uvol=%f vvol=%f\n",u,v,x,y,uvol,vvol);
      for(int j=0;j<vcell->cSize;j++)
	{
	  double x;
	  if(fabs(val[j] - vcell->theMean[j]) > 1.e-5 * fabs(val[j] + vcell->theMean[j]) &&fabs( vcell->theMean[j]) > 1.e-13 )	    
	    {
	      if(val[j] > Max[j])
		{
		  x = Max[j];
		  phin = (x-vcell->theMean[j])/(val[j]-vcell->theMean[j]);
		  vcell->limit=3;
		  //    printf("val>Max val[%d]=%e Max=%e Mean=%e\n",j,val[j],Max[j],vcell->theMean[j]);
		}
	      else if(val[j] < Min[j])
		{
		  x = Min[j];
		  phin = (x-vcell->theMean[j])/(val[j]-vcell->theMean[j]);
		  vcell->limit=3;
		  //printf("val<Min val[%d]=%e Min=%e Mean=%e\n",j,val[j],Min[j],vcell->theMean[j]);
		}
	      else 
		{
		  phin = 1.0;
		}
	    }
	  else
	    phin = 1.0;
	  if(phin<0.0) {phin = 0.0;vcell->limit=2;}
	  vcell->limSlope[j] = (vcell->limSlope[j]<phin)?vcell->limSlope[j]:phin;
	}
    }
}
}

/*********************************************************************************/


void DGMomentLimiter::limit(DGCell *vcell,double time)
{
  int d,k,i,q;
  int k_im1,km1_i,i_km1,im1_k;
  DGCell *left,*right,*top,*bottom;
  left=right=top=bottom=0;
  mPoint p1 = vcell->pMin;
  mPoint p2 = vcell->pMax;
  mPoint p=(p1+p2);
  p(0)/=2.;p(1)/=2.;p(2)/=2.;
  mPoint omin,omax;
  const int cSize = vcell->cSize;
  const int n = vcell->theMeshEntity->getLevel();
  const double small_number = 1.0e-10;
  for(i=0;i<vcell->theMeshEntity->size(n-1);i++)
    {
      mEntity *bound = vcell->theMeshEntity->get(n-1,i);
      DGCell *other=0;
      if(bound->get(n,0) != vcell->theMeshEntity)
	other = (DGCell*)bound->get(n,0)->getCell();
        else if (  bound->get(n,1)) other = (DGCell*)bound->get(n,1)->getCell();
      if (other)
	{
	  omin = other -> pMin;
	  omax = other -> pMax;
	  if (p(0) > omax(0)) left = other; 
	  if (p(0) < omin(0)) right = other; 
	  if (p(1) > omax(1)) bottom = other; 
	  if (p(1) < omin(1)) top = other; 
	}	 
    }
   
  double tmp1,tmp2;
  bool isLimited;
  bool stopLimiter = 0;
  
  vcell->theFieldsCoefficients->copy();

  for(q=0; q<cSize; q++) {
    d=(vcell->fOrder+1)*(vcell->fOrder+1)-1;
    vcell->limit = d;
    for (k=vcell->fOrder; k>0; k--)
      {
	if (stopLimiter) break;
	
	if (fabs(vcell->theFieldsCoefficients->getCopy(q,d))>small_number) 
	  {
	    if (top && bottom)
	      tmp2 = minmod(vcell->theFieldsCoefficients->getCopy(q,d),
			    (top->theFieldsCoefficients->get(q,d-2)-vcell->theFieldsCoefficients->get(q,d-2))/sqrt(double (2*k+1)/(2*k-1)),
			    (vcell->theFieldsCoefficients->get(q,d-2)-bottom->theFieldsCoefficients->get(q,d-2))/sqrt(double (2*k+1)/(2*k-1)),isLimited);
	    else if (top && !bottom) 
	      tmp2 = minmod(vcell->theFieldsCoefficients->getCopy(q,d),
			    (top->theFieldsCoefficients->get(q,d-2)-vcell->theFieldsCoefficients->get(q,d-2))/sqrt(double (2*k+1)/(2*k-1)), isLimited);
	    else if (!top&&bottom)
	      tmp2 = minmod(vcell->theFieldsCoefficients->getCopy(q,d),
			    (vcell->theFieldsCoefficients->get(q,d-2)-bottom->theFieldsCoefficients->get(q,d-2))/sqrt(double (2*k+1)/(2*k-1)),isLimited);
	    
	    bool flag = isLimited;
	    
	    if (left && right)
	      tmp1 = minmod(vcell->theFieldsCoefficients->getCopy(q,d),
			    (right->theFieldsCoefficients->get(q,d-1)-vcell->theFieldsCoefficients->get(q,d-1))/sqrt(double (2*k+1)/(2*k-1)),
			    (vcell->theFieldsCoefficients->get(q,d-1)-left->theFieldsCoefficients->get(q,d-1))/sqrt(double (2*k+1)/(2*k-1)),isLimited);
	    else if (left && !right)
	      tmp1 = minmod(vcell->theFieldsCoefficients->getCopy(q,d),
			    (vcell->theFieldsCoefficients->get(q,d-1)-left->theFieldsCoefficients->get(q,d-1))/sqrt(double (2*k+1)/(2*k-1)),isLimited);
	    else if (!left && right)
	      tmp1 = minmod(vcell->theFieldsCoefficients->getCopy(q,d),
			    (right->theFieldsCoefficients->get(q,d-1)-vcell->theFieldsCoefficients->get(q,d-1))/sqrt(double (2*k+1)/(2*k-1)), isLimited);
	    
	    if (flag || isLimited) {
	      vcell->limit=d-1;
	      vcell->theFieldsCoefficients->getCopy(q,d) = absmin(tmp1,tmp2);
	      //vcell->getMeshEntity()->print();
	    } 
	    else {
	      stopLimiter=1;
	      break;
	    }
	  }
	//	else printf("A coeff smaller than small number\n");
	d--;
	
	for (i=k-1; i>=0; i--) 
	  {
	    if (fabs(vcell->theFieldsCoefficients->getCopy(q,d))>small_number)
	      {
		double tmp1,tmp2;
		if (i<k-1) i_km1 =(k-1)*(k-1)+2*i+1; else i_km1 = i*i+(k-1)*2; 
		if (i>0) im1_k = k*k+2*(i-1)+1;
		
		if (top && bottom)
		  tmp2 = minmod(vcell->theFieldsCoefficients->getCopy(q,d),
				(top->theFieldsCoefficients->get(q,i_km1)-vcell->theFieldsCoefficients->get(q,i_km1))/sqrt(double (2*k+1)/(2*k-1)),
				(vcell->theFieldsCoefficients->get(q,i_km1)-bottom->theFieldsCoefficients->get(q,i_km1))/sqrt(double (2*k+1)/(2*k-1)),isLimited);
		else if (top && !bottom) 
		  tmp2 = minmod(vcell->theFieldsCoefficients->getCopy(q,d),
				(top->theFieldsCoefficients->get(q,i_km1)-vcell->theFieldsCoefficients->get(q,i_km1))/sqrt(double (2*k+1)/(2*k-1)), isLimited);
		else if (!top&&bottom)
		  tmp2 = minmod(vcell->theFieldsCoefficients->getCopy(q,d),
				(vcell->theFieldsCoefficients->get(q,i_km1)-bottom->theFieldsCoefficients->get(q,i_km1))/sqrt(double (2*k+1)/(2*k-1)),isLimited);
		
		bool flag = isLimited;
		
		if (left && right&&i)
		  tmp1 = minmod(vcell->theFieldsCoefficients->getCopy(q,d),
				(right->theFieldsCoefficients->get(q,im1_k)-vcell->theFieldsCoefficients->get(q,im1_k))/sqrt(double (2*i+1)/(2*i-1)),
				(vcell->theFieldsCoefficients->get(q,im1_k)-left->theFieldsCoefficients->get(q,im1_k))/sqrt(double (2*i+1)/(2*i-1)),isLimited);
		else if (left && !right&&i)
		  tmp1 = minmod(vcell->theFieldsCoefficients->getCopy(q,d),
				(vcell->theFieldsCoefficients->get(q,im1_k)-left->theFieldsCoefficients->get(q,im1_k))/sqrt(double (2*i+1)/(2*i-1)),isLimited);
		else if (!left && right&&i)
		  tmp1 = minmod(vcell->theFieldsCoefficients->getCopy(q,d),
				(right->theFieldsCoefficients->get(q,im1_k)-vcell->theFieldsCoefficients->get(q,im1_k))/sqrt(double (2*i+1)/(2*i-1)), isLimited);
		
		if (flag || isLimited) {
		  vcell->limit=d-1;
		  if (i) vcell->theFieldsCoefficients->getCopy(q,d) = absmin(tmp1,tmp2);
		  else vcell->theFieldsCoefficients->getCopy(q,d) = tmp2;
		} 
	      }
	    d--;
	    
	    if (fabs(vcell->theFieldsCoefficients->getCopy(q,d))>small_number)
	      {
		isLimited = 0;
		k_im1 = k*k+2*(i-1);
		km1_i = (k-1)*(k-1)+2*i;
		if (top && bottom &&i)
		  tmp2 = minmod(vcell->theFieldsCoefficients->getCopy(q,d),
				(top->theFieldsCoefficients->get(q,k_im1)-vcell->theFieldsCoefficients->get(q,k_im1))/sqrt(double (2*i+1)/(2*i-1)),
				(vcell->theFieldsCoefficients->get(q,k_im1)-bottom->theFieldsCoefficients->get(q,k_im1))/sqrt(double (2*i+1)/(2*i-1)),isLimited);
		else if (top && !bottom && i) 
		  tmp2 = minmod(vcell->theFieldsCoefficients->getCopy(q,d),
				(top->theFieldsCoefficients->get(q,k_im1)-vcell->theFieldsCoefficients->get(q,k_im1))/sqrt(double (2*i+1)/(2*i-1)), isLimited);
		else if (!top&&bottom && i)
		  tmp2 = minmod(vcell->theFieldsCoefficients->getCopy(q,d),
				(vcell->theFieldsCoefficients->get(q,k_im1)-bottom->theFieldsCoefficients->get(q,k_im1))/sqrt(double (2*i+1)/(2*i-1)),isLimited);
		
		bool flag = isLimited;
		
		if (left && right)
		  tmp1 = minmod(vcell->theFieldsCoefficients->getCopy(q,d),
				(right->theFieldsCoefficients->get(q,km1_i)-vcell->theFieldsCoefficients->get(q,km1_i))/sqrt(double (2*k+1)/(2*k-1)),
				(vcell->theFieldsCoefficients->get(q,km1_i)-left->theFieldsCoefficients->get(q,km1_i))/sqrt(double (2*k+1)/(2*k-1)),isLimited);
		else if (left && !right)
		  tmp1 = minmod(vcell->theFieldsCoefficients->getCopy(q,d),
				(vcell->theFieldsCoefficients->get(q,km1_i)-left->theFieldsCoefficients->get(q,km1_i))/sqrt(double (2*k+1)/(2*k-1)),isLimited);
		else if (!left && right)
		  tmp1 = minmod(vcell->theFieldsCoefficients->getCopy(q,d),
				(right->theFieldsCoefficients->get(q,km1_i)-vcell->theFieldsCoefficients->get(q,km1_i))/sqrt(double (2*k+1)/(2*k-1)), isLimited);
		
		if (flag || isLimited) {
		  vcell->limit=d-1;
		  if (i) vcell->theFieldsCoefficients->getCopy(q,d) = absmin(tmp1,tmp2);
		  else vcell->theFieldsCoefficients->getCopy(q,d) = tmp1;
		} 
		else {
		  if (vcell->limit != d) 
		    {
		      stopLimiter=1;
		      break;
		    }
		}
	      }
	    d--;
	  }
      }
  }
}

inline double dotprod(int q,double *ev,double *x)
{	
  int cSize = 4;
  double sum =0.0;
  for(int s=0; s<4; s++) 
	  sum +=ev[q*cSize+s]*x[s];
  return sum;
}

inline double rightprod(int q,double *ev,double *x)
{
  int cSize = 4;
  double sum =0.0;
  for(int s=0; s<4; s++) sum +=ev[s*cSize+q]*x[s];
  return sum;
}

////*********************LIMITER FOR EULER EQUATIONS***********************************/

void DGMomentLimiterEuler::limit(DGCell *vcell,double time)
{
  const double small_number = 1.0e-05;
  int d,k,i,q;
  int k_im1,km1_i,i_km1,im1_k;
  DGCell *left,*right,*top,*bottom;
  left=right=top=bottom=0;
  mPoint p1 = vcell->pMin;
  mPoint p2 = vcell->pMax;
  mPoint p=(p1+p2);
  p(0)/=2.;p(1)/=2.;p(2)/=2.;
  mPoint omin,omax;
  const int cSize = vcell->cSize;
  const int n = vcell->theMeshEntity->getLevel();
  for(i=0;i<vcell->theMeshEntity->size(n-1);i++)
    {
      mEntity *bound = vcell->theMeshEntity->get(n-1,i);
      DGCell *other=0;
      if(bound->get(n,0) != vcell->theMeshEntity)
	other = (DGCell*)bound->get(n,0)->getCell();
      else if (  bound->get(n,1)) other = (DGCell*)bound->get(n,1)->getCell();
      if (other)
	{
	  omin = other -> pMin;
	  omax = other -> pMax;
	  
	  if (p(0) > omax(0)) if (bound->getClassification()->getId() != 510) left = other; else right =other;
	  if (p(0) < omin(0)) if (bound->getClassification()->getId() != 510) right = other; else left =other;
	  if (p(1) > omax(1)) if (bound->getClassification()->getId() != 610) bottom = other; else top =other;
	  if (p(1) < omin(1)) if (bound->getClassification()->getId() != 610) top = other; else bottom =other;
	}	 
    }
  
  double tmp[4];
  bool isLimited[4], wasLimited[4], continue_limiting[4], stopLimiter,previous_coeff_limited[4];

  stopLimiter = 0;
  for (q=0;q<4;q++) continue_limiting[q] = 1;
  
  vcell->theFieldsCoefficients->copy();
 
  double lev_x[16],lev_y[16],rev_x[16],rev_y[16];
  mVector pp1(1,0,0);
  vcell->getConservationLaw()->compute_left_eigenvectors(vcell->theMean,(mVector&) pp1,(double *)lev_x);
  mVector pp2(0,1,0);
  vcell->getConservationLaw()->compute_left_eigenvectors(vcell->theMean,pp2,(double*) lev_y);
   vcell->getConservationLaw()->compute_right_eigenvectors(vcell->theMean,(mVector&)pp1,(double*)rev_x);
  vcell->getConservationLaw()->compute_right_eigenvectors(vcell->theMean,pp2,(double*)rev_y);
  
  d=(vcell->fOrder+1)*(vcell->fOrder+1)-1; //the highest index in the coefficient
  vcell->limit = d; //for plotting - the highest coefficient that was not limited
  
  double slope[4],bottom_slope[4],top_slope[4],left_slope[4],right_slope[4];

  for (k=vcell->fOrder; k>0; k--) {// limiting (p,p), (p,p-1), (p-1,p) etc
    if (stopLimiter) break;
    
    for(q=0; q<cSize; q++) 
      {
	slope[q] =vcell->theFieldsCoefficients->getCopy(q,d);
	if (top) top_slope[q] = top->theFieldsCoefficients->get(q,d-2)-vcell->theFieldsCoefficients->get(q,d-2);
	if (bottom) bottom_slope[q] = vcell->theFieldsCoefficients->get(q,d-2)-bottom->theFieldsCoefficients->get(q,d-2);
	if (left) left_slope[q] = vcell->theFieldsCoefficients->get(q,d-1)-left->theFieldsCoefficients->get(q,d-1);
	if (right) right_slope[q] = right->theFieldsCoefficients->get(q,d-1)-vcell->theFieldsCoefficients->get(q,d-1);
      } 
    
    double ev;
    bool slope_limited=0;
    for(q=0; q<cSize; q++) {                 //limit coef (k,k) in vertical direction
      ev = dotprod(q,lev_y,slope);
      if (continue_limiting[q] &&(fabs(ev) > small_number))
	{
	  if (top && bottom)
	    tmp[q] = minmod(ev, dotprod(q,lev_y,top_slope)/sqrt(double (2*k+1)/(2*k-1)), 
			      dotprod(q,lev_y,bottom_slope)/sqrt(double (2*k+1)/(2*k-1)), wasLimited[q]);
	  else if (top && !bottom) 
	    tmp[q] = minmod(ev, dotprod(q,lev_y,top_slope)/sqrt(double (2*k+1)/(2*k-1)), wasLimited[q]);
	  else if (!top&&bottom)
	    tmp[q] = minmod(ev, dotprod(q,lev_y,bottom_slope)/sqrt(double (2*k+1)/(2*k-1)), wasLimited[q]);
	}
      else {
		  tmp[q] = ev;
		  wasLimited[q] = 0;
	  }
    }
    
    slope_limited = 0;
    for (q=0; q<cSize; q++) if (wasLimited[q]) slope_limited = 1;
    if (slope_limited) // if a component was limited, update the slope on (k,k)
      for(q=0; q<cSize; q++) 
	slope[q]=vcell->theFieldsCoefficients->getCopy(q,d) = dotprod(q,rev_y,tmp);
    
    for(q=0; q<cSize; q++) {           //limit coef (k,k) in vertical direction
      ev = dotprod(q,lev_x,slope);
      if (continue_limiting[q] &&(fabs(ev) > small_number))
	{
	  if (left && right)
	    tmp[q] = minmod(ev, dotprod(q,lev_x,right_slope)/sqrt(double (2*k+1)/(2*k-1)),
			    dotprod(q,lev_x,left_slope)/sqrt(double (2*k+1)/(2*k-1)),isLimited[q]);
	  else if (left && !right)
	    tmp[q] = minmod(ev,dotprod(q,lev_x,left_slope)/sqrt(double (2*k+1)/(2*k-1)),isLimited[q]);
	  else if (!left && right)
	    tmp[q] = minmod(ev,dotprod(q,lev_x,right_slope)/sqrt(double (2*k+1)/(2*k-1)),isLimited[q]);
	  
	  if (wasLimited[q] || isLimited[q]) {
	    if (q==0)
		  vcell->limit=d-1;
	  } 
	  else  continue_limiting[q] = 0;
	} 
      else  {
		  tmp[q] = ev; 
		  isLimited[q] = 0;
	  }
    }
    
    stopLimiter = 1;
    for(q=0; q<cSize; q++) if (continue_limiting[q]) stopLimiter = 0;  
    if (stopLimiter) break;
    
    slope_limited  =0;
    for (q=0; q<cSize; q++) if (isLimited[q]) slope_limited = 1;
    if (slope_limited) 
      for(q=0; q<cSize; q++) 
	vcell->theFieldsCoefficients->getCopy(q,d) = dotprod(q,rev_x,tmp);
    
    d--; // put the counter on next coef
  
	for (i=k-1; i>=0; i--) //limit (i,k)
      {
	if (i<k-1) i_km1 =(k-1)*(k-1)+2*i+1; else i_km1 = i*i+(k-1)*2; 
	if (i>0) im1_k = k*k+2*(i-1)+1; 
	for(q=0; q<cSize; q++) 
	  {
	    slope[q] =vcell->theFieldsCoefficients->getCopy(q,d);
	    if (top) top_slope[q] = top->theFieldsCoefficients->get(q,i_km1)-vcell->theFieldsCoefficients->get(q,i_km1);
	    if (bottom) bottom_slope[q] = vcell->theFieldsCoefficients->get(q,i_km1)-bottom->theFieldsCoefficients->get(q,i_km1);
	    if (left && i) left_slope[q] = vcell->theFieldsCoefficients->get(q,im1_k)-left->theFieldsCoefficients->get(q,im1_k);
	    if (right && i) right_slope[q] = right->theFieldsCoefficients->get(q,im1_k)-vcell->theFieldsCoefficients->get(q,im1_k);
	  } 
	if (i) {//limit (i,k) in x dir
	  for(q=0; q<cSize; q++) {
	    ev = dotprod(q,lev_x,slope);
	    if (continue_limiting[q] && (fabs(ev) > small_number)) 
	      {if (left && right)
		tmp[q] = minmod(ev,dotprod(q,lev_x,right_slope)/sqrt(double (2*i+1)/(2*i-1)),
				dotprod(q,lev_x,left_slope)/sqrt(double (2*i+1)/(2*i-1)),wasLimited[q]);
	      else if (left && !right)
		tmp[q] = minmod(ev,dotprod(q,lev_x,left_slope)/sqrt(double (2*i+1)/(2*i-1)),wasLimited[q]);
	      else if (!left && right)
		tmp[q] = minmod(ev,dotprod(q,lev_x,right_slope)/sqrt(double (2*i+1)/(2*i-1)), wasLimited[q]);
	      }
	    else tmp[q] = ev;
	  }
	  slope_limited = 0;
	  for (q=0; q<cSize; q++) if (wasLimited[q]) slope_limited = 1;
	    for(q=0; q<cSize; q++) 
	      vcell->theFieldsCoefficients->getCopy(q,d) = dotprod(q,rev_x,tmp);
	}
	else for(q=0; q<cSize; q++) wasLimited[q] = 0;
	
	for(q=0; q<cSize; q++) {//limit (i,k) in y dir
	  ev = dotprod(q,lev_y,slope);
	  if (continue_limiting[q] && (fabs(ev) > small_number)) {
	    if (top && bottom)
	      tmp[q] = minmod(ev, dotprod(q,lev_y,top_slope)/sqrt(double (2*k+1)/(2*k-1)),
			      dotprod(q,lev_y,bottom_slope)/sqrt(double (2*k+1)/(2*k-1)),isLimited[q]);
	    else if (top && !bottom) 
	      tmp[q] = minmod(ev,dotprod(q,lev_y,top_slope)/sqrt(double (2*k+1)/(2*k-1)), isLimited[q]);
	    else if (!top&&bottom)
	      tmp[q] = minmod(ev,dotprod(q,lev_y,bottom_slope)/sqrt(double (2*k+1)/(2*k-1)),isLimited[q]); 
	    
	    if (wasLimited[q] || isLimited[q]) 
	      {
		if (q==0) 
			vcell->limit=d-1;
			previous_coeff_limited[q]=1;
	      } 
	  else previous_coeff_limited[q] = 0;
	  }
	  else tmp[q] = ev;
	}
	
	slope_limited  =0;
	for (q=0; q<cSize; q++) if (isLimited[q]) slope_limited = 1;
	if (slope_limited) 
	  for(q=0; q<cSize; q++) 
	    vcell->theFieldsCoefficients->getCopy(q,d) = dotprod(q,rev_y,tmp);
	
	d--;
	
	k_im1 = k*k+2*(i-1);
	km1_i = (k-1)*(k-1)+2*i;
	
	for(q=0; q<cSize; q++)
	  {
	    slope[q] =vcell->theFieldsCoefficients->getCopy(q,d);
	    if (top && i) top_slope[q] = top->theFieldsCoefficients->get(q,k_im1)-vcell->theFieldsCoefficients->get(q,k_im1);
	    if (bottom&& i) bottom_slope[q] = vcell->theFieldsCoefficients->get(q,k_im1)-bottom->theFieldsCoefficients->get(q,k_im1);
	    if (left) left_slope[q] = vcell->theFieldsCoefficients->get(q,km1_i)-left->theFieldsCoefficients->get(q,km1_i);
	    if (right) right_slope[q] = right->theFieldsCoefficients->get(q,km1_i)-vcell->theFieldsCoefficients->get(q,km1_i);
	  } 
	
	if (i) { //limit (k,i) in y direction
	  for(q=0; q<cSize; q++) {
	    ev = dotprod(q,lev_y,slope);
	    if (continue_limiting[q] && (ev>small_number))
	      {
		if (top && bottom &&i)
		  tmp[q] = minmod(ev,dotprod(q,lev_y,top_slope)/sqrt(double (2*i+1)/(2*i-1)),
				  dotprod(q,lev_y,bottom_slope)/sqrt(double (2*i+1)/(2*i-1)),wasLimited[q]);
		else if (top && !bottom && i) 
		  tmp[q] = minmod(ev,dotprod(q,lev_y,top_slope)/sqrt(double (2*i+1)/(2*i-1)), wasLimited[q]);
		else if (!top&&bottom && i)
		tmp[q] = minmod(ev,dotprod(q,lev_y,bottom_slope)/sqrt(double (2*i+1)/(2*i-1)),wasLimited[q]);
	      }
	    else tmp[q] = ev;
	  }
	  slope_limited = 0;
	  for (q=0; q<cSize; q++) if (wasLimited[q]) slope_limited = 1;
	  if (slope_limited) 
	    for(q=0; q<cSize; q++) 
	      vcell->theFieldsCoefficients->getCopy(q,d) = dotprod(q,rev_y,tmp);
	  
	  if (slope_limited)
	    for(q=0; q<cSize; q++) 
	      slope[q] =vcell->theFieldsCoefficients->getCopy(q,d);
	}
	else for (q=0; q<cSize; q++) wasLimited[q] = 0;
	
	for(q=0; q<cSize; q++) { //limit (k,i) in x direction)
	  ev = dotprod(q,lev_x,slope);
	  if (continue_limiting[q] &&(fabs(ev) > small_number))
	    {
	      if (left && right)
		tmp[q] = minmod(ev,dotprod(q,lev_x,right_slope)/sqrt(double (2*k+1)/(2*k-1)),
				dotprod(q,lev_x,left_slope)/sqrt(double (2*k+1)/(2*k-1)),isLimited[q]);
	      else if (left && !right)
		tmp[q] = minmod(ev,dotprod(q,lev_x,left_slope)/sqrt(double (2*k+1)/(2*k-1)),isLimited[q]);
	      else if (!left && right)
		tmp[q] = minmod(ev,dotprod(q,lev_x,right_slope)/sqrt(double (2*k+1)/(2*k-1)), isLimited[q]);
	      
     	      if (wasLimited[q] || isLimited[q]) {
		if (q==0)
			vcell->limit=d-1;
	      } 
	      else {
		if (!previous_coeff_limited[q]) 
		  continue_limiting[q] = 0;
	      }
	    } 
	  else tmp[q] = ev;
	}
	
	stopLimiter = 1;
	for(q=0; q<cSize; q++) if (continue_limiting[q]) stopLimiter = 0;  
	if (stopLimiter) break;
	slope_limited  =0;
	for (q=0; q<cSize; q++) if (isLimited[q]) slope_limited = 1;
	if (slope_limited) 
	  for(q=0; q<cSize; q++) 
	    vcell->theFieldsCoefficients->getCopy(q,d) = dotprod(q,rev_x,tmp);
	
	d--;
      }
  }
}


#define SIGN(a)	  (a < 0.0 ? -1 : 1)
double minmod(double a,double b,double c, bool &isLimited)
{
if (SIGN(a)==SIGN(b) && SIGN(b)==SIGN(c))
{
double abs_a = fabs(a);
double abs_b = fabs(b);
double abs_c = fabs(c);
if (abs_a<abs_b && abs_a < abs_c) {isLimited=0;return a;}
else if (abs_b<abs_c) {isLimited=1;return b;}
else {isLimited=1; return c;}

}
else {isLimited = 1; return 0.0;}
}

double minmod(double a,double b,bool &isLimited)
{
if (SIGN(a)==SIGN(b))
{
	double tol = 1e-7;
if( abs(a - b) < tol*abs(a) ){isLimited = 0; return a;}
 if (fabs(a)<fabs(b)) {isLimited = 0; return a;}
 else  {
	isLimited = 1;return b;
	}
}
else {isLimited = 1;return 0.0;}
}


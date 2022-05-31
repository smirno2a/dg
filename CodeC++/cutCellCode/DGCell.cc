#include "DGCell.h"
#include "ConservationLaw.h"
#include "FunctionSpace.h"
#include "Integrator.h"
#include "Mapping.h"
#include "mEntity.h"
#include "Face.h"
#include "Edge.h"
#include "Vertex.h"
#include "FieldEvaluator.h"
#include "mMesh.h"
#include "Geometry.h"
#include <stdio.h>
#include <math.h>
#include <iomanip>

#define PI 3.14159265359
extern void invmat (double **, double **, int);

double ** allocateMatrix(int n)
{
  int i;
  double ** v = new double *[n];
  double *help = new double[n*n];
  for (i=0; i!= n; ++i)
    v[i] = &help[i*n];
  return v;
}



void freeMatrix( double **v)
{
  if(v)
    {
      if(v[0]) delete [] v[0];
      delete [] v;
    }
}


DGCell::DGCell (ConservationLaw*l, 
		mEntity*e, 
		FunctionSpace*f,
		FunctionSpace *er,
		Mapping *m,
		GaussIntegrator *I)
  : theConservationLaw(l),theMeshEntity(e),theFunctionSpace(f),
    theMapping(m),theGaussIntegrator(I)
{
  //theGeometry = g;
  int i;
  fSize                 = f->size();
  fOrder                = f->order();
  cSize                 = l->getNbFields();
  theFieldsCoefficients = new DG_Dofs(2,cSize,fSize);
  theMean               = new double[cSize];
  limSlope              = new double[cSize];
  theRightHandSide      = new double [fSize*cSize];
  deltaUp0              = new double[cSize];
  Up0                   = new double[cSize];
  Rp0                   = new double[cSize];
  
  for(i=0;i<cSize;++i)
    for(int j=0;j<fSize;++j) theFieldsCoefficients->get(i,j) = 0.;
  limit = 0;
  order = 2 * fOrder + theMapping->order()-1;
  
  theInvertMassMatrix = 0;
  init();
  
  pMin = pMax = ((Vertex*)theMeshEntity->get(0,0))->point();
  int Size0 = theMeshEntity->size(0);
  for(i=1;i<Size0;++i)
    {
      mPoint p1 = ((Vertex*)theMeshEntity->get(0,i))->point();
      if (p1(0) < pMin(0)) pMin(0) = p1(0);
      if (p1(1) < pMin(1)) pMin(1) = p1(1);
      if (p1(2) < pMin(2)) pMin(2) = p1(2);
      if (p1(0) > pMax(0)) pMax(0) = p1(0);
      if (p1(1) > pMax(1)) pMax(1) = p1(1);
      if (p1(2) > pMax(2)) pMax(2) = p1(2);
    }

  computeCellSize();
  computeVolume();
  computePerimeter();
  ZeroError();
  computeMaxH();

  for(int i = 0; i < 8 ; ++i){
    Neighbours[i] = NULL;
    neighbour_weights[i] = 0.0;
    Neighbours_cm[i] = NULL;
    neighbour_weights_cm[i] = 0.0;
  }
  for (size_t i = 8; i < 12; i++)
  {
    Neighbours[i] = NULL;
  }
  
  for(int i =0 ; i< 4; ++i){
    dist_intersection[i] = 1.0;
    dist_intersection_cm[i] = 1.0;
  }
}

DGCell::~DGCell ()
{
  delete theFieldsCoefficients;
  delete [] theMean;
  delete [] limSlope;
  delete [] theRightHandSide;
  delete theFunctionSpace;
  freeMatrix(theInvertMassMatrix);
}

void DGCell::cleanup()
{
  delete theMapping;
  theMapping = 0;
}

volumeGaussPoint::volumeGaussPoint (const volumeGaussPoint &other)
{
  int i;
  JacTimesWeight = other.JacTimesWeight;
  x = other.x;
  y = other.y;
  z = other.z;
  nbF = other.nbF;
  grads.reserve(nbF);
  grads.resize(nbF);
  fcts = other.fcts;
  for(i=0;i<nbF;++i)
    grads[i] = other.grads[i];
}

void DGCell::init ()
{
  int i,j,k;
  double** theMassMatrix;
  if (!theFunctionSpace->isOrthogonal())
    {
      freeMatrix(theInvertMassMatrix);
      theMassMatrix       = allocateMatrix(fSize);
      theInvertMassMatrix = allocateMatrix(fSize);
      for(j=0;j<fSize;++j)
	for(k=0;k<fSize;++k)
	  theMassMatrix[j][k] = 0.0;
    }
  
  nbPtGauss = theGaussIntegrator->nbIntegrationPoints(theMeshEntity,order);
  volumeGaussPoints.reserve(nbPtGauss);
  //Computing the integration points, weights, Jacobian.
  //x,y,z - physical coordinates
  double u,v,w,weight;
  for(i=0;i<nbPtGauss;++i)
    {
      volumeGaussPoint pg;
      pg.nbF = fSize;
      theGaussIntegrator->iPoint(theMeshEntity,i,order,u,v,w,weight); 
	 // u=0;v=1;w=0;
      theMapping->eval(u,v,w,pg.x,pg.y,pg.z);
      pg.grads.reserve(fSize);
      pg.grads.resize(fSize);
      detJac = theFunctionSpace->grads(u,v,w,theMapping,pg.grads);
      pg.JacTimesWeight = detJac*weight;
      for(j=0;j<fSize;++j) pg.grads[j] *= (pg.JacTimesWeight);
      //values of basis functions at (u,v,w)
      //theFunctionSpace->fcts(u,v,w,pg.fcts);
      pg.fcts=theGaussIntegrator->iFct(theMeshEntity,theFunctionSpace,i,order); 
      volumeGaussPoints.push_back(pg);
 
      ///Computing the Mass Matrix
      if (!theFunctionSpace->isOrthogonal())
	{
	  for( j=0;j<fSize;++j) {
	    for( k=0;k<fSize;++k){
	      theMassMatrix [j][k] += pg.fcts[k] * pg.fcts[j] * pg.JacTimesWeight;
	      //printf("%6.3f",theMassMatrix [j][k]);
            }
	  }
	}
    }
  if (!theFunctionSpace->isOrthogonal())
    {
      invmat(theMassMatrix,theInvertMassMatrix,fSize);
      freeMatrix(theMassMatrix);
    }
}

void DGCell::computeVolumeContribution (double t)
{
  // we compute int_{Cell} F grad w dCell
  int i,j,k;
  double rhs[MaxNbEqn];
  double *RHS;
  double val[MaxNbEqn];
  mVector fluxes[MaxNbEqn];
 
  for(i=0;i<nbPtGauss;++i)
    {
      volumeGaussPoint  &pg = pt(i);
      interpolate(pg.fcts, val);
      
      if(!theConservationLaw->isPhysical(val))
	{
	  //printf("non physical field found in VOLUME integral at point (%f,%f,%f) \n",pg.x, pg.y, pg.z ); 
	  //  printf("%f %f %f %f\n",val[0],val[1],val[2],val[3]);
	  double p;
	  int dim = theMeshEntity->getLevel();
	  if (dim==2)
	    p = (0.4)*(val[3]-0.5*(val[1]*val[1]+val[2]*val[2])/val[0]);
	  else 
	    p = (0.4)*(val[3]-0.5*(val[1]*val[1]+val[2]*val[2]+val[4]*val[4])/val[0]);
	  //printf("rho=%e p=%e \n",val[0],p);
	  
	 
	  for(j=0;j<cSize;++j)
	    for (k=1;k<fSize;++k)
	      theFieldsCoefficients->get(j,k)= 0.0;
	  interpolate(pg.fcts, val);
	  if (dim==2)
	    p = (0.4)*(val[3]-0.5*(val[1]*val[1]+val[2]*val[2])/val[0]);
	  else 
	    p = (0.4)*(val[3]-0.5*(val[1]*val[1]+val[2]*val[2]+val[4]*val[4])/val[0]);
	  //printf("Corrected values rho=%e p=%e \n",val[0],p);
	  if (p<0) exit(0);
	  ZeroRHS();
	  computeVolumeContribution(t);
	  break;
	}
      else
	{
	  const mPoint p(pg.x,pg.y,pg.z);  
	  theConservationLaw->RHS(val,p,t,rhs);
	  theConservationLaw->Fi(p,val,fluxes);
	  RHS = theRightHandSide;
	  const double *fcts = pg.fcts;
	  const double *rhs_const =rhs;
	  double tmp=0;
	  for(j=0;j<cSize;++j)
	    for(k=0;k<fSize;++k)		    
	      (*RHS++)+= fluxes[j] *pg.grads[k]  + rhs_const[j] * fcts[k]*pg.JacTimesWeight;
	}
    }
  //printf("\n");
  //for(j=0;j<9;++j) printf("%e\n",theRightHandSide[j]);
}

void DGCell::reverseVolumeContribution (double t)
{
  // we compute int_{Cell} F grad w dCell
  int i,j,k;
  double rhs[MaxNbEqn];
  double *RHS;
  double val[MaxNbEqn];
  mVector fluxes[MaxNbEqn];
 
  for(i=0;i<nbPtGauss;++i)
    {
      volumeGaussPoint  &pg = pt(i);
      interpolate(pg.fcts, val);
      
      const mPoint p(pg.x,pg.y,pg.z);  
      theConservationLaw->RHS(val,p,t,rhs);
      theConservationLaw->Fi(p,val,fluxes);
      RHS = theRightHandSide;
      const double *fcts = pg.fcts;
      const double *rhs_const =rhs;
      double tmp=0;
      for(j=0;j<cSize;++j)
	for(k=0;k<fSize;++k)		    
	  (*RHS++)-= fluxes[j] *pg.grads[k]  + rhs_const[j] * fcts[k]*pg.JacTimesWeight;
    }
  //printf("\n");
  //for(j=0;j<9;++j) printf("%e\n",theRightHandSide[j]);
}

void DGCell::dinterpolate(double u, double v, double w, mVector *grads)
{

  //  if(theMeshEntity->isAdjacencyCreated(2))printf("oups ...\n");

  vector<mVector>  gr;
  gr.reserve(fSize);  
  theFunctionSpace->grads(u,v,w,theMapping,gr);

  mVector g = gr[0];
  for(int i=0;i<cSize;++i)
    {
      // grads[i] = g * (theFieldsCoefficients[i])[0];
      grads[i] = g * theFieldsCoefficients->get(i,0);
    }
  
  for(int j=1;fSize;++j)
    {
      g = gr[j];
      for(int i=0;i<cSize;++i)
	{
	  grads[i] += (g * theFieldsCoefficients->get(i,j));
	}
    }
}
    
double DGCell::getError() const
{
	return sqrt(error);
}


void DGCell::computeMaxH()
{
double h = 0.;
mPoint p0 =((Vertex*)theMeshEntity->get(0,0))->point();
  
  for(int i=1;i<theMeshEntity->size(0);++i)
    {
      mPoint p1 = ((Vertex*)theMeshEntity->get(0,i))->point();
      mVector v(p0,p1);
	  if (h * h < v*v) h = sqrt(v*v);
	  p0 = p1;
	 } 	  
  maxH = h;
}

void DGCell::computeMinH()
{
double h = 0.;
mPoint p0 =((Vertex*)theMeshEntity->get(0,0))->point();
  
  for(int i=1;i<theMeshEntity->size(0);++i)
    {
      mPoint p1 = ((Vertex*)theMeshEntity->get(0,i))->point();
      mVector v(p0,p1);
	  if (h > v*v) h = v*v;
	  p0 = p1;
	 } 
  minH = h;
}

double DGCell::getMaxH() const
{
	return maxH;
}

double DGCell::getMinH() const
{
	return minH;
}

void DGCell::find_center_of_mass(mPoint *cm, mPoint *vertex)
{
  mPoint c1, c2, c3, c4;
  double line1[3], line2[3], det;
  // cell center of triangle 124
  c1(0) = (vertex[0](0) + vertex[1](0) + vertex[3](0))/3.;
  c1(1) = (vertex[0](1) + vertex[1](1) + vertex[3](1))/3.;

  // cell center of triangle 123
  c2(0) = (vertex[0](0) + vertex[1](0) + vertex[2](0))/3.;
  c2(1) = (vertex[0](1) + vertex[1](1) + vertex[2](1))/3.;

  // cell center of triangle 234
  c3(0) = (vertex[1](0) + vertex[2](0) + vertex[3](0))/3.;
  c3(1) = (vertex[1](1) + vertex[2](1) + vertex[3](1))/3.;

  // cell center of triangle 341
  c4(0) = (vertex[0](0) + vertex[2](0) + vertex[3](0))/3.;
  c4(1) = (vertex[0](1) + vertex[2](1) + vertex[3](1))/3.;
  
  // line 1 a1*x+b1*y + c1 = 0 // c1-c3
  line1[0] = c3(1)-c1(1); line1[1] = -(c3(0) - c1(0)); line1[2] = -c1(1)*line1[1] - c1(0)*line1[0];

  // line 2 a2*x+b2*y + c2 = 0 // c2-c4
  line2[0] = c4(1)-c2(1); line2[1] = -(c4(0) - c2(0)); line2[2] = -c2(1)*line2[1] - c2(0)*line2[0];

  det = line1[0]*line2[1] - line1[1]*line2[0];
  mPoint point_interesct((line1[1]*line2[2] - line1[2]*line2[1])/det, (line1[2]*line2[0] - line1[0]*line2[2])/det);

  (*cm)(0) = point_interesct(0);
  (*cm)(1) = point_interesct(1);        

  return;
}

int DGCell::createReconstructionNeighbourhood()
{
   //printf("createReconstrcutionNeighbourhood is being called\n");
  list<mEntity*>::const_iterator it;
  vector<mPoint>::const_iterator mit;
  vector<DGCell *>::const_iterator cit;
    
  list<mEntity *> allVertneigh;
  vector<DGCell *> Vertneigh;
  vector<mPoint> neighCenter;
  vector<mPoint> neighCentermass;
  mEntity *ent;

  int no_vertex = theMeshEntity->size(0);
	int no_edges = theMeshEntity->size(1);
  int neighbourhood_size = 0;
  int flag = 0;
  mPoint v[4], v_neigh[4], centroid(0.,0.0), temp_c(0.0,0.0), temp_n(0.0,0.0);
  mPoint centermass(0.,0.), Cellcentermass(0.,0.);
  double angle[20];
  int index[20];
  int temp_reconstruct_index[8] = {-1,-1,-1,-1,-1,-1,-1,-1};
  for(int i = 0; i < 20; ++i){
    index[i] = i;
    angle[i] = 0;
  } 

  for(int i = 0; i < 8 ; ++i){
    Neighbours[i] = NULL;
    neighbour_weights[i] = 0.0;
  }
  for(int i =0 ; i< 4; ++i){
    dist_intersection[i] = 0.0;
  }

  for(int k=0;k<no_vertex;k++)
      {
        ent = theMeshEntity->get(0,k); 
        v[k] = ((Vertex*)ent)->point();
        centroid += v[k];
        for(int i=0;i<ent->size(2);i++)
          {
            allVertneigh.push_back(ent->get(2,i));
            //ent->get(n,i)->print();
          }
      }
      centroid = centroid * 0.25;

      cellCentroid(0) = centroid(0); cellCentroid(1) = centroid(1);

      find_center_of_mass(&Cellcentermass, v);
  
  // Setting up the coefficients of determinant of Jacobian
  //double coeffx[4], coeff_neighx[4], coeffy[4], coeff_neighy[4];
  double coeff_neighx[4], coeff_neighy[4];
	double Jaccoeff[3], Jaccoeff_neigh[3];
	coeffx[0] = 0.25*(v[1](0) + v[2](0) + (v[0](0) + v[3](0)));
	coeffx[1] = 0.25*(v[1](0) + v[2](0) - (v[0](0) + v[3](0)));
	coeffx[2] = 0.25*(v[3](0) + v[2](0) - (v[0](0) + v[1](0)));
	coeffx[3] = 0.25*(v[0](0) + v[2](0) - (v[3](0) + v[1](0)));
	coeffy[0] = 0.25*(v[1](1) + v[2](1) + (v[0](1) + v[3](1)));
	coeffy[1] = 0.25*(v[1](1) + v[2](1) - (v[0](1) + v[3](1)));
	coeffy[2] = 0.25*(v[3](1) + v[2](1) - (v[0](1) + v[1](1)));
	coeffy[3] = 0.25*(v[0](1) + v[2](1) - (v[3](1) + v[1](1)));
	Jaccoeff[0] = coeffx[1]*coeffy[2] - coeffx[2]*coeffy[1];
	Jaccoeff[1] = coeffx[1]*coeffy[3] - coeffx[3]*coeffy[1];
	Jaccoeff[2] = coeffx[3]*coeffy[2] - coeffx[2]*coeffy[3];
  scaling_derivative[0] = sqrt(coeffx[1]*coeffx[1] + coeffy[1]*coeffy[1] );
  scaling_derivative[1] = sqrt(coeffx[2]*coeffx[2] + coeffy[2]*coeffy[2] );

  dist_intersection[0] = dist_intersection[1] = scaling_derivative[0];
  dist_intersection[2] = dist_intersection[3] = scaling_derivative[1];

  alpha[0] = Jaccoeff[1]/Jaccoeff[0]; alpha[1] = Jaccoeff[2]/Jaccoeff[0]; 
  
  for(it = allVertneigh.begin();it!=allVertneigh.end();++it)
      {
		DGCell *cell = (DGCell*)(*it)->getCell();
    flag = 1;
    for(int j = 0; j < neighbourhood_size; ++j){
      if (Vertneigh[j] == cell) {
        flag = 0; break;
      }
    }
    if(cell != this && flag){

      Vertneigh.push_back(cell);
      temp_c = temp_c*0.0;

      // Setting up detector for Quads
      for(int k=0;k<theMeshEntity->size(0);k++)
      {
        v_neigh[k] = ((Vertex*)cell->theMeshEntity->get(0,k))->point();
        temp_c += v_neigh[k];
      }
      neighCenter.push_back(temp_c*0.25);
      find_center_of_mass(&centermass, v_neigh);
      neighCentermass.push_back(centermass);

      temp_n = neighCenter[neighbourhood_size] - centroid;
      temp_n = temp_n * (1.0/sqrt(temp_n(0)*temp_n(0) + temp_n(1)*temp_n(1)));
      if (temp_n(0) > 0 && temp_n(1) > 0) {
        angle[neighbourhood_size] = atan(temp_n(1)/temp_n(0));
      }
      else if (temp_n(0) < 0 && temp_n(1) > 0) {
        angle[neighbourhood_size] = atan(temp_n(1)/temp_n(0)) + PI;
      }
      else if (temp_n(0) < 0 && temp_n(1) < 0) {
        angle[neighbourhood_size] = atan(temp_n(1)/temp_n(0)) + PI;
      }
      else if (temp_n(0) > 0 && temp_n(1) < 0) {
        angle[neighbourhood_size] = atan(temp_n(1)/temp_n(0)) + 2.*PI;
      }
      else {
        if (-1e-8 < temp_n(1) && temp_n(1) < 1e-8 && temp_n(0) > 0){
          angle[neighbourhood_size] = 0;
        }
        else if (-1e-8 < temp_n(0) && temp_n(0) < 1e-8 && temp_n(1) > 0){
          angle[neighbourhood_size] = 0.5*PI;
        }
        else if (-1e-8 < temp_n(1) && temp_n(1) < 1e-8 && temp_n(0) < 0){
          angle[neighbourhood_size] = PI;
        }
        else{
          angle[neighbourhood_size] = 1.5*PI;
        }
      }
      neighbourhood_size++;
    }
      }

      // Sorting the angles
      double min_angle = angle[0];
      int temp_index = 0;
      for(int i = 0 ; i < neighbourhood_size-1; ++i){
        for(int j= i+1 ; j < neighbourhood_size; ++j){
          if (angle[j] < angle[i]){
            min_angle = angle[j];
            angle[j] = angle[i];
            angle[i] = min_angle;

            temp_index = index[j]; index[j] = index[i]; index[i] = temp_index;
          }
        }
      }
      
       double Jac_sum = 0.0;
      int Jac_flag = 0;
      //if(-0.2 < centroid(0) && centroid(0) < 0.2 &&  -0.2 < centroid(1) && centroid(1) < 0.2)
      {
        Jac_sum = fabs(Jaccoeff[1]/Jaccoeff[0]) + fabs(Jaccoeff[2]/Jaccoeff[0]);
        if (Jac_sum > 0.4){
          Jac_flag = 4;
        //printf("Centroid : (%.12e,%.12e)  ,",centroid(0),centroid(1));
        //printf("Jacobian sum : %.12e\n",Jac_sum);
        }
        else if(Jac_sum > 0.2){
          Jac_flag = 2;
        }
        else if(Jac_sum > 0.1){
          Jac_flag = 1;
        }
        else {
          Jac_flag = 0;
        }
        
      }

      int neigh_num_edges = 0; 
      /*for (cit = Vertneigh.begin(); cit != Vertneigh.end(); cit++)
      {
        // To make sure element is not boundary adjacent
        DGCell *cell = (DGCell*)(*cit);
        neigh_num_edges = cell->theMeshEntity->size(1);
         for(int k =0 ; k < neigh_num_edges; ++k){
          ent = cell->theMeshEntity->get(1,k);
          if(ent->getClassification()->getId() == 2 ||ent->getClassification()->getId() == 30000 || ent->getClassification()->getId() == 610){
            return Jac_flag;
          }
        }

      }*/
      

      double  line1[3] = {0,0,0}, line2[3] = {0,0,0}, dx = 0.0, dy = 0.0, det = 0.;
      double dist1 = 0., dist2 = 0, dist3 = 0., dist = 0.0;
      double  linev1[3] = {0,0,0}, linev2[3] = {0,0,0}, dxv = 0.0, dyv = 0.0, detv = 0.;
      double temp_linev1[3] = {0,0,0}, temp_linev2[3] = {0,0,0}, detv1 = 0;
      double distv1 = 0., distv2 = 0, distv3 = 0., distv = 0.0;
      mPoint edgeMidpoints[4], V[2];
      // Defining edge midpoints
      for(int i = 0; i < 4; ++i){
        edgeMidpoints[i] = v[i] * 0.5 + v[(i+1)%4] * 0.5;
      }
      V[0](0) = coeffx[1]; V[0](1) = coeffy[1]; //;edgeMidpoints[1] - edgeMidpoints[3];
      V[1](0) = coeffx[2]; V[1](1) = coeffy[2]; //;edgeMidpoints[2] - edgeMidpoints[0];
      V[0] = V[0] * (1.0/sqrt(V[0](0)*V[0](0) + V[0](1)*V[0](1) ));
      V[1] = V[1] * (1.0/sqrt(V[1](0)*V[1](0) + V[1](1)*V[1](1) ));

      // Finding interpolation points
      for(int i = 0; i < neighbourhood_size; ++i){
        // line 1 a1*x+b1*y + c1 = 0
        dx = neighCenter[index[(i+1)%neighbourhood_size]](0) - neighCenter[index[i]](0);
        dy = neighCenter[index[(i+1)%neighbourhood_size]](1) - neighCenter[index[i]](1);
        line1[0] = -dy; line1[1] = dx; line1[2] = neighCenter[index[i]](0)*dy - neighCenter[index[i]](1)*dx;

        // line 2 a2*x+b2*y + c2 = 0
        line2[0] = -V[0](1); line2[1] = V[0](0); line2[2] = centroid(0) * V[0](1) - centroid(1) * V[0](0);

        det = line1[0]*line2[1] - line1[1]*line2[0];

        // centermass line 1 a1*x+b1*y + c1 = 0
        dxv = neighCentermass[index[(i+1)%neighbourhood_size]](0) - neighCentermass[index[i]](0);
        dyv = neighCentermass[index[(i+1)%neighbourhood_size]](1) - neighCentermass[index[i]](1);
        linev1[0] = -dyv; linev1[1] = dxv; linev1[2] = neighCentermass[index[i]](0)*dyv - neighCentermass[index[i]](1)*dxv;

        // centermass line 2 a2*x+b2*y + c2 = 0
        linev2[0] = -V[0](1); linev2[1] = V[0](0); linev2[2] = Cellcentermass(0) * V[0](1) - Cellcentermass(1) * V[0](0);

        detv = linev1[0]*linev2[1] - linev1[1]*linev2[0];
     
        if (sqrt(det*det) > 1e-10){
           mPoint point_interesct((line1[1]*line2[2] - line1[2]*line2[1])/det, (line1[2]*line2[0] - line1[0]*line2[2])/det);
           dist = centroid.pointdistance(point_interesct);
           dist1 = point_interesct.pointdistance(neighCenter[index[i]]);
           dist2 = point_interesct.pointdistance(neighCenter[index[(i+1)%neighbourhood_size]]);
           dist3 = neighCenter[index[i]].pointdistance(neighCenter[index[(i+1)%neighbourhood_size]]);
           if(dist1 + dist2 < dist3+1e-10){
             if ( (point_interesct(0) - centroid(0))*V[0](0) + (point_interesct(1) - centroid(1))*V[0](1) > 0  && dist > dist_intersection[0] ) {

              mPoint point_interesctv((linev1[1]*linev2[2] - linev1[2]*linev2[1])/detv, (linev1[2]*linev2[0] - linev1[0]*linev2[2])/detv);
              distv = Cellcentermass.pointdistance(point_interesctv);
              distv1 = point_interesctv.pointdistance(neighCentermass[index[i]]);
              distv2 = point_interesctv.pointdistance(neighCentermass[index[(i+1)%neighbourhood_size]]);
              distv3 = neighCentermass[index[i]].pointdistance(neighCentermass[index[(i+1)%neighbourhood_size]]);
              if(0.043 < centroid(0) && centroid(0) < 0.047 &&  -0.095 < centroid(1) && centroid(1) < -0.091621){
              //if(0.64 < centroid(0) && centroid(0) < 0.651 &&  0.613 < centroid(1) && centroid(1) < 0.621){
              //  printf("Center mass intersection distances : %.12e, %.12e, %.12e\n",distv1, distv2, distv3);
              }
              
              if(distv1 + distv2 < distv3+1e-10)
              {
                  Neighbours[0] = Vertneigh[index[i]]; Neighbours[1] = Vertneigh[index[(i+1)%neighbourhood_size]]; 
                  Neighbours[8] = Vertneigh[index[(i+1)%neighbourhood_size]]; 
                  temp_reconstruct_index[0] = index[i]; temp_reconstruct_index[1] = index[(i+1)%neighbourhood_size];
                  intersection_compvar[0] = (point_interesct(0) + point_interesct(1) - coeffx[0] - coeffy[0])/(coeffx[1] + coeffy[1]);   

                  dist_intersection[0] = dist;
                  neighbour_weights[0] = dist2 / (dist1 + dist2);
                  neighbour_weights[1] = dist1 / (dist1 + dist2);
                  
                  dist_intersection_cm[0] = distv;
                  neighbour_weights_cm[0] = distv2 / (distv1 + distv2);
                  neighbour_weights_cm[1] = distv1 / (distv1 + distv2);    
         
              }
              else {
                int temp_i = i, temp_j = (i+1)%neighbourhood_size;
                 if (distv1 > distv2){

                  Neighbours[1] = Vertneigh[index[i]]; Neighbours[0] = Vertneigh[index[(i+1)%neighbourhood_size]]; 
                  Neighbours[8] = Vertneigh[index[(i+2)%neighbourhood_size]]; 
                  temp_reconstruct_index[1] = index[i]; temp_reconstruct_index[0] = index[(i+1)%neighbourhood_size];
                  intersection_compvar[0] = (point_interesct(0) + point_interesct(1) - coeffx[0] - coeffy[0])/(coeffx[1] + coeffy[1]);   

                  dist_intersection[0] = dist;
                  neighbour_weights[0] = dist1 / (dist1 + dist2);
                  neighbour_weights[1] = dist2 / (dist1 + dist2);

                  temp_i = temp_j;
                  temp_j = (i+2)%neighbourhood_size;   
                 }
                 else{

                  Neighbours[0] = Vertneigh[index[i]]; Neighbours[1] = Vertneigh[index[(i+1)%neighbourhood_size]]; 
                  Neighbours[8] = Vertneigh[index[(i>0)?(i-1)%neighbourhood_size:neighbourhood_size-1]]; 
                  temp_reconstruct_index[0] = index[i]; temp_reconstruct_index[1] = index[(i+1)%neighbourhood_size];
                  intersection_compvar[0] = (point_interesct(0) + point_interesct(1) - coeffx[0] - coeffy[0])/(coeffx[1] + coeffy[1]);   

                  dist_intersection[0] = dist;
                  neighbour_weights[0] = dist2 / (dist1 + dist2);
                  neighbour_weights[1] = dist1 / (dist1 + dist2);

                  temp_i = i;
                  temp_j = (i>0)?i-1:neighbourhood_size-1;
                 }

                 // centermass line 1 a1*x+b1*y + c1 = 0
                dxv = neighCentermass[index[temp_j]](0) - neighCentermass[index[temp_i]](0);
                dyv = neighCentermass[index[temp_j]](1) - neighCentermass[index[temp_i]](1);
                temp_linev1[0] = -dyv; temp_linev1[1] = dxv; temp_linev1[2] = neighCentermass[index[temp_i]](0)*dyv - neighCentermass[index[temp_i]](1)*dxv;

                // centermass line 2 a2*x+b2*y + c2 = 0
                temp_linev2[0] = -V[0](1); temp_linev2[1] = V[0](0); temp_linev2[2] = Cellcentermass(0) * V[0](1) - Cellcentermass(1) * V[0](0);

                detv1 = temp_linev1[0]*temp_linev2[1] - temp_linev1[1]*temp_linev2[0];

                 mPoint point_interesctv1((temp_linev1[1]*temp_linev2[2] - temp_linev1[2]*temp_linev2[1])/detv1, (temp_linev1[2]*temp_linev2[0] - temp_linev1[0]*temp_linev2[2])/detv1);
                distv = Cellcentermass.pointdistance(point_interesctv1);
                distv1 = point_interesctv1.pointdistance(neighCentermass[index[temp_i]]);
                distv2 = point_interesctv1.pointdistance(neighCentermass[index[temp_j]]);
                distv3 = neighCentermass[index[temp_i]].pointdistance(neighCentermass[index[temp_j]]);
                
                dist_intersection_cm[0] = distv;
                neighbour_weights_cm[0] = distv2 / (distv1 + distv2);
                neighbour_weights_cm[1] = distv1 / (distv1 + distv2);


              }
            
              
               for(int j = 0; j < 2; ++j){
                  for(int k=0;k<theMeshEntity->size(0);k++)
                  {
                    v_neigh[k] = ((Vertex*)Neighbours[j]->theMeshEntity->get(0,k))->point();
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
                  Jaccoeff_neigh[2] = coeff_neighx[3]*coeff_neighy[2] - coeff_neighx[2]*coeff_neighy[3];
                  coeff_alpha[j] = (1./Jaccoeff_neigh[0]) * (coeff_neighy[2]* (coeffx[2] + coeffx[3]*intersection_compvar[0]) - coeff_neighx[2] * (coeffy[2] + coeffy[3]*intersection_compvar[0]) );
                  coeff_beta[j] = (1./Jaccoeff_neigh[0]) * (-coeff_neighy[1]* (coeffx[2] + coeffx[3]*intersection_compvar[0]) + coeff_neighx[1] * (coeffy[2] + coeffy[3]*intersection_compvar[0]) );
               }
                }
             else if( (point_interesct(0) - centroid(0))*V[0](0) + (point_interesct(1) - centroid(1))*V[0](1) < 0  && dist > dist_intersection[1]) {

                mPoint point_interesctv((linev1[1]*linev2[2] - linev1[2]*linev2[1])/detv, (linev1[2]*linev2[0] - linev1[0]*linev2[2])/detv);
              distv = Cellcentermass.pointdistance(point_interesctv);
              distv1 = point_interesctv.pointdistance(neighCentermass[index[i]]);
              distv2 = point_interesctv.pointdistance(neighCentermass[index[(i+1)%neighbourhood_size]]);
              distv3 = neighCentermass[index[i]].pointdistance(neighCentermass[index[(i+1)%neighbourhood_size]]);
              //if(0.043 < centroid(0) && centroid(0) < 0.047 &&  -0.095 < centroid(1) && centroid(1) < -0.091621){
              //if(0.64 < centroid(0) && centroid(0) < 0.651 &&  0.613 < centroid(1) && centroid(1) < 0.621){
               // printf("Center mass intersection distances : %.12e, %.12e, %.12e\n",distv1, distv2, distv3);
              //}
              
              if(distv1 + distv2 < distv3+1e-10)
              {  
                Neighbours[2] = Vertneigh[index[i]]; Neighbours[3] = Vertneigh[index[(i+1)%neighbourhood_size]]; 
                Neighbours[9] = Vertneigh[index[(i+1)%neighbourhood_size]];
                temp_reconstruct_index[2] = index[i]; temp_reconstruct_index[3] = index[(i+1)%neighbourhood_size];
                intersection_compvar[1] = (point_interesct(0) + point_interesct(1) - coeffx[0] - coeffy[0])/(coeffx[1] + coeffy[1]);

                dist_intersection[1] = dist;
                neighbour_weights[2] = dist2 / (dist1 + dist2);
                neighbour_weights[3] = dist1 / (dist1 + dist2);
               

                dist_intersection_cm[1] = distv;
                neighbour_weights_cm[2] = distv2 / (distv1 + distv2);
                neighbour_weights_cm[3] = distv1 / (distv1 + distv2);
               
                }
              else{
                int temp_i = i, temp_j = (i+1)%neighbourhood_size;
                 if (distv1 > distv2){

                  Neighbours[3] = Vertneigh[index[i]]; Neighbours[2] = Vertneigh[index[(i+1)%neighbourhood_size]]; 
                  Neighbours[9] = Vertneigh[index[(i+2)%neighbourhood_size]]; 
                  temp_reconstruct_index[3] = index[i]; temp_reconstruct_index[2] = index[(i+1)%neighbourhood_size];
                  intersection_compvar[1] = (point_interesct(0) + point_interesct(1) - coeffx[0] - coeffy[0])/(coeffx[1] + coeffy[1]);   

                  dist_intersection[1] = dist;
                  neighbour_weights[2] = dist1 / (dist1 + dist2);
                  neighbour_weights[3] = dist2 / (dist1 + dist2);

                  temp_i = temp_j;
                  temp_j = (i+2)%neighbourhood_size;   
                 }
                 else{

                  Neighbours[2] = Vertneigh[index[i]]; Neighbours[3] = Vertneigh[index[(i+1)%neighbourhood_size]]; 
                  Neighbours[9] = Vertneigh[index[(i>0)?(i-1)%neighbourhood_size:neighbourhood_size-1]]; 
                  temp_reconstruct_index[2] = index[i]; temp_reconstruct_index[3] = index[(i+1)%neighbourhood_size];
                  intersection_compvar[1] = (point_interesct(0) + point_interesct(1) - coeffx[0] - coeffy[0])/(coeffx[1] + coeffy[1]);   

                  dist_intersection[1] = dist;
                  neighbour_weights[2] = dist2 / (dist1 + dist2);
                  neighbour_weights[3] = dist1 / (dist1 + dist2);

                  temp_i = i;
                  temp_j = (i>0)?i-1:neighbourhood_size-1;
                 }

                 // centermass line 1 a1*x+b1*y + c1 = 0
                dxv = neighCentermass[index[temp_j]](0) - neighCentermass[index[temp_i]](0);
                dyv = neighCentermass[index[temp_j]](1) - neighCentermass[index[temp_i]](1);
                temp_linev1[0] = -dyv; temp_linev1[1] = dxv; temp_linev1[2] = neighCentermass[index[temp_i]](0)*dyv - neighCentermass[index[temp_i]](1)*dxv;

                // centermass line 2 a2*x+b2*y + c2 = 0
                temp_linev2[0] = -V[0](1); temp_linev2[1] = V[0](0); temp_linev2[2] = Cellcentermass(0) * V[0](1) - Cellcentermass(1) * V[0](0);

                detv1 = temp_linev1[0]*temp_linev2[1] - temp_linev1[1]*temp_linev2[0];

                 mPoint point_interesctv1((temp_linev1[1]*temp_linev2[2] - temp_linev1[2]*temp_linev2[1])/detv1, (temp_linev1[2]*temp_linev2[0] - temp_linev1[0]*temp_linev2[2])/detv1);
                distv = Cellcentermass.pointdistance(point_interesctv1);
                distv1 = point_interesctv1.pointdistance(neighCentermass[index[temp_i]]);
                distv2 = point_interesctv1.pointdistance(neighCentermass[index[temp_j]]);
                distv3 = neighCentermass[index[temp_i]].pointdistance(neighCentermass[index[temp_j]]);
                
                dist_intersection_cm[1] = distv;
                neighbour_weights_cm[2] = distv2 / (distv1 + distv2);
                neighbour_weights_cm[3] = distv1 / (distv1 + distv2);


              }
                
               
              
               for(int j = 2; j < 4; ++j){
                  for(int k=0;k<theMeshEntity->size(0);k++)
                  {
                    v_neigh[k] = ((Vertex*)Neighbours[j]->theMeshEntity->get(0,k))->point();
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
                  Jaccoeff_neigh[2] = coeff_neighx[3]*coeff_neighy[2] - coeff_neighx[2]*coeff_neighy[3];
                  coeff_alpha[j] = (1./Jaccoeff_neigh[0]) * (coeff_neighy[2]* (coeffx[2] + coeffx[3]*intersection_compvar[1]) - coeff_neighx[2] * (coeffy[2] + coeffy[3]*intersection_compvar[1]) );
                  coeff_beta[j] = (1./Jaccoeff_neigh[0]) * (-coeff_neighy[1]* (coeffx[2] + coeffx[3]*intersection_compvar[1]) + coeff_neighx[1] * (coeffy[2] + coeffy[3]*intersection_compvar[1]) );
               }
             }
           }
        }

        // line 2 a2*x+b2*y + c2 = 0
        line2[0] = -V[1](1); line2[1] = V[1](0); line2[2] = centroid(0) * V[1](1) - centroid(1) * V[1](0);

        det = line1[0]*line2[1] - line1[1]*line2[0];

        // centermass line 2 a2*x+b2*y + c2 = 0
        linev2[0] = -V[1](1); linev2[1] = V[1](0); linev2[2] = Cellcentermass(0) * V[1](1) - Cellcentermass(1) * V[1](0);

        detv = linev1[0]*linev2[1] - linev1[1]*linev2[0];



        if (sqrt(det*det) >  1e-10){
           mPoint point_interesct((line1[1]*line2[2] - line1[2]*line2[1])/det, (line1[2]*line2[0] - line1[0]*line2[2])/det);
           dist = centroid.pointdistance(point_interesct);
           dist1 = point_interesct.pointdistance(neighCenter[index[i]]);
           dist2 = point_interesct.pointdistance(neighCenter[index[(i+1)%neighbourhood_size]]);
           dist3 = neighCenter[index[i]].pointdistance(neighCenter[index[(i+1)%neighbourhood_size]]);
           if(dist1 + dist2 < dist3+1e-10){
             if ( (point_interesct(0) - centroid(0))*V[1](0) + (point_interesct(1) - centroid(1))*V[1](1) > 0 && dist > dist_intersection[2] ) {

              mPoint point_interesctv((linev1[1]*linev2[2] - linev1[2]*linev2[1])/detv, (linev1[2]*linev2[0] - linev1[0]*linev2[2])/detv);
              distv = Cellcentermass.pointdistance(point_interesctv);
              distv1 = point_interesctv.pointdistance(neighCentermass[index[i]]);
              distv2 = point_interesctv.pointdistance(neighCentermass[index[(i+1)%neighbourhood_size]]);
              distv3 = neighCentermass[index[i]].pointdistance(neighCentermass[index[(i+1)%neighbourhood_size]]);
              //if(0.64 < centroid(0) && centroid(0) < 0.651 &&  0.613 < centroid(1) && centroid(1) < 0.621){
              //  printf("Center mass intersection distances : %.12e, %.12e, %.12e\n",distv1, distv2, distv3);
              //  printf("Center mass : (%.12e,%.12e), Neighbouring  center mass : (%.12e,%.12e), (%.12e,%.12e)\n ",Cellcentermass(0), Cellcentermass(1), neighCentermass[index[i]](0),neighCentermass[index[i]](1), neighCentermass[index[(i+1)%neighbourhood_size]](0), neighCentermass[index[(i+1)%neighbourhood_size]](1));
              //  printf("Centeroid : (%.12e,%.12e), Neighbouring  centeroid : (%.12e,%.12e), (%.12e,%.12e)\n ",centroid(0), centroid(1), neighCenter[index[i]](0),neighCenter[index[i]](1), neighCenter[index[(i+1)%neighbourhood_size]](0), neighCenter[index[(i+1)%neighbourhood_size]](1));
              //}
              //if(0.043 < centroid(0) && centroid(0) < 0.047 &&  -0.095 < centroid(1) && centroid(1) < -0.091621){
              //if(0.64 < centroid(0) && centroid(0) < 0.651 &&  0.613 < centroid(1) && centroid(1) < 0.621){
              //  printf("Center mass intersection distances : %.12e, %.12e, %.12e\n",distv1, distv2, distv3);
              //}
                
              if(distv1 + distv2 < distv3+1e-10)
              {
                Neighbours[4] = Vertneigh[index[i]]; Neighbours[5] = Vertneigh[index[(i+1)%neighbourhood_size]];
                Neighbours[10] = Vertneigh[index[(i+1)%neighbourhood_size]];
                temp_reconstruct_index[4] = index[i]; temp_reconstruct_index[5] = index[(i+1)%neighbourhood_size];
                intersection_compvar[2] = (point_interesct(0) + point_interesct(1) - coeffx[0] - coeffy[0])/(coeffx[2] + coeffy[2]);
                
                dist_intersection[2] = dist;
                neighbour_weights[4] = dist2 / (dist1 + dist2);
                neighbour_weights[5] = dist1 / (dist1 + dist2);
               
                 dist_intersection_cm[2] = distv;
                 neighbour_weights_cm[4] = distv2 / (distv1 + distv2);
                 neighbour_weights_cm[5] = distv1 / (distv1 + distv2);    
          
               
               }
              else{
                int temp_i = i, temp_j = (i+1)%neighbourhood_size;
                 if (distv1 >= distv2){

                  Neighbours[5] = Vertneigh[index[i]]; Neighbours[4] = Vertneigh[index[(i+1)%neighbourhood_size]]; 
                  Neighbours[10] = Vertneigh[index[(i+2)%neighbourhood_size]]; 
                  temp_reconstruct_index[5] = index[i]; temp_reconstruct_index[4] = index[(i+1)%neighbourhood_size];
                  intersection_compvar[2] = (point_interesct(0) + point_interesct(1) - coeffx[0] - coeffy[0])/(coeffx[2] + coeffy[2]);   

                  dist_intersection[2] = dist;
                  neighbour_weights[4] = dist1 / (dist1 + dist2);
                  neighbour_weights[5] = dist2 / (dist1 + dist2);

                  temp_i = temp_j;
                  temp_j = (i+2)%neighbourhood_size;   
                 }
                 else{

                  Neighbours[4] = Vertneigh[index[i]]; Neighbours[5] = Vertneigh[index[(i+1)%neighbourhood_size]]; 
                  Neighbours[10] = Vertneigh[index[(i>0)?(i-1)%neighbourhood_size:neighbourhood_size-1]]; 
                  temp_reconstruct_index[4] = index[i]; temp_reconstruct_index[5] = index[(i+1)%neighbourhood_size];
                  intersection_compvar[2] = (point_interesct(0) + point_interesct(1) - coeffx[0] - coeffy[0])/(coeffx[2] + coeffy[2]);   

                  dist_intersection[2] = dist;
                  neighbour_weights[4] = dist2 / (dist1 + dist2);
                  neighbour_weights[5] = dist1 / (dist1 + dist2);

                  temp_i = i;
                  temp_j = (i>0)?i-1:neighbourhood_size-1;
                 }

                 // centermass line 1 a1*x+b1*y + c1 = 0
                dxv = neighCentermass[index[temp_j]](0) - neighCentermass[index[temp_i]](0);
                dyv = neighCentermass[index[temp_j]](1) - neighCentermass[index[temp_i]](1);
                temp_linev1[0] = -dyv; temp_linev1[1] = dxv; temp_linev1[2] = neighCentermass[index[temp_i]](0)*dyv - neighCentermass[index[temp_i]](1)*dxv;

                // centermass line 2 a2*x+b2*y + c2 = 0
                temp_linev2[0] = -V[1](1); temp_linev2[1] = V[1](0); temp_linev2[2] = Cellcentermass(0) * V[1](1) - Cellcentermass(1) * V[1](0);

                detv1 = temp_linev1[0]*temp_linev2[1] - temp_linev1[1]*temp_linev2[0];

                 mPoint point_interesctv1((temp_linev1[1]*temp_linev2[2] - temp_linev1[2]*temp_linev2[1])/detv1, (temp_linev1[2]*temp_linev2[0] - temp_linev1[0]*temp_linev2[2])/detv1);
                distv = Cellcentermass.pointdistance(point_interesctv1);
                distv1 = point_interesctv1.pointdistance(neighCentermass[index[temp_i]]);
                distv2 = point_interesctv1.pointdistance(neighCentermass[index[temp_j]]);
                distv3 = neighCentermass[index[temp_i]].pointdistance(neighCentermass[index[temp_j]]);
                
                dist_intersection_cm[2] = distv;
                neighbour_weights_cm[4] = distv2 / (distv1 + distv2);
                neighbour_weights_cm[5] = distv1 / (distv1 + distv2);


              }

              //if(0.64 < centroid(0) && centroid(0) < 0.651 &&  0.613 < centroid(1) && centroid(1) < 0.621){
                            //}
              //if(0.043 < centroid(0) && centroid(0) < 0.047 &&  -0.095 < centroid(1) && centroid(1) < -0.091621){
              //if(0.64 < centroid(0) && centroid(0) < 0.651 &&  0.613 < centroid(1) && centroid(1) < 0.621){
                //printf("Center mass intersection distances : %.12e, %.12e, %.12e\n",distv1, distv2, distv3);
                //printf("Center mass intersection distances : %.12e, %.12e, %.12e\n",distv1, distv2, distv3);
                //printf("Center mass : (%.12e,%.12e), Neighbouring  center mass : (%.12e,%.12e), (%.12e,%.12e)\n ",Cellcentermass(0), Cellcentermass(1), neighCentermass[index[i]](0),neighCentermass[index[i]](1), neighCentermass[index[(i+1)%neighbourhood_size]](0), neighCentermass[index[(i+1)%neighbourhood_size]](1));
                //printf("Centeroid : (%.12e,%.12e), Neighbouring  centeroid : (%.12e,%.12e), (%.12e,%.12e)\n ",centroid(0), centroid(1), neighCenter[index[i]](0),neighCenter[index[i]](1), neighCenter[index[(i+1)%neighbourhood_size]](0), neighCenter[index[(i+1)%neighbourhood_size]](1));

              //}
               
               
               for(int j = 4; j < 6; ++j){
                  for(int k=0;k<theMeshEntity->size(0);k++)
                  {
                    v_neigh[k] = ((Vertex*)Neighbours[j]->theMeshEntity->get(0,k))->point();
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
                  Jaccoeff_neigh[2] = coeff_neighx[3]*coeff_neighy[2] - coeff_neighx[2]*coeff_neighy[3];
                  coeff_alpha[j] = (1./Jaccoeff_neigh[0]) * (coeff_neighy[2]* (coeffx[1] + coeffx[3]*intersection_compvar[2]) - coeff_neighx[2] * (coeffy[1] + coeffy[3]*intersection_compvar[2]) );
                  coeff_beta[j] = (1./Jaccoeff_neigh[0]) * (-coeff_neighy[1]* (coeffx[1] + coeffx[3]*intersection_compvar[2]) + coeff_neighx[1] * (coeffy[1] + coeffy[3]*intersection_compvar[2]) );
               }
             }
             else if( (point_interesct(0) - centroid(0))*V[1](0) + (point_interesct(1) - centroid(1))*V[1](1) < 0 &&  dist > dist_intersection[3]){

               mPoint point_interesctv((linev1[1]*linev2[2] - linev1[2]*linev2[1])/detv, (linev1[2]*linev2[0] - linev1[0]*linev2[2])/detv);
              distv = Cellcentermass.pointdistance(point_interesctv);
              distv1 = point_interesctv.pointdistance(neighCentermass[index[i]]);
              distv2 = point_interesctv.pointdistance(neighCentermass[index[(i+1)%neighbourhood_size]]);
              distv3 = neighCentermass[index[i]].pointdistance(neighCentermass[index[(i+1)%neighbourhood_size]]);
              //if(0.043 < centroid(0) && centroid(0) < 0.047 &&  -0.095 < centroid(1) && centroid(1) < -0.091621){
              //if(0.64 < centroid(0) && centroid(0) < 0.651 &&  0.613 < centroid(1) && centroid(1) < 0.621){
              //  printf("Center mass intersection distances : %.12e, %.12e, %.12e\n",distv1, distv2, distv3);
              //}
               
              if(distv1 + distv2 < distv3+1e-10)
              {
                 Neighbours[6] = Vertneigh[index[i]]; Neighbours[7] = Vertneigh[index[(i+1)%neighbourhood_size]];
                 Neighbours[11] = Vertneigh[index[(i+1)%neighbourhood_size]];
                 temp_reconstruct_index[6] = index[i]; temp_reconstruct_index[7] = index[(i+1)%neighbourhood_size]; 
                 intersection_compvar[3] = (point_interesct(0) + point_interesct(1) - coeffx[0] - coeffy[0])/(coeffx[2] + coeffy[2]);

                 dist_intersection[3] = dist;
                 neighbour_weights[6] = dist2 / (dist1 + dist2);
                 neighbour_weights[7] = dist1 / (dist1 + dist2);
               
              
                 dist_intersection_cm[3] = distv;
                 neighbour_weights_cm[6] = distv2 / (distv1 + distv2);
                 neighbour_weights_cm[7] = distv1 / (distv1 + distv2);
               
               }
              else{
                int temp_i = i, temp_j = (i+1)%neighbourhood_size;
                 if (distv1 > distv2){

                  Neighbours[7] = Vertneigh[index[i]]; Neighbours[6] = Vertneigh[index[(i+1)%neighbourhood_size]]; 
                  Neighbours[11] = Vertneigh[index[(i+2)%neighbourhood_size]]; 
                  temp_reconstruct_index[7] = index[i]; temp_reconstruct_index[6] = index[(i+1)%neighbourhood_size];
                  intersection_compvar[3] = (point_interesct(0) + point_interesct(1) - coeffx[0] - coeffy[0])/(coeffx[2] + coeffy[2]);   

                  dist_intersection[3] = dist;
                  neighbour_weights[6] = dist1 / (dist1 + dist2);
                  neighbour_weights[7] = dist2 / (dist1 + dist2);

                  temp_i = temp_j;
                  temp_j = (i+2)%neighbourhood_size;   
                 }
                 else{

                  Neighbours[6] = Vertneigh[index[i]]; Neighbours[7] = Vertneigh[index[(i+1)%neighbourhood_size]]; 
                  Neighbours[11] = Vertneigh[index[(i>0)?(i-1)%neighbourhood_size:neighbourhood_size-1]]; 
                  temp_reconstruct_index[6] = index[i]; temp_reconstruct_index[7] = index[(i+1)%neighbourhood_size];
                  intersection_compvar[3] = (point_interesct(0) + point_interesct(1) - coeffx[0] - coeffy[0])/(coeffx[2] + coeffy[2]);   

                  dist_intersection[3] = dist;
                  neighbour_weights[6] = dist2 / (dist1 + dist2);
                  neighbour_weights[7] = dist1 / (dist1 + dist2);

                  temp_i = i;
                  temp_j = (i>0)?i-1:neighbourhood_size-1;
                 }

                 // centermass line 1 a1*x+b1*y + c1 = 0
                dxv = neighCentermass[index[temp_j]](0) - neighCentermass[index[temp_i]](0);
                dyv = neighCentermass[index[temp_j]](1) - neighCentermass[index[temp_i]](1);
                temp_linev1[0] = -dyv; temp_linev1[1] = dxv; temp_linev1[2] = neighCentermass[index[temp_i]](0)*dyv - neighCentermass[index[temp_i]](1)*dxv;

                // centermass line 2 a2*x+b2*y + c2 = 0
                temp_linev2[0] = -V[1](1); temp_linev2[1] = V[1](0); temp_linev2[2] = Cellcentermass(0) * V[1](1) - Cellcentermass(1) * V[1](0);

                detv1 = temp_linev1[0]*temp_linev2[1] - temp_linev1[1]*temp_linev2[0];

                 mPoint point_interesctv1((temp_linev1[1]*temp_linev2[2] - temp_linev1[2]*temp_linev2[1])/detv1, (temp_linev1[2]*temp_linev2[0] - temp_linev1[0]*temp_linev2[2])/detv1);
                distv = Cellcentermass.pointdistance(point_interesctv1);
                distv1 = point_interesctv1.pointdistance(neighCentermass[index[temp_i]]);
                distv2 = point_interesctv1.pointdistance(neighCentermass[index[temp_j]]);
                distv3 = neighCentermass[index[temp_i]].pointdistance(neighCentermass[index[temp_j]]);
                
                dist_intersection_cm[3] = distv;
                neighbour_weights_cm[6] = distv2 / (distv1 + distv2);
                neighbour_weights_cm[7] = distv1 / (distv1 + distv2);
              }
                
               
               
               for(int j = 6; j < 8; ++j){
                  for(int k=0;k<theMeshEntity->size(0);k++)
                  {
                    v_neigh[k] = ((Vertex*)Neighbours[j]->theMeshEntity->get(0,k))->point();
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
                  Jaccoeff_neigh[2] = coeff_neighx[3]*coeff_neighy[2] - coeff_neighx[2]*coeff_neighy[3];
                  coeff_alpha[j] = (1./Jaccoeff_neigh[0]) * (coeff_neighy[2]* (coeffx[1] + coeffx[3]*intersection_compvar[3]) - coeff_neighx[2] * (coeffy[1] + coeffy[3]*intersection_compvar[3]) );
                  coeff_beta[j] = (1./Jaccoeff_neigh[0]) * (-coeff_neighy[1]* (coeffx[1] + coeffx[3]*intersection_compvar[3]) + coeff_neighx[1] * (coeffy[1] + coeffy[3]*intersection_compvar[3]) );
               }
             }
           }
        }


      }



      // Finding interpolation points using the reconstruction stencil made of cell center mass
     // dist_intersection_cm[0] = dist_intersection_cm[1] = scaling_derivative[0];
     // dist_intersection_cm[2] = dist_intersection_cm[3] = scaling_derivative[1];

     /* for(int i = 0; i < neighbourhood_size; ++i){
        // line 1 a1*x+b1*y + c1 = 0
        dxv = neighCentermass[index[(i+1)%neighbourhood_size]](0) - neighCentermass[index[i]](0);
        dyv = neighCentermass[index[(i+1)%neighbourhood_size]](1) - neighCentermass[index[i]](1);
        linev1[0] = -dyv; linev1[1] = dxv; linev1[2] = neighCentermass[index[i]](0)*dyv - neighCentermass[index[i]](1)*dxv;

        // line 2 a2*x+b2*y + c2 = 0
        linev2[0] = -V[0](1); linev2[1] = V[0](0); linev2[2] = Cellcentermass(0) * V[0](1) - Cellcentermass(1) * V[0](0);

        detv = line1[0]*line2[1] - line1[1]*line2[0];
        if (sqrt(detv*detv) > 1e-10){
           mPoint point_interesctv((linev1[1]*linev2[2] - linev1[2]*linev2[1])/detv, (linev1[2]*linev2[0] - linev1[0]*linev2[2])/detv);
           distv = Cellcentermass.pointdistance(point_interesctv);
           distv1 = point_interesctv.pointdistance(neighCentermass[index[i]]);
           distv2 = point_interesctv.pointdistance(neighCentermass[index[(i+1)%neighbourhood_size]]);
           distv3 = neighCentermass[index[i]].pointdistance(neighCentermass[index[(i+1)%neighbourhood_size]]);
           if(distv1 + distv2 < distv3+1e-10){
             if ( (point_interesct(0) - Cellcentermass(0))*V[0](0) + (point_interesct(1) - Cellcentermass(1))*V[0](1) > 0  && dist > dist_intersection_cm[0] ) {
               Neighbours_cm[0] = Vertneigh[index[i]]; Neighbours_cm[1] = Vertneigh[index[(i+1)%neighbourhood_size]]; 
               dist_intersection_cm[0] = dist;
               neighbour_weights_cm[0] = dist2 / (dist1 + dist2);
               neighbour_weights_cm[1] = dist1 / (dist1 + dist2); 
                }
             else if( (point_interesct(0) - Cellcentermass(0))*V[0](0) + (point_interesct(1) - Cellcentermass(1))*V[0](1) < 0  && dist > dist_intersection_cm[1]) {
               Neighbours_cm[2] = Vertneigh[index[i]]; Neighbours_cm[3] = Vertneigh[index[(i+1)%neighbourhood_size]]; 
               dist_intersection_cm[1] = dist;
               neighbour_weights_cm[2] = dist2 / (dist1 + dist2);
               neighbour_weights_cm[3] = dist1 / (dist1 + dist2); 
              }
           }
        }


        // line 2 a2*x+b2*y + c2 = 0
        linev2[0] = -V[1](1); linev2[1] = V[1](0); linev2[2] = Cellcentermass(0) * V[1](1) - Cellcentermass(1) * V[1](0);

        detv = linev1[0]*linev2[1] - linev1[1]*linev2[0];
        if (sqrt(detv*detv) >  1e-10){
           mPoint point_interesct((line1[1]*line2[2] - line1[2]*line2[1])/det, (line1[2]*line2[0] - line1[0]*line2[2])/det);
           dist = Cellcentermass.pointdistance(point_interesct);
           dist1 = point_interesct.pointdistance(neighCentermass[index[i]]);
           dist2 = point_interesct.pointdistance(neighCentermass[index[(i+1)%neighbourhood_size]]);
           dist3 = neighCentermass[index[i]].pointdistance(neighCentermass[index[(i+1)%neighbourhood_size]]);
           if(dist1 + dist2 < dist3+1e-10){
             if ( (point_interesct(0) - Cellcentermass(0))*V[1](0) + (point_interesct(1) - Cellcentermass(1))*V[1](1) > 0 && dist > dist_intersection_cm[2] ) {
               Neighbours_cm[4] = Vertneigh[index[i]]; Neighbours_cm[5] = Vertneigh[index[(i+1)%neighbourhood_size]];
               dist_intersection_cm[2] = dist;
               neighbour_weights_cm[4] = dist2 / (dist1 + dist2);
               neighbour_weights_cm[5] = dist1 / (dist1 + dist2);
             }
             else if( (point_interesct(0) - Cellcentermass(0))*V[1](0) + (point_interesct(1) - Cellcentermass(1))*V[1](1) < 0 &&  dist > dist_intersection_cm[3]){
               Neighbours_cm[6] = Vertneigh[index[i]]; Neighbours_cm[7] = Vertneigh[index[(i+1)%neighbourhood_size]];
               dist_intersection_cm[3] = dist;
               neighbour_weights_cm[6] = dist2 / (dist1 + dist2);
               neighbour_weights_cm[7] = dist1 / (dist1 + dist2);
             }
           }
        }
      }*/
      
      // Copying the orignial reconstruction stencil over the cm reconstruction stencil
      for(int i = 0;i < 8; ++i){
      //  Neighbours_cm[i] = Neighbours[i];
       // neighbour_weights_cm[i] = neighbour_weights[i];
      }

      

    // Checking if the neighoburing cells have been identified?
    int id_flag = 0;
    for(int i =0 ; i< 8; ++i){
      if(Neighbours[i] != NULL) id_flag++;
    }
    //if(!id_flag) printf("Neigbourhood size : %d, center of mass L (%.12e,%.12e)\n",neighbourhood_size, centroid(0),centroid(1));

    //if( -0.2 < centroid(0) && centroid(0) < 0.2 && -.2 < centroid(1) && centroid(1) < 0.2)
   /* if( -0.2 < centroid(0) && centroid(0) < -0.15 && -.2 < centroid(1) && centroid(1) < -0.15)
        {
      //if(id_flag == 6)
      // {
      //  printf("Neigbourhood size : %d, %d,  center of mass L (%.12e,%.12e)\n",neighbourhood_size, id_flag, centroid(0),centroid(1));
        //printf("Recosntruction indices : ");
        //for(int i = 0; i < 8; ++i){
        //  printf("%d ",temp_reconstruct_index[i]);
        //}
        //printf("\n");
      //}

      printf("No of neighbours : %d\n Centers : ", id_flag);
      for(int i =0 ; i < 8; ++i){
        printf("(%.12e, %.12e),  [%.12e],  ", neighCenter[temp_reconstruct_index[i]](0),neighCenter[temp_reconstruct_index[i]](1), neighbour_weights[i]);
      }
      printf("\n");
      printf("Ordeed neihgbour centroids : ");
      for(int i= 0;i< neighbourhood_size;++i){
        printf("(%.12e, %.12e), ", neighCenter[index[i]](0),neighCenter[index[i]](1));
      }
      printf("\n");
      printf("Ordeed neihgbour angles : ");
      for(int i= 0;i< neighbourhood_size;++i){
        printf("(%f), ", angle[i]*(180/PI));
      }
      printf("\n");
      //for(int i =0 ; i< neighbourhood_size; ++i){
      //  printf("Neighbour %d center : (%.12e,%.12e)\n",i,neighCenter[index[i]](0),neighCenter[index[i]](1));
     // }
    }*/
    //if(neighbourhood_size >= 20) printf("Element with insufficent memeory!!\n");

  //if( -0.14 < centroid(0) && centroid(0) < -0.13 &&  0.16 < centroid(1) && centroid(1) < 0.167){
	//	if( -0.15 < centroid(0) && centroid(0) < -0.14 &&  0.16 < centroid(1) && centroid(1) < 0.167){
		//if( -0.284 < centroid(0) && centroid(0) < -0.26 &&  -0.8 < centroid(1) && centroid(1) < -0.76){
			/*	if( -0.167 < centroid(0) && centroid(0) < -0.155 &&  -0.175 < centroid(1) && centroid(1) < -0.167){
        printf("c3 coeff : %.12e, %.12e, %.12e \n", coeffx[3], coeffy[3],sqrt(coeffx[3]*coeffx[3] + coeffy[3]*coeffy[3] ));
        printf("c2 coeff : %.12e, %.12e, %.12e \n", coeffx[2], coeffy[2],sqrt(coeffx[2]*coeffx[2] + coeffy[2]*coeffy[2] ));
        printf("c1 coeff : %.12e, %.12e, %.12e \n", coeffx[1], coeffy[1], sqrt(coeffx[1]*coeffx[1] + coeffy[1]*coeffy[1] ));
        //printf("Neigbourhood size : %d, %d, center of mass L (%.12e,%.12e)\n",id_flag, neighbourhood_size, centroid(0),centroid(1));
         printf("Centroid : (%.12e,%.12e)\n",centroid(0),centroid(1));
         printf("Jacobian : %.12e, %.12e, %.12e\n",Jaccoeff[0], Jaccoeff[1], Jaccoeff[2]);
				//printf("Intersection ref elem variables : %.12e, %.12e, %.12e, %.12e\n",intersection_compvar[0],intersection_compvar[1],intersection_compvar[2],intersection_compvar[3]);
				}*/
       
    return Jac_flag;
		
}
/************ B O U N D A R Y   T E R M S *************/

void DGBoundaryCell::check()
{
  mEntity *b[2] = {0,0};
  int k = 0;
  int n = theBoundaryEntity->getLevel();
  
  if(!theBoundaryEntity->isAdjacencyCreated(n+1))
    {
      printf("weird face because no upward adj :");theBoundaryEntity->print();
      if(theBoundaryEntity->parent())
	{
	  printf("this face has parent : ");
	  theBoundaryEntity->parent()->print();
	}      
      if(theBoundaryEntity->isAdjacencyCreated(2))
	{
	  printf("this face has childeren and should not be there\n ");
	}      
    }
mUpwardAdjacencyContainer::iter it;
mUpwardAdjacencyContainer::iter end_np1 = theBoundaryEntity->end(n+1);

  for( it= theBoundaryEntity->begin(n+1); it !=  end_np1; ++it)
    {
      b[k++] = (*it);
    }

  if(mleft != b[0] || mright != b[1])
    { 
      mleft = b[0];
      mright = b[1];
      right = 0;
      left = (DGCell*)mleft->getCell();
      if(mright)right = (DGCell*)mright->getCell();
      init();
    }
}

DGBoundaryCell::DGBoundaryCell (mEntity *ent, Mapping *m,GaussIntegrator *i)
  : theBoundaryEntity(ent),mleft(0),mright(0),theMapping(m),
    theGaussIntegrator(i) 
{
  left=right=0;
  check();
}


DGBoundaryCell::~DGBoundaryCell ()
{
  delete theMapping;
  for(int i=0;i<nbPtGauss;++i) delete boundaryGaussPoints[i];
}

int DGBoundaryCell::computeOrder() const
{
  if(right) return 
      2 * ((left->fOrder>right->fOrder) ?
	  left->fOrder : right->fOrder) +  theMapping->order();
  else return 2*(left->fOrder) + theMapping->order();
}


boundaryGaussPoint::boundaryGaussPoint()
{}

boundaryGaussPoint::~boundaryGaussPoint()
{}

void DGBoundaryCell::init()
{
  int i,j,k,qleft,qright;
  int order = computeOrder()-1;
  cSize = left->theConservationLaw->getNbFields();
  lSize = left->fSize;
  if (right) rSize = right->fSize; else rSize = 0;
  //theBoundaryEntity->print();
  if(boundaryGaussPoints.size())
    {
      for(i=0;i<boundaryGaussPoints.size();++i)delete boundaryGaussPoints[i];
      boundaryGaussPoints.clear();
    }

  nbPtGauss = theGaussIntegrator->nbIntegrationPoints(theBoundaryEntity,order);
  boundaryGaussPoints.reserve(nbPtGauss);
  double u,v,w,weight,uleft,vleft,wleft,urght,vrght,wrght;
   
  for(j=0;j<left->theMeshEntity->size(0);j++)
	{
	  Vertex* vV = (Vertex*) left->theMeshEntity->get(0,j);
	  int found = 0;
	  for(k=0;k<theBoundaryEntity->size(0);k++) 
	    {
	      Vertex* vB = (Vertex*) theBoundaryEntity->get(0,k);
	      if (vV==vB) {found = 1; break;}
	    }
	  if (found == 0) {qleft=j-1; if (qleft<0) qleft=3; break;} 
	}
  if (right) {
     for(j=0;j<right->theMeshEntity->size(0);j++)
	{
	  Vertex* vV = (Vertex*) right->theMeshEntity->get(0,j);
	  int found = 0;
	  for(k=0;k<theBoundaryEntity->size(0);k++) 
	    {
	      Vertex* vB = (Vertex*) theBoundaryEntity->get(0,k);
	      if (vV==vB) {found = 1; break;}
	    }
	  if (found == 0) {qright=j-1; if (qright<0) qright=3; break;} 
	}
  }

  for(i=0;i<nbPtGauss;++i)
    {
      boundaryGaussPoint *pg = new boundaryGaussPoint();
      theGaussIntegrator->iPoint(theBoundaryEntity,i,order,u,v,w,weight);
      // compute parametric coordinates on both sides
	 
      theMapping->eval(u,v,w,pg->x,pg->y,pg->z);
      
      if(!left->theMapping->invert(pg->x,pg->y,pg->z,uleft,vleft,wleft))
	{
	  printf("l : impossible to compute invert mapping in ");theBoundaryEntity->print();
	  for(int ip=0;ip<theBoundaryEntity->size(0);ip++)
	    theBoundaryEntity->get(0,ip)->print();
	}
      
      if(right)
	if(!right->theMapping->invert(pg->x,pg->y,pg->z,urght,vrght,wrght))
	  {
	    printf("r : impossible to compute invert mapping in ");theBoundaryEntity->print();
	    for(int ip=0;ip<theBoundaryEntity->size(0);ip++)
	      theBoundaryEntity->get(0,ip)->print();
	  }
      
	
      // get the normal vector
      
      if(left->theMeshEntity->find(theBoundaryEntity))
	left->theMapping->normalVector(theBoundaryEntity,uleft,vleft,wleft,n);
      else
	{
	  //	  printf("computing normal in the complex way  :  ");
	  bool found = false;
	  for(j=0;j<mleft->size(theBoundaryEntity->getLevel());++j)
	    {
	      mEntity *initial = mleft->get(theBoundaryEntity->getLevel(),j); 
	      list<mEntity*> leaves;	  
	      initial->getLeaves(leaves);
	      for(list<mEntity*>::const_iterator it = leaves.begin();it != leaves.end();++it)
		{
		  if(*it == theBoundaryEntity)
		    {
		      left->theMapping->normalVector(initial,uleft,vleft,wleft,n);
		     		     // printf("%f %f %f %f\n",n(0),n(1),pg->x,pg->y);
		      found = true;
		      break;
		    }
		}
	    }
	  if(!found)
	    {
	      printf("unable to compute normal vector\n");
	    }
	}
	  
      // get det of jacobian
      detJac = theMapping->detJac(u,v,w);
      pg->JacTimesWeight = detJac * weight;
      // get shape functions on both sides
      pg->fctleft = new double[lSize];
      left->theFunctionSpace->fcts(uleft,vleft,wleft,pg->fctleft);
      //printf("Nb of Gauss points= %d, nb of functions= %d \n", nbPtGauss,lSize); 
    /*  printf("{");
      for (q=0; q<lSize; q++) 
	if (q<lSize-1) printf("%17.16e, ", pg->fctleft[q]);
	else printf("%17.16e}, \n", pg->fctleft[q]);*/
      
	 double* ff=0;
	 ff = left->theFunctionSpace->fcts(left->fOrder,qleft,i);
     //pg->fctleft = left->theFunctionSpace->fcts(left->fOrder,qleft,i);
	//for(j=0;j<rSize; j++) if (fabs (pg->fctleft[j]-ff[j])>0.000000001) printf (" left %d  %d %d %e %e \n",q, i, j, pg->fctleft[j],ff[j]);
	//for(q=0;q<2; q++) printf (" left %e %e \n",pg->fctleft[q],ff[q]);
	//printf("\n");
      if(right)
	{
		pg->fctrght = new double [rSize];
	 right->theFunctionSpace->fcts(urght,vrght,wrght,pg->fctrght); 
	
      
     //pg->fctrght = right->theFunctionSpace->fcts(right->fOrder,qright,i);
	//  ff = right->theFunctionSpace->fcts(right->fOrder,qright,i);
	 // for(j=0;j<rSize; j++) if (fabs (pg->fctrght[j]-ff[j])>0.000000001) printf (" right %d %d %d %e %e \n",q, i, j, pg->fctrght[j],ff[j]);
    /*for(q=0;q<2; q++) printf (" right %e %e \n",pg->fctright[q],ff[q]);
	printf("\n");*/
	}
      boundaryGaussPoints.push_back(pg);
    }
  computeSize();
  //printf("\n");
}

void DGBoundaryCell::setPeriodicBC()
{
  mEntity* ent=theBoundaryEntity->getAttachedEntity(TYPE_CELL);
  DGBoundaryCell *symm = (DGBoundaryCell *)ent->getCell();
  DGCell* another = symm->left; 
  rSize = another->fSize; 	
  if (theBoundaryEntity->getClassification()->getId()==610 || theBoundaryEntity->getClassification()->getId()==510)
    for(int i=0;i<nbPtGauss;++i)
      {
	double urght,vrght,wrght;
	boundaryGaussPoint *pg = pt(i);      
	mPoint p(pg->x,pg->y,pg->z);
	if (theBoundaryEntity->getClassification()->getId()==510)
	  another->theMapping->invert(-pg->x,pg->y,pg->z,urght,vrght,wrght);
	if (theBoundaryEntity->getClassification()->getId()==610)
	  another->theMapping->invert(pg->x,-pg->y,pg->z,urght,vrght,wrght);

	pg->fctrght = new double [rSize];
	another->theFunctionSpace->fcts(urght,vrght,wrght,pg->fctrght); 
	}
}

void DGBoundaryCell::computeError()
{
  if(!mright)return;
  double valleft[MaxNbEqn],valrght[256];
  double error = 0;
  boundaryGaussPoint *pg;
  for(int i=0;i<nbPtGauss;++i)
    {
	pg = pt(i);
      left  ->interpolate(pg->fctleft,valleft);
      right ->interpolate(pg->fctrght,valrght);
      double jump = left->theConservationLaw->jumpQuantity(valleft)
	- right->theConservationLaw->jumpQuantity(valrght); 
      error += (jump*jump) * pg->JacTimesWeight;
    }
  left->error += error;
  right->error += error;
}

void DGBoundaryCell::computeJump()
{
/*  if(!mright)return;
  DGCell *left = (DGCell*)mleft->getCell();
  DGCell *right = (DGCell*)mright->getCell();
  int order = computeOrder(left,right);
  double u,v,w,weight,valleft[MaxNbEqn],valrght[MaxNbEqn];
  double jump = 0;
  mVector n;
  
  for(int i=0;i< theGaussIntegrator.nbIntegrationPoints(theBoundaryEntity,order);++i)
    {
      theGaussIntegrator.iPoint(theBoundaryEntityi,order,u,v,w,weight);
      double detJac = theMapping->detJac(u,v,w);

      left  ->interpolate(uleft[i],vleft[i],wleft[i],valleft);
	 
	  right ->interpolate(urght[i],vrght[i],wrght[i],valrght);
      //printf("%e %e \n", valleft, valrght);
  if(left->theMeshEntity->find(theBoundaryEntity))
	left->theMapping->normalVector(theBoundaryEntity,uleft_e[i],vleft_e[i],wleft_e[i],n);
     else
	{
	  if(!right)printf("aaargh\n");
	  right->theMapping->normalVector(theBoundaryEntity,urght_e[i],vrght_e[i],wrght_e[i],n);
	  n *= -1.0;
	}
      mPoint p;
      theMapping->eval(u,v,w,p(0),p(1),p(2));
      int orient = left->theConservationLaw->getEdgeOrientation(n,p,valleft);
	  if (orient == 1) 
	  {
		  left->jump+= (valrght[0] - valleft[0]) * detJac * weight;
	  }
	  else if(right) right->jump +=(valrght[0] - valleft[0]) * detJac * weight;
	  //printf(" %e %e %e %e %e \n", valrght[0], valleft[0], detJac, weight,(valrght - valleft) * detJac * weight);
	  
  }
  */
}

int DGBoundaryCell::computeBoundaryContributions(double T)
{
  int i,j,k;
  int nonPhysical = 0;
  int temp = 0;
  double riemann[MaxNbEqn];
  boundaryGaussPoint *pg;
  double valleft[MaxNbEqn],valrght[MaxNbEqn];
   
  for(i=0;i<nbPtGauss;++i)
    {
      pg = pt(i);      
		//	printf("%f %f\n",pg->N[0],pg->N[1]);
      mPoint p(pg->x,pg->y,pg->z);
      left ->interpolate(pg->fctleft, valleft);
      if(!left->theConservationLaw->isPhysical(valleft))
	{
	  //printf("non physical field found on an EDGE/FACE at point (%f,%f,%f)\n",pg->x,pg->y,pg->z); 
	  double p;
	  int dim = left->theMeshEntity->getLevel();
	  if (dim==2)
	    p = (0.4)*(valleft[3]-0.5*(valleft[1]*valleft[1]+valleft[2]*valleft[2])/valleft[0]);
	  else 
	    p = (0.4)*(valleft[3]-0.5*(valleft[1]*valleft[1]+valleft[2]*valleft[2]+valleft[4]*valleft[4])/valleft[0]);
	  //printf("rho=%e p=%e \n",valleft[0],p);
	  //printf("Boundary Id =  %d \n",theBoundaryEntity->getClassification()->getId());
 	  nonPhysical =1;

     for(j=0;j<cSize;++j)
     for(k=1;k<lSize;++k) 
      left->theFieldsCoefficients->get(j,k) = 0.0;

      left ->interpolate(pg->fctleft, valleft);
	}

      if (right)
	  {
	right ->interpolate(pg->fctrght, valrght);
	
	if(!right->theConservationLaw->isPhysical(valrght))
	  {
	    //printf("non physical field found (right) on an edge at point (%f,%f,%f)\n",pg->x,pg->y,pg->z); 
	    double p;
	    int dim = left->theMeshEntity->getLevel();
	    if (dim==2)
	      p= (0.4)*(valrght[3]-0.5*(valrght[1]*valrght[1] +valrght[2]*valrght[2])/valrght[0]);
	    else
	      p= (0.4)*(valrght[3]-0.5*(valrght[1]*valrght[1] +valrght[2]*valrght[2]+valrght[4]*valrght[4])/valrght[0]);
	    //printf("rho = %e p = %e \n",valrght[0],p);
	    //printf("Id =  %d \n",theBoundaryEntity->getClassification()->getId());
	    if (nonPhysical==1) nonPhysical = 3;
	    else nonPhysical = 2;

      for(j=0;j<cSize;++j)
     for(k=1;k<rSize;++k) 
      right->theFieldsCoefficients->get(j,k) = 0.0;

      right ->interpolate(pg->fctrght, valrght);

	  }
	  }

    nonPhysical = 0;
  //Periodic boundary conditions
  if (theBoundaryEntity->getClassification()->getId()==610 || theBoundaryEntity->getClassification()->getId()==510)
    {
      mEntity* ent=theBoundaryEntity->getAttachedEntity(TYPE_CELL);
      DGBoundaryCell *symm = (DGBoundaryCell *)ent->getCell();
      DGCell* another = symm->left; 
	  another->interpolate(pg->fctrght, valrght);
    }	

      if(right)      //inner cell
	left->theConservationLaw->
	  riemannSolver (n,p,valleft,valrght,riemann);
      // or boundary conditions
      else
	if(theBoundaryEntity->getClassification()->getId()==20000) 
	  left->theConservationLaw->boundary (n,pg->N,theBoundaryEntity->getClassification()->getId(),  p, valleft,riemann,T);
	else if (theBoundaryEntity->getClassification()->getId()==510 || theBoundaryEntity->getClassification()->getId()==610) 
	  left->theConservationLaw->riemannSolver(n,p,valleft,valrght,riemann);
	else left->theConservationLaw->boundaryFlux (n,theBoundaryEntity->getClassification()->getId(),  p, valleft,riemann,T);
      
	//printf("%f %f %f %f %f \n",riemann[0],riemann[1],riemann[2],riemann[3],riemann[4]); 
      // left and right contributions
      double *RHS_LEFT = left->theRightHandSide;
      const double *RIEMANN = riemann;
      const double *FUNCTION_LEFT  = pg->fctleft;
      for(j=0;j<cSize;++j)
	for(k=0;k<lSize;++k)
	  (*(RHS_LEFT++)) -= RIEMANN[j]*FUNCTION_LEFT[k]*pg->JacTimesWeight;
       
      if (right)
	{
	  double  *RHS_RIGHT = right->theRightHandSide;
	  const double *FUNCTION_RIGHT  = pg->fctrght;  
	  for(j=0;j<cSize;++j)
	    for(k=0;k<rSize;++k)
	      (*(RHS_RIGHT++)) += RIEMANN[j]*FUNCTION_RIGHT[k]*pg->JacTimesWeight;
	}
    }
  return nonPhysical;
}

void DGBoundaryCell::reverseBoundaryContributions(double T)
{
  int i,j,k;
  double riemann[MaxNbEqn];
  boundaryGaussPoint *pg;
  double valleft[MaxNbEqn],valrght[MaxNbEqn];
  
  for(i=0;i<nbPtGauss;++i)
    {
      pg= pt(i);      
      left ->interpolate(pg->fctleft, valleft);
      if(right) right ->interpolate(pg->fctrght, valrght);
      
      //Periodic boundary conditions
      if (theBoundaryEntity->getClassification()->getId()==610 || theBoundaryEntity->getClassification()->getId()==510)
	{
	  mEntity* ent=theBoundaryEntity->getAttachedEntity(TYPE_CELL);
	  DGBoundaryCell *symm = (DGBoundaryCell *)ent->getCell();
	  DGCell* another = symm->left; 
	  another->interpolate(pg->fctrght, valrght);
	}
      
      mPoint p(pg->x,pg->y,pg->z);
      if(right)      //inner cell
	left->theConservationLaw->riemannSolver(n,p,valleft,valrght,riemann);
      // or boundary conditions
      else
	if(theBoundaryEntity->getClassification()->getId()==20000) 
	  left->theConservationLaw->boundary (n,pg->N,theBoundaryEntity->getClassification()->getId(),  p, valleft,riemann,T);
	else if (theBoundaryEntity->getClassification()->getId()==510 || theBoundaryEntity->getClassification()->getId()==610) 
	  left->theConservationLaw->riemannSolver(n,p,valleft,valrght,riemann);
	else left->theConservationLaw->boundaryFlux (n,theBoundaryEntity->getClassification()->getId(),  p, valleft,riemann,T);
      
      
      // left and right contributions
      double *RHS_LEFT = left->theRightHandSide;
      const double *RIEMANN = riemann;
      const double *FUNCTION_LEFT  = pg->fctleft;
      for(j=0;j<cSize;++j)
	for(k=0;k<lSize;++k)
	  (*(RHS_LEFT++)) += RIEMANN[j]*FUNCTION_LEFT[k]*pg->JacTimesWeight;
      
      if (right)
	{
	  double  *RHS_RIGHT = right->theRightHandSide;
	  const double *FUNCTION_RIGHT  = pg->fctrght;  
	  for(j=0;j<cSize;++j)
	    for(k=0;k<rSize;++k)
	      (*(RHS_RIGHT++)) -= RIEMANN[j]*FUNCTION_RIGHT[k]*pg->JacTimesWeight;
	}
    }
}

void DGCell::setToZero(double t)
{
  int i,j,k;
  int level = theMeshEntity->getLevel();
  int nbBound = theMeshEntity->size(level-1);
  mEntity* ent;
  for(i=0;i<nbBound;i++) 
    {
      ent = theMeshEntity->get(level-1,i);
      DGBoundaryCell* cell=(DGBoundaryCell *)ent->getCell();
      cell->reverseBoundaryContributions(t);
    }
  reverseVolumeContribution(t);

  for(j=0;j<cSize;++j)
    for(k=1;k<fSize;++k) 
      theFieldsCoefficients->get(j,k) = 0.0;

  computeVolumeContribution(t);

  for(i=0;i<nbBound;i++) 
    {
      ent = theMeshEntity->get(level-1,i);
      DGBoundaryCell* cell=(DGBoundaryCell *)ent->getCell();
      cell->computeBoundaryContributions(t);
    }
}

void DGBoundaryCell::computeSize()
{
  switch(theBoundaryEntity->getType())
    {
    case 0:
      size = 0.0;
      break;
    case 1:        /*Edge*/
      {
	Vertex* v1 = (Vertex* )theBoundaryEntity->get(0,0);
	Vertex* v2 = (Vertex* )theBoundaryEntity->get(0,1);
	mPoint p1 = v1->point();
	mPoint p2 = v2->point();
	size = sqrt((p1(0)-p2(0))*(p1(0)-p2(0))+(p1(1)-p2(1))*(p1(1)-p2(1))+(p1(2)-p2(2))*(p1(2)-p2(2)));
     	}
      break;
    case 2:        /*Triangle*/
      size = 0.5*detJac;
      //printf("Needs to be done");
      break;
    default: printf("no size for diemnsions greater than 2\n");
    }
}

 void DGBoundaryCell::normalToCircle(mPoint &p1, mPoint &p2,mPoint &p3)
{
  double x,y,r;
  double c1 = p1(0)*p1(0) + p1(1)*p1(1);
  double c2 = c1 - p2(0)*p2(0) - p2(1)*p2(1);
  double c3 = c1 - p3(0)*p3(0) - p3(1)*p3(1);
  double d1 = p1(0)-p3(0);
  double d2 = p1(0)-p2(0);
  double d3 = p1(1)-p3(1);
  double d4 = p1(1)-p2(1);
  double center_x = 0.5*(c2*d3-c3*d4)/(d2*d3-d1*d4);
  double center_y = 0.5*(c2*d1-c3*d2)/(d4*d1-d3*d2);
  for(int i=0;i<nbPtGauss;++i)
    {
      x = boundaryGaussPoints[i]->x - center_x;
      y = boundaryGaussPoints[i]->y - center_y;
      r = sqrt(x*x +y*y);
	  //printf("radius %e \n", r);
	   if (n(0)*x +n(1)*y>0)
	   {
      boundaryGaussPoints[i]->N(0) = x/r;  
      boundaryGaussPoints[i]->N(1) = y/r;  
      boundaryGaussPoints[i]->N(2) = 0.0;
	   } else 
	   {
      boundaryGaussPoints[i]->N(0) = -x/r;  
      boundaryGaussPoints[i]->N(1) = -y/r;  
      boundaryGaussPoints[i]->N(2) = 0.0;
	   }
    }
}

 void DGBoundaryCell::normalToCircle(double center_x, double center_y, double center_z)
{
  double x,y,z,r;
  //printf(" %e %e \n",center_x,center_y);
  for(int i=0;i<nbPtGauss;++i)
    {
      x = boundaryGaussPoints[i]->x - center_x;
      y = boundaryGaussPoints[i]->y - center_y;
      z = boundaryGaussPoints[i]->z - center_z;
      r = sqrt(x*x +y*y + z*z);
//	  printf("radius %e \n",r);
	  if (n(0)*x +n(1)*y + n(2)*z>0)
	  {
      boundaryGaussPoints[i]->N(0) = x/r;  
      boundaryGaussPoints[i]->N(1) = y/r;  
      boundaryGaussPoints[i]->N(2) = z/r;
	  } 
	  else
	  {
      boundaryGaussPoints[i]->N(0) = -x/r;  
      boundaryGaussPoints[i]->N(1) = -y/r;  
      boundaryGaussPoints[i]->N(2) = -z/r;
	  } 
    }
}

double DGBoundaryCell::computeRadius(mPoint &p1, mPoint &p2,mPoint &p3)
{
  double x,y;
  double c1 = p1(0)*p1(0) + p1(1)*p1(1);
  double c2 = c1 - p2(0)*p2(0) - p2(1)*p2(1);
  double c3 = c1 - p3(0)*p3(0) - p3(1)*p3(1);
  double d1 = p1(0)-p3(0);
  double d2 = p1(0)-p2(0);
  double d3 = p1(1)-p3(1);
  double d4 = p1(1)-p2(1);
  double center_x = 0.5*(c2*d3-c3*d4)/(d2*d3-d1*d4);
  double center_y = 0.5*(c2*d1-c3*d2)/(d4*d1-d3*d2);
  x = p1(0)-center_x;
  y = p1(1)-center_y;
  return sqrt(x*x+y*y);
 }

// Calculates the determinant of an nxn square matrix a
// Defined recursively using expansion by minors
double DGBoundaryCell::determinant(double **a, int n)
{
       double det=0;
       double **m=NULL;
       
       if (n==1)
       {
            det=a[0][0];   
       }
       else if(n==2)
       {
            det=a[0][0]*a[1][1]-a[1][0]*a[0][1];  
       }
       else if(n>2)
       {
            for(int k=0;k<n;k++)
            {
                    // Allocate space for minor
                    m=(double**)malloc((n-1)*sizeof(double *));
                    for(int i=0;i<n-1;i++)
                    {
                             // Allocate space for each column in the minor
                             m[i]=(double*)malloc((n-1)*sizeof(double));
                    }
                    // Create the minor matrix
                    for(int i=1;i<n;i++)
                    {
                            int c=0;
                            for (int j=0;j<n;j++)
                            {
                                if (j!=k)
                                {
                                         m[i-1][c]=a[i][j];
                                         c++;
                                }
                            }
                    }
                    // Calculate the determinant recursively
                    det+=pow(-1.,k)*a[0][k]*determinant(m,n-1);
                    // Free space used for the minor
                    for(int i=0; i<n-1;i++)
                    {
                            free(m[i]);
                    }
                    free(m);
            }
       }
       return det;
}

// Returns the radius of the sphere formed by four points.  If the four points
// are coplaner then there are zero or infinitly many solutions so return 0.
// The centre is returned as p[4].
double DGBoundaryCell::computeSphere(mPoint p[4])
{
       double r=0;
       
       // Allocate space for the minors
       double **a=(double**)malloc(4*sizeof(double *));
       for (int i=0;i<4;i++)
       {
           a[i]=(double*)malloc(4*sizeof(double));
       }
       
       // Find determinant M11
       for (int i=0;i<4;i++)
       {
           a[i][0]=p[i](0);
           a[i][1]=p[i](1);
           a[i][2]=p[i](2);
           a[i][3]=1;
       }
       double m11=determinant(a,4);
       
       // Find determinant M12
       for (int i=0;i<4;i++)
       {
           a[i][0]=p[i](0)*p[i](0)+p[i](1)*p[i](1)+p[i](2)*p[i](2);
           a[i][1]=p[i](1);
           a[i][2]=p[i](2);
           a[i][3]=1;
       }
       double m12=determinant(a,4);
       
       // Find determinant M13
       for (int i=0;i<4;i++)
       {
           a[i][0]=p[i](0);
           a[i][1]=p[i](0)*p[i](0)+p[i](1)*p[i](1)+p[i](2)*p[i](2);
           a[i][2]=p[i](2);
           a[i][3]=1;
       }
       double m13=determinant(a,4);
       
       // Find determinant M14
       for (int i=0;i<4;i++)
       {
           a[i][0]=p[i](0);
           a[i][1]=p[i](1);
           a[i][2]=p[i](0)*p[i](0)+p[i](1)*p[i](1)+p[i](2)*p[i](2);
           a[i][3]=1;
       }
       double m14=determinant(a,4);
       
       // Find determinant M15
       for (int i=0;i<4;i++)
       {
           a[i][0]=p[i](0)*p[i](0)+p[i](1)*p[i](1)+p[i](2)*p[i](2);
           a[i][1]=p[i](0);
           a[i][2]=p[i](1);
           a[i][3]=p[i](2);
       }
       double m15=determinant(a,4);
      // Calculate the centre and the radius
       // If M11=0 then points define one or infinitly many spheres - take r=0
       if (m11>=0.00000001||m11<=-0.00000001)
       {
           p[4](0)=0.5*m12/m11;
           p[4](1)=0.5*m13/m11;
           p[4](2)=0.5*m14/m11;
           r=sqrt(p[4](0)*p[4](0)+p[4](1)*p[4](1)+p[4](2)*p[4](2)-m15/m11);
       }
       else
       {
           p[4](0)=0.0;
           p[4](1)=0.0;
           p[4](2)=0.0;
       }
       
       // Free space used for minors
       for (int i=0;i<4;i++)
       {
           free(a[i]);
       }
       free(a);
       
       return r;
} 

// Compute of maximum eigenvalue of the Jacobian of the Flux
// through the boundary
double DGBoundaryCell::computeMaxEVal()
{
	if ( right)return left->theConservationLaw->maximumEigenValue(n,left->theMean,right->theMean);
	else return left->theConservationLaw->maximumEigenValue(n,left->theMean,left->theMean);
}

// Compute the Jacabian Matrix of the Flux through the boundary
void DGBoundaryCell::computeJacobian(int Side, double* Up0, double** Jac)
{
    mVector nn(n);
    if (Side==0) nn*=(-1.);
	left->theConservationLaw->computeJacobian(nn,Up0,Jac);
}


// Computes change in Flux along boundary cell.
void DGBoundaryCell::computeDeltaF(int Side, double* Up0, double* deltaUp0,double* flux)
{       
   mPoint point;
   mVector flux1[MaxNbEqn];
   mVector flux2[MaxNbEqn];
   double Q1[MaxNbEqn];
   
   for (int i=0; i<cSize; i++) Q1[i]=Up0[i]+deltaUp0[i];
   left->theConservationLaw->Fi(point,Q1,flux1);
   left->theConservationLaw->Fi(point,Up0,flux2);
  
   // Calculate the DELTA FLUX
    for (int i=0; i<5; i++)
    {
        flux[i]=(flux1[i](0)-flux2[i](0))*n(0)+(flux1[i](1)-flux2[i](1))*n(1)+(flux1[i](2)-flux2[i](2))*n(2);
        if (Side==0) flux[i]=flux[i]*(-1);
    }
}

void DGCell::ZeroError ()
{
  error = 0;
}

void DGCell::ZeroJump ()
{
  jump = 0;
  jumpError=0;
}

/*********************** U T I L S *****************************/

// compute  of variables and their derivatives on the cell
void DGCell::computeMean ()
{
  
 if(theFunctionSpace->isOrthogonal())
 {
	 double firstBasisFunction = pt(0).fcts[0];
	for(int j=0;j<cSize;++j)theMean[j] = theFieldsCoefficients->get(j,0)* firstBasisFunction;
	}
 else 
 {
	 int i,j; 
  for(i=0;i<cSize;++i)theMean[i] = 0.0;
  double vol = 0.0;
  double val[MaxNbEqn];
    for(i=0;i<nbPtGauss;++i)
    {
		volumeGaussPoint  &pg = pt(i);
      interpolate(pg.fcts, val);
      vol += pg.JacTimesWeight;
      for(j=0;j<cSize;++j) theMean[j] +=pg.JacTimesWeight*val[j];
    }
  for(j=0;j<cSize;++j)theMean[j] /= vol;
 }
 }

void DGCell::adaptTimeStep (double CFLMAX, double &DT)
{
  // CFL = a DT / DX
  // CFL < CFLMAX -> DT < DX CFLMAX / a
  double a = computeMAXEV();
  double deltat = cellSize*CFLMAX / (a * (2.*fOrder + 1.0));

  if (DT > deltat)DT = deltat;
}

double DGCell::computeMAXEV()
{
   switch(theMeshEntity->getType())
	{
	case 2:        /*Triangle*/
		{
  volumeGaussPoint &pg = pt(0);
  mPoint p(pg.x,pg.y,pg.z);  
  double a = theConservationLaw->maximumEigenValue(theMean,p);
  double val[MaxNbEqn];
  double u,v,w,b;
  u=0;v=0;w=0;
  interpolate(u,v,w,val);
  b  = theConservationLaw->maximumEigenValue(val,p);
  a = (b >a ? b : a);
  u=1.;v=0;w=0;
  interpolate(u,v,w,val);
  b  = theConservationLaw->maximumEigenValue(val,p);
  a = (b >a ? b : a);
   u=0;v=1.;w=0;
  interpolate(u,v,w,val);
  b  = theConservationLaw->maximumEigenValue(val,p);
  a = (b >a ? b : a);
  return a;
		}
		break;
	case 3:			/*Square*/
		{
			mPoint p((pMax(0)+pMin(0))*0.5,(pMax(1)+pMin(1))*0.5);
			return theConservationLaw->maxEigenValueCube(pMax(0)-pMin(0),pMax(1)-pMin(1),0.0,theMean,p);
		}
		break;
	case 4:
		{
			mPoint p((pMax(0)+pMin(0))*0.5,(pMax(1)+pMin(1))*0.5,(pMax(2)+pMin(2))*0.5);
		return theConservationLaw->maximumEigenValue(theMean,p);
		}
		break;
	default:
		printf("computeMAXEV: this type has not been coded yet\n");
		return 0.0;
   }
}

void DGCell::computeCellSize()
{	
  switch(theMeshEntity->getType())
    {
    case 2:        /*Triangle*/
      {
	mPoint p1 = ((Vertex*)theMeshEntity->get(0,0))->point();
	mPoint p2 = ((Vertex*)theMeshEntity->get(0,1))->point();
	mPoint p3 = ((Vertex*)theMeshEntity->get(0,2))->point();
	mVector A(p1,p2);
	mVector B(p2,p3);
	mVector C(p3,p1);
	double a = sqrt(A*A);
	double b = sqrt(B*B);
	double c = sqrt(C*C);
	double halfP = 0.5 * (a+b+c); 
	cellSize = sqrt(halfP*(halfP-a)*(halfP-b)*(halfP-c))/halfP*2.;
      }
      break;
    case 3:         /*Square*/
      {
	mPoint p1 = ((Vertex*)theMeshEntity->get(0,0))->point();
	mPoint p2 = ((Vertex*)theMeshEntity->get(0,1))->point();
	mVector a(p1,p2);
	cellSize = a.L2norm(); 
      }
    case 4:        /*Tet*/
      {
	mPoint p1 = ((Vertex*)theMeshEntity->get(0,0))->point();
	mPoint p2 = ((Vertex*)theMeshEntity->get(0,1))->point();
	mPoint p3 = ((Vertex*)theMeshEntity->get(0,2))->point();
	mPoint p4 = ((Vertex*)theMeshEntity->get(0,3))->point();
	mVector A(p1,p2);
	mVector B(p2,p3);
	mVector C(p3,p1);
	mVector D(p1,p4);
	mVector E(p2,p4);
	mVector F(p3,p4);
	cellSize = A.L2norm();
	if (cellSize > B.L2norm() )cellSize = B.L2norm();
	if (cellSize > C.L2norm() )cellSize = C.L2norm();
	if (cellSize > D.L2norm() )cellSize = D.L2norm();
	if (cellSize > E.L2norm() )cellSize = E.L2norm();
	if (cellSize > F.L2norm() )cellSize = F.L2norm();
	break;
      }
    default : printf("size function for this case is not coded yet\n");
      cellSize=0.0;
    }
}

void DGCell::computeVolume()
{	
  switch(theMeshEntity->getType())
    {
    case 2:        /*Triangle*/
      volume = 0.5*detJac;
      //volumeRatio = 1./pt(0).detJac;
      
      break;
    case 3:         /*Square*/
      {
	mPoint p1 = ((Vertex*)theMeshEntity->get(0,0))->point();
	mPoint p2 = ((Vertex*)theMeshEntity->get(0,1))->point();
	mPoint p3 = ((Vertex*)theMeshEntity->get(0,2))->point();
	mVector a(p1,p2);
	mVector b(p2,p3);
	volume = a.L2norm()*b.L2norm();
	//volumeRatio = 4./volume;
      }
      break;
    case 4:			/*Tetrahedron*/
      {
	volume = detJac/6.;
	// volumeRatio = 1./pt(0).detJac;
      }
      break;
    default : printf("volume function for this case is not coded yet\n");
    }
}

void DGCell::computePerimeter()
{	
  switch(theMeshEntity->getType())
    {
    case 2:        /*Triangle*/
      {
	mPoint p1 = ((Vertex*)theMeshEntity->get(0,0))->point();
	mPoint p2 = ((Vertex*)theMeshEntity->get(0,1))->point();
	mPoint p3 = ((Vertex*)theMeshEntity->get(0,2))->point();
	mVector a(p1,p2);
	mVector b(p2,p3);
	mVector c(p3,p1);
	perimeter = a.L2norm()+b.L2norm()+c.L2norm();
      }
      break;
    case 3:        /*Quad*/       
      {
	mPoint p1 = ((Vertex*)theMeshEntity->get(0,0))->point();
	mPoint p2 = ((Vertex*)theMeshEntity->get(0,1))->point();
	mPoint p3 = ((Vertex*)theMeshEntity->get(0,2))->point();
	mPoint p4 = ((Vertex*)theMeshEntity->get(0,3))->point();
	mVector a(p1,p2);
	mVector b(p2,p3);
	mVector c(p3,p4);
	mVector d(p4,p1);
	perimeter = a.L2norm()+b.L2norm()+c.L2norm()+d.L2norm();
	  }
      break;
    default : //printf("perimeter function for this case is not coded yet\n");
		perimeter=0.0;
    }
}

void DGCell::L2Proj (FieldEvaluator *f, double *proj)
{
  double val[MaxNbEqn];
  int i, j, k;
  const int s = cSize * fSize;

  for( i=0;i<s;++i) proj[i] = 0.0;
  
  for(i=0;i<nbPtGauss;++i)
    {
      volumeGaussPoint &pg=pt(i);
      f->eval(mPoint(pg.x,pg.y,pg.z),0.0,val);
     
      for( j=0;j<cSize;++j)
	for( k=0;k<fSize;++k)
	  proj[k+fSize*j] += val[j] * pg.fcts[k] * pg.JacTimesWeight;
    }
}

void DGCell::L2Proj (double* dU)
{
  double val[MaxNbEqn];
  int i, j, k;
  const int dof = cSize * fSize;

  for( i=0;i<dof;++i) theRightHandSide[i] = 0.0;
  
  for(i=0;i<nbPtGauss;++i)
    {
      volumeGaussPoint &pg=pt(i);
      interpolate(pg.fcts,val);
     for( j=0;j<cSize;++j)	val[j] += dU[j];
      for( j=0;j<cSize;++j)
	for( k=0;k<fSize;++k)
	  theRightHandSide[k+fSize*j] += val[j] * pg.fcts[k] * pg.JacTimesWeight;
    }

  if (theFunctionSpace->isOrthogonal())
    {
	 double inv_Jac = 1./detJac;
	 double* coeff = theFieldsCoefficients->get();
      for(i=0;i<dof;++i) coeff[i] = theRightHandSide[i]*inv_Jac;
    }
  else printf("Multigrid projection is not done for not orthogonal spaces\n");
}


void DGCell::L2ProjInitial (FieldEvaluator *f)
{
  int i;
  L2Proj(f,theRightHandSide);		 
  if (theFunctionSpace->isOrthogonal())
    {
	 double inv_Jac = 1./detJac;
	 int dof = cSize*fSize;
	 double* coeff = theFieldsCoefficients->get();
      for(i=0;i<dof;++i)
	  coeff[i] = theRightHandSide[i]*inv_Jac;
    }
  else
    {

		for(int k = 0;k<cSize;++k)		
	for(int i = 0;i<fSize;++i)
	  {
	    double dxi = 0.0;
		for(int j = 0;j<fSize;++j)
	      dxi += theInvertMassMatrix[i][j] * theRightHandSide[j+fSize*k];
	    
		theFieldsCoefficients->get(k,i) = dxi;
	  }
    }
}




void DGCell::getBox(double &XminBox,double &YminBox,double &ZminBox, 
		    double &XMaxBox,double &YMaxBox,double &ZMaxBox)
{
  XminBox = pMin(0);
  YminBox = pMin(1);
  ZminBox = pMin(2);
  XMaxBox = pMax(0);
  YMaxBox = pMax(1);
  ZMaxBox = pMax(2);
}

bool DGCell::inCell (const mPoint &p, double *val)
{
  /*
  printf("%f %f %f %f %f %f %f %f %f\n",p(0),p(1),p(2),pMin(0),pMin(1),pMin(2),
	 pMax(0),pMax(1),pMax(2));
  */
  if(p(0) > pMax(0) ||
     p(1) > pMax(1) ||
     p(2) > pMax(2) ||
     p(0) < pMin(0) ||
     p(1) < pMin(1) ||
     p(2) < pMin(2))return false;

  //  printf("in box...\n");
  
  double u,v,w;
  theMapping->invert(p(0),p(1),p(2),u,v,w);
  //  printf("%f %f %f\n",u,v,w);
  if(!theMapping->inReferenceElement(u,v,w))return false;
  //  printf("ok\n");
  interpolate(u,v,w,val);
  return true;
}

void DGCell::write(ostream &o)
{
  o << theFunctionSpace->order() << " ";
  o.precision(16);
  for(int i=0;i<cSize;++i)
    for(int j=0;j<fSize;++j)
      o << theFieldsCoefficients->get(i,j) << " ";
  
}

void DGCell::read(istream &is)
{
  for(int i=0;i<cSize;++i)
    {
      for(int j=0;j<fSize;++j) 
        is >> theFieldsCoefficients->get(i,j);
	}  
}

void DGCell::adaptOrder(int newOrder)
{
  int i,j;
  delete theFunctionSpace;
  switch(theMeshEntity->getType())
    {
    case mEntity::TRI  : theFunctionSpace = new OrthogonalTriangleFunctionSpace(newOrder); break;
    case mEntity::QUAD : theFunctionSpace = new QuadFunctionSpace(newOrder); break;
    case mEntity::TET  : theFunctionSpace = new OrthogonalTetFunctionSpace(newOrder); break;
    case mEntity::HEX  : theFunctionSpace = new HexFunctionSpace(newOrder); break;
    }
  int oldFSize = fSize;
  fSize                 = theFunctionSpace->size();
  fOrder                = theFunctionSpace->order();
  order = 2 * fOrder + theMapping->order()-1;

  DG_Dofs *oldFieldsCoefficients;
  oldFieldsCoefficients = new DG_Dofs(2,cSize,oldFSize);
  oldFieldsCoefficients->operator =(theFieldsCoefficients);
  delete theFieldsCoefficients;
  theFieldsCoefficients = new DG_Dofs(2,cSize,fSize);
  if (fSize>oldFSize)
    {
      for(i=0;i<cSize;++i)
	{
    for(j=0;j<oldFSize;++j)
	{
      theFieldsCoefficients->get(i,j) = oldFieldsCoefficients->get(i,j);
	  //printf("%e \n", oldFieldsCoefficients->get(i,j));
	}
    for(j=oldFSize;j<fSize;++j) 
      theFieldsCoefficients->get(i,j) = 0.;
	}
    }
else theFieldsCoefficients->operator =(oldFieldsCoefficients);
  delete [] theRightHandSide;
  theRightHandSide = new double [fSize*cSize];
  ZeroRHS();
	volumeGaussPoints.clear();
  freeMatrix(theInvertMassMatrix);
  theInvertMassMatrix = 0;
  init();
  int Sizenm1 = theMeshEntity->size(theMeshEntity->getLevel()-1);
  for(j=0;j<Sizenm1;j++)
    {
      mEntity *b = theMeshEntity->get(1,j);
      DGBoundaryCell *bc = (DGBoundaryCell*)b->getCell();
      bc->init();
    }
}


/*-------- Coeffs Management ------------*/

DG_Dofs::DG_Dofs (short n, short NbUnknowns, short FunctionSpaceSize)
  : Nb (n), NU (NbUnknowns), FS(FunctionSpaceSize)
{
  theFieldsCoefficients = new double* [Nb];
  for(int i=0;i<Nb;++i) theFieldsCoefficients[i] = new double [NbUnknowns * FunctionSpaceSize];
}

DG_Dofs* DG_Dofs::operator= (const DG_Dofs *other)
{
  Nb=other->Nb;
  NU=other->NU;
  FS=other->FS;
  //theFieldsCoefficients = new double *[Nb];
  //for(int i=0;i<Nb;++i) theFieldsCoefficients[i] = new double [NU * FS];
  for(int i=0;i<NU;++i) 
	  for (int j=0;j<FS;++j)
	  {
		  theFieldsCoefficients[0][i*FS+j] = other->theFieldsCoefficients[0][i*FS+j];
		  //printf(" %d %e %e \n",FS,theFieldsCoefficients[0][i*FS+j],other->theFieldsCoefficients[0][i*FS+j]);
	  }
  return this;
}



DG_Dofs::~DG_Dofs ()
{
  for(int i=0;i<Nb;++i) delete [] theFieldsCoefficients[i];
  delete [] theFieldsCoefficients;
}

void DG_Dofs::copy()   //copies solution coefficients kept in the fist vector to the second vector in DG_Dofs
{
const int N = NU*FS;
for (int i=0;i<N;++i)
theFieldsCoefficients[1][i] = theFieldsCoefficients[0][i];  
}

void DG_Dofs::swap()   //copies solution coefficients kept in the fist vector to the second vector in DG_Dofs
{
  const int N = NU*FS;
  double tmp;
  for (int i=0;i<N;++i)
    {
      tmp = theFieldsCoefficients[0][i] ; 
      theFieldsCoefficients[0][i] = theFieldsCoefficients[1][i];
      theFieldsCoefficients[1][i] = tmp; 
    }
}

void DG_Dofs::copyBack() //copies solution coeff from the second vector to the first.
{
const int N = NU*FS;
for (int i=0;i<N;++i)
theFieldsCoefficients[0][i] = theFieldsCoefficients[1][i];  
}

void DG_Dofs::copyLimitedCoeffs(int k) //copies solution coeff from the second vector to the first after limiting
{
const int N = NU*FS;
 for (int q=0;q<NU;q++)
   for (int i=q*FS+k;i<FS*(q+1);++i)
     theFieldsCoefficients[0][i] = theFieldsCoefficients[1][i];  
}

void DG_Dofs:: eval (vector<double> &ff, double *field, double t) 
{
  //  double *b = 0;
  double *a,A,F;
 
  //double frac = (t-t0)/(t1-t0);
  //  printf("%f %f %f %f %d\n",t,t0,t1,frac,frac != 0.0);
  
  for(int i=0;i<NU;++i)
    {
      a = & theFieldsCoefficients[0][i*FS];
      //  if(frac != 0.0) 
      //	{
      //  b = &theFieldsCoefficients[1][i*FS];
//}
      field[i] = 0.0; 
      for(int j=0;j<FS;++j)
	{
	  //  if(b)field[i] += ff[j] * ((1.-frac) * a[j] + (frac)*b[j]);
	  // else 
	  F = ff[j]; A = a[j];
	  //field[i] += ff[j] * a[j];
	  field[i] += F * A;
	 	}
    }      
}

void DG_Dofs:: eval (vector<double> &ff, double *field)
{
  double *a, A,F;
  //  double *k = &ff.front();
  for(int i=0;i<NU;++i)
    {
      a = & theFieldsCoefficients[0][i*FS];
      field[i] = 0.0;
      for(int j=0;j<FS;++j)
	{
	  F=ff[j]; A=a[j];
	  //field[i] += ff[j] * a[j];
	  field[i] += F * A;
	  //	  printf("%f ",a[j]);
	}
      //      printf("\n");
    }      
}

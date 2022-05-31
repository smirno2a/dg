// Vertex.cpp: implementation of the Vertex class.
//
//////////////////////////////////////////////////////////////////////
#include <stdio.h>
#include "Vertex.h"
#include "mException.h"
#include "mVector.h"
#include "mIdGenerator.h"
#include "mEntity.h"
#include <stdlib.h>
#include <math.h>
Vertex::Vertex(mIdGenerator &theIdGenerator, const mPoint &pt, gEntity *classif)
  :mEntity(),p(pt)
{
  iD = theIdGenerator.generateId();
  theClassification = classif;
}


Vertex::Vertex(int theId, const mPoint &pt, gEntity *classif)
  :mEntity(),p(pt)
{
  iD = theId;
  theClassification = classif;
}

unsigned long Vertex::getRAND()
{
  srand(iD);
  unsigned long RAND = rand();
  return RAND; 
}

void Vertex::deleteId(mIdGenerator &theIdGenerator)
{
  theIdGenerator.addId(iD);
  iD = 0;
}

void Vertex :: print() const
{
  printf("vertex Id %d = (%12.5E,%12.5E,%12.5E)\n",getId(),p(0),p(1),p(2));
}



double Vertex:: distance(Vertex* v2){

double		dist=pow(p(0)-v2->p(0),2)+pow(p(1)-v2->p(1),2);
dist = sqrt(dist);
		return dist;
}

void Vertex::moveposition(mPoint* np){
    p(0)=(*np)(0);
		p(1)=(*np)(1);
		p(2)=(*np)(2);
}
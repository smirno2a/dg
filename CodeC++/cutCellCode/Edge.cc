// Edge.cpp: implementation of the Edge class.
#include <stdio.h>
#include <math.h>
#include "Edge.h"
#include "mEntity.h"
#include "Vertex.h"

Edge::Edge(Vertex *v1, Vertex *v2, gEntity *classification)
{
  theAdjacencies[0] = new mDownwardAdjacencyContainer (2);
  theAdjacencies[0]->add(v1);
  theAdjacencies[0]->add(v2);
  theClassification = classification;
  unsigned long int i1 = v1->getId();
  unsigned long int i2 = v2->getId();
  iD = v1->getRAND() + v2->getRAND();
}

Edge::~Edge()
{
//printf("Deleting edge with id = %d \n",iD);
int j, nb;
if (theAdjacencies[0])
{
nb = theAdjacencies[0]->size();
for (j=0;j<nb;++j)                   //deleting pointers to this element from its vertices
theAdjacencies[0]->get(j)->del(this);
}
if (theAdjacencies[1])
{
nb = theAdjacencies[1]->size();     //deleting pointer to this element from its parent(?) or child(?)
for (j=0;j<nb;++j)
theAdjacencies[1]->get(j)->del(this);
}
if (theAdjacencies[2])
{
printf("faces should be destructed first \n");
exit(0);
}
if (theAdjacencies[3])
{
printf("~Edge: don't know how to destruct in3D\n");
exit(0);
}
}

int Edge::getNbTemplates(int what)const
{
  if(what == 0) return 2;
  else return 0;
}

Vertex *Edge::vertex(int i) const
{
  return (Vertex*)theAdjacencies[0]->get(i);
}

Vertex *Edge::commonVertex(Edge *other) const
{
  if(other->vertex(0) == vertex(0))return vertex(0);
  if(other->vertex(1) == vertex(0))return vertex(0);
  if(other->vertex(0) == vertex(1))return vertex(1);
  if(other->vertex(1) == vertex(1))return vertex(1);
  return 0;
}

void Edge::print()const
{
  printf("edge %p with Id %d and vertices %d %d\n",this,iD,get(0,0)->getId(),
	 get(0,1)->getId());
}

// Face.cpp: implementation of the Face class.
#include "Face.h"
#include "mEntity.h"
#include "Vertex.h"
#include "Edge.h"
#include <stdio.h>

Face::Face(mDownwardAdjacencyContainer *lis ,gEntity *classification)
  : mEntity(lis,classification,2)
{
  theType = TRI;
}


Face::Face(Vertex *v1, Vertex *v2, Vertex *v3,gEntity *classification)
{
  theType = TRI;
  theAdjacencies[0] = new mDownwardAdjacencyContainer (3);
  theAdjacencies[0]->add(v1);
  theAdjacencies[0]->add(v2);	
  theAdjacencies[0]->add(v3);
  theClassification = classification;
  iD = v1->getRAND()+ v2->getRAND() + v3->getRAND();
}

Face::Face(Edge *e1, Edge *e2, Edge *e3,gEntity *classification)
{
  theType = TRI;
  theAdjacencies[1] = new mDownwardAdjacencyContainer (3);
  theAdjacencies[1]->add(e1);
  theAdjacencies[1]->add(e2);	
  theAdjacencies[1]->add(e3);	
  theClassification = classification;
  /*this will give the double of the result*/
  iD = e1->getId() + e2->getId() + e3->getId(); 
  iD /=2;
}

Face::Face(Vertex *v1, Vertex *v2, Vertex *v3, Vertex *v4, gEntity *classification)
{
  theType = QUAD;
  theAdjacencies[0] = new mDownwardAdjacencyContainer (4);
  theAdjacencies[0]->add(v1);
  theAdjacencies[0]->add(v2);	
  theAdjacencies[0]->add(v3);	
  theAdjacencies[0]->add(v4);	
  theClassification = classification;
  iD = v1->getRAND()+ v2->getRAND() + v3->getRAND() + v4->getRAND();
}

Face::Face(Edge *e1, Edge *e2, Edge *e3, Edge *e4, gEntity *classification)
{
  theType = QUAD;
  theAdjacencies[1] = new mDownwardAdjacencyContainer (4);
  theAdjacencies[1]->add(e1);
  theAdjacencies[1]->add(e2);	
  theAdjacencies[1]->add(e3);	
  theAdjacencies[1]->add(e4);	
  theClassification = classification;
  /*this will give the double of the result*/
  iD = e1->getId() + e2->getId() + e3->getId() + e4->getId(); 
  iD /=2;
}

Face::Face(Vertex *v1, Vertex *v2, Vertex *v3, Vertex *v4, Vertex *v5, gEntity *classification)
{
  theType = PANT;
  theAdjacencies[0] = new mDownwardAdjacencyContainer (5);
  theAdjacencies[0]->add(v1);
  theAdjacencies[0]->add(v2);	
  theAdjacencies[0]->add(v3);	
  theAdjacencies[0]->add(v4);
  theAdjacencies[0]->add(v5);
  theClassification = classification;
  iD = v1->getRAND()+ v2->getRAND() + v3->getRAND() + v4->getRAND() + v5->getRAND();
}

Face::Face(Edge *e1, Edge *e2, Edge *e3, Edge *e4, Edge *e5, gEntity *classification)
{
  theType = PANT;
  theAdjacencies[1] = new mDownwardAdjacencyContainer (5);
  theAdjacencies[1]->add(e1);
  theAdjacencies[1]->add(e2);	
  theAdjacencies[1]->add(e3);	
  theAdjacencies[1]->add(e4);
  theAdjacencies[1]->add(e5);
  theClassification = classification;
  /*this will give the double of the result*/
  iD = e1->getId() + e2->getId() + e3->getId() + e4->getId() + e5->getId(); 
  iD /=2;
}
int Face:: getNbTemplates (int what) const
{
	switch(getType()){
			case TRI:
				return 3;
			case QUAD:
				return 4;
			case PANT:
				return 5;
			default:
				return 0;
	}
  //return(getType() == TRI)?3:4;
}

mEntity *Face::getTemplate(int ith, int what, int with)const
{
  switch(what)
    {
    case 0:
      if(theAdjacencies[0])return (*theAdjacencies[0])[ith];
      if(theAdjacencies[1])
	{
	  Edge *e1 = (Edge*)get(1,ith);
	  Edge *e2 = (Edge*)get(1,(ith+1)%getNbTemplates(1));
	  return e1->commonVertex(e2);
	}
      break;
    case 1:
      return new Edge ((Vertex*)theAdjacencies[0]->get(ith),
			(Vertex*)theAdjacencies[0]->get((ith+1)%(theAdjacencies[0]->size())),
			theClassification);
      break;
    }
  return (mEntity*) 0;
}

Vertex * Face::commonVertex (Face *f1, Face *f2)
{
  set<Vertex*,EntityLessThanKey> thisVertices;
  set<Vertex*,EntityLessThanKey> f1Vertices;
  set<Vertex*,EntityLessThanKey> f2Vertices;
  getVertices(thisVertices);
  f1->getVertices(f1Vertices);
  f2->getVertices(f2Vertices);
  for( set<Vertex*,EntityLessThanKey>::const_iterator iter = thisVertices.begin();
       iter != thisVertices.end();
       ++iter)
    {
      if(f1Vertices.find(*iter) != f1Vertices.end() &&
	 f2Vertices.find(*iter) != f2Vertices.end())return *iter;
    }
  return 0;
}

Edge * Face::commonEdge (Face *f1)
{
  if(!isAdjacencyCreated(1))return 0;

  for(int i=0;i<size(1);i++)
    {
      mEntity *e1 = get(1,i);
      for(int j=0;j<f1->size(1);j++)
	{
	  if(e1 == f1->get(1,j))return (Edge*)e1;
	}      
    }
  return 0;
}

Face::~Face()
{
//  if(theClassification)
    for(int i=0;i<4;i++)
      if(i <= getLevel() && theAdjacencies[i])
	for(int j=0;j<theAdjacencies[i]->size();j++)
	  theAdjacencies[i]->get(j)->del(this);
}


void Face::print()const
{ 
	if(size(0) == 5)
    printf("face %d with vertices %d %d %d %d %d\n",iD,get(0,0)->getId(),
	   get(0,1)->getId(),get(0,2)->getId(),get(0,3)->getId(),get(0,4)->getId());
  else if(size(0) == 4)
    printf("face %d with vertices %d %d %d %d\n",iD,get(0,0)->getId(),
	   get(0,1)->getId(),get(0,2)->getId(),get(0,3)->getId());
  else
    printf("face %d with vertices %d %d %d\n",iD,get(0,0)->getId(),
	   get(0,1)->getId(),get(0,2)->getId());    
}


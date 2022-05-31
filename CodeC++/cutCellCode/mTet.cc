// mTet.cpp: implementation of the mTet class.
#include "mTet.h"
#include "mEntity.h"
#include "mException.h"
#include "Vertex.h"
#include "Edge.h"
#include "Face.h"


mTet::mTet(Vertex *v1, Vertex *v2, Vertex *v3, Vertex *v4, gEntity *classification)
{
  for(int i=1;i<4;i++)theAdjacencies[i] = (mAdjacencyContainer*)0;
  theAdjacencies[0] = new mDownwardAdjacencyContainer (4);
  theAdjacencies[0]->add(v1);
  theAdjacencies[0]->add(v2);	
  theAdjacencies[0]->add(v3);	
  theAdjacencies[0]->add(v4);	
  theClassification = classification;
  iD = v1->getRAND() + v2->getRAND() + v3->getRAND() + v4->getRAND();
  //iD = v1->getId() + v2->getId() + v3->getId() + v4->getId();
  theAttachable=0;
}

mTet::mTet(Face *f1, Face *f2, Face *f3, Face *f4, gEntity *classification)
{
  for(int i=1;i<4;i++)theAdjacencies[i] = (mAdjacencyContainer*)0;
  theAdjacencies[2] = new mDownwardAdjacencyContainer (4);
  theAdjacencies[2]->add(f1);
  theAdjacencies[2]->add(f2);	
  theAdjacencies[2]->add(f3);	
  theAdjacencies[2]->add(f4);	
  theClassification = classification;
  computeId();
  theAttachable=0;
}

int mTet:: getNbTemplates (int what) const
{
  switch(what)
    {
    case 0: return 4; break;
    case 1: return 6; break;
    case 2: return 4; break;
    default : throw new mException (__LINE__,__FILE__,"");
    }
}

static int Tev[6][2] = {{0,1},{0,2},{0,3},{1,2},
			{1,3},{2,3}};
static int Tfv[4][3] =   {{0,1,2},
			  {0,1,3},
			  {0,2,3},
			  {1,2,3}};
static int Tve[4][2] = {{0,1},{0,3},{1,3},{2,4}};
static int Tvf[4][2] = {{0,1},{0,3},{1,3},{2,4}};

mEntity *mTet::getTemplate(int ith, int what, int with)const
{
  switch(what)
    {
    case 0:
      if(theAdjacencies[1])
	return (((Edge*)theAdjacencies[1]->get(Tve[ith][0]))->commonVertex((Edge*)theAdjacencies[1]
									    ->get(Tve[ith][1])));
      else if (theAdjacencies[2])
	{
	  return 0;
	}
      break;
    case 1:
      return new Edge ((Vertex*)theAdjacencies[0]->get(Tev[ith][0]),
			(Vertex*)theAdjacencies[0]->get(Tev[ith][1]),
			theClassification);
      break;
    case 2:
      if(with == 0)
	return new Face ((Vertex*)theAdjacencies[0]->get(Tfv[ith][0]),
			  (Vertex*)theAdjacencies[0]->get(Tfv[ith][1]),
			  (Vertex*)theAdjacencies[0]->get(Tfv[ith][2]),
			  theClassification);
      else if (with == 1) throw new mException (__LINE__,__FILE__,"not done yet");
    default : throw new mException (__LINE__,__FILE__,"this tet template is not done");
    }
  return 0;
}



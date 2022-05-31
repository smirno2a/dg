// Edge.h: interface for the Edge class.
//
//////////////////////////////////////////////////////////////////////

#ifndef _EDGE_H_
#define _EDGE_H_
#include "mEntity.h"
class gEntity;
class Vertex;

class Edge : public mEntity 
{
public:
  Edge (Vertex *v1, Vertex *v2,gEntity *classification);
  ~Edge();
  Vertex* vertex(int) const;
  Vertex* commonVertex (Edge *) const;
  virtual int getLevel() const {return 1;} 
  virtual mType getType() const{return EDGE;} 
  virtual void print() const;
  virtual int getNbTemplates (int what) const;
};

#endif 

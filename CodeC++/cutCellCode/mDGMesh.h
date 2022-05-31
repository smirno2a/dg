#ifndef _DGMESH_H_
#define _DGMESH_H_

#include <list>
#include "mMesh.h"
#include "mPoint.h"

class mEntity;
using namespace std;
 
class mDGMesh : public mMesh
{
  Vertex* splitEdge (Edge *);
  void splitTriangle (Face *);
  void unsplitTriangle(Face *e);
  void splitQuad (Face *);
  void reconnect (int dim);
  protected:
   mMeshEntityContainer allSplittedEntities;
  public :
    iter beginsplit(int what){return allSplittedEntities.begin(what);} 
    iter endsplit(int what){return allSplittedEntities.end(what);}
    void split  (int dim, list<mEntity *> &toSplit);
	void split (mEntity* );
    void unsplit(int dim, list<mEntity *> &toUnSplit);
	void unsplit (mEntity* );
    void getLeaves(mEntity *, list<mEntity*> &leaves);
    int getRefinementDepth(mEntity *);
    int getRefinementLevel(mEntity *);

		void adjustcc(int dim, double percent);
		Face* findnb(Face* self, Edge *e);
		Edge* findbd(Face* self, int ID);
		Vertex* findend(Edge* e, Vertex* p);
		Vertex* leftend(Face* self, Edge* e);
		Vertex* rightend(Face* self, Edge* e);
		mPoint intersection(mPoint* A, Vertex* B);
		void splitPant(Face* self, Vertex* A, Vertex* B);
		void merge(int dim);
		double longest(Face* self);
		void combineTQ(Face* q1, Face* t1, Edge* e1, Face* q2, Face* t2, Edge* e2);
		void resplitPant(int dim);
		void resplitQuad(int dim);
		void splitQuaddial(Face* self, Vertex* A);
		void createConnections2(int to);
        int lineinside(mPoint a, mPoint b, mPoint c);
		void dealQuad(int dim);
    
		void attachRefinementLevel( mEntity *e, int ref);
};

#endif

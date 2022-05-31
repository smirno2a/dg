#include "mImportExport.h"
#include "mDGMesh.h"
#include "mEntity.h"
#include "Vertex.h"
#include "Face.h"
#include "Edge.h"
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <map>
#include <vector>+
#define EPS 0.5
#define SIZEN 10000


Face* mDGMesh::findnb(Face* self, Edge *e){
	//find the neigbor face for edge e other than face self
	if(e->getClassification()->getId()==40000||e->getClassification()->getId()==20000){
		printf("boundary edge has only one neighbor\n");
	  return NULL;
	}
	else{
		Face* fc=(Face*)e->get(2,0);
    if  (fc!=self )
        return fc;
		else return (Face*)e->get(2,1);}
}

Edge* mDGMesh::findbd(Face* self, int ID){
	//return the boundary edge for face self
	int ii=0;
	Edge* e;
	for(ii=0;ii!=5;ii++) {
	  e=(Edge*)self->get(1,ii);
		if(e->getClassification()->getId()==ID){
			return e;
			break;
		}
	}
	return NULL;
}

Vertex* mDGMesh::findend(Edge* e, Vertex* p){
	//return the othen end of an edge
	if ((Vertex*)e->get(0,0)!=p)
     return (Vertex*)e->get(0,0);
	else return (Vertex*)e->get(0,1);
}

double mDGMesh::longest(Face* self){
    switch(self->getType()){
		case mEntity::PANT:{
			    //Edge* ee[5];
					Vertex* w[5];
					for(int i=0;i!=5;i++){
					    //ee[i]=(Edge*)self->get(1,i);
							w[i]=(Vertex*)self->get(0,i);
					}
					double d[5];
					double maxedge=0;
					for(int i=0;i!=5;i++){
						d[i]=w[i]->distance(w[(i+1)%5]);
						if (d[i]>maxedge)
							maxedge=d[i];
					}
					return maxedge;
											 }break;

			case mEntity::QUAD:{
			    //Edge* ee[5];
					Vertex* w[4];
					for(int i=0;i!=4;i++){
					    //ee[i]=(Edge*)self->get(1,i);
							w[i]=(Vertex*)self->get(0,i);
					}
					double d[4];
					double maxedge=0;
					for(int i=0;i!=4;i++){
						d[i]=w[i]->distance(w[(i+1)%4]);
						if (d[i]>maxedge)
							maxedge=d[i];
					}return maxedge;
											 }break;

				case mEntity::TRI:{
			    //Edge* ee[5];
					Vertex* w[3];
					for(int i=0;i!=3;i++){
					    //ee[i]=(Edge*)self->get(1,i);
							w[i]=(Vertex*)self->get(0,i);
					}
					double d[3];
					double maxedge=0;
					for(int i=0;i!=3;i++){
						d[i]=w[i]->distance(w[(i+1)%3]);
						if (d[i]>maxedge)
							maxedge=d[i];
					}return maxedge;
											 }break;
				default:return 0;
		}
}

Vertex* mDGMesh::leftend(Face* self, Edge* e){
	switch(self->getType()){
		case mEntity::PANT:{
			    Edge* ee[5];
					Vertex* w[5];
					for(int i=0;i!=5;i++){
					    ee[i]=(Edge*)self->get(1,i);
							w[i]=(Vertex*)self->get(0,i);
					}
					for(int i=0;i!=5;i++){
						if(e==ee[i]){
					        return w[i];
									break;
						}
					}
		}
											 break;
		case mEntity::QUAD:{
					Edge* ee[4];
					Vertex* w[4];
					for(int i=0;i!=4;i++){
					    ee[i]=(Edge*)self->get(1,i);
							w[i]=(Vertex*)self->get(0,i);
					}
					for(int i=0;i!=4;i++){
						if(e==ee[i]){
					        return w[i];
									break;
						}
					}
											 }
											 break;
		case mEntity::TRI:{
					Edge* ee[3];
					Vertex* w[3];
					for(int i=0;i!=3;i++){
					    ee[i]=(Edge*)self->get(1,i);
							w[i]=(Vertex*)self->get(0,i);
					}
					for(int i=0;i!=3;i++){
						if(e==ee[i]){
					        return w[i];
									break;
						}
					}
											}
											break;
		default:return NULL;
	}
}

Vertex* mDGMesh::rightend(Face* self, Edge* e)
{
	Vertex* left=leftend(self,e);
	return findend(e, left);
}

void mDGMesh::adjustcc(int dim, double percent){
    iter it;
		iter end_m1 = end (dim-1);
		for (it= begin(dim-1);it !=end_m1;++it)
		{
			mEntity *e=*it;
			if(e->getClassification()->getId()==40000||e->getClassification()->getId()==20000)
			{ 
				Vertex* v1=(Vertex*)e->get(0,0);
				Vertex* v2=(Vertex*)e->get(0,1);
				double length=v1->distance(v2);
				Face *f0=(Face*)e->get(2,0);
				switch(f0->getType()){
				case mEntity::TRI :
					{
					}
					break;
				case mEntity::PANT :{
					Edge* e[5];
					Vertex* w[5];
					for(int i=0;i!=5;i++){
					    e[i]=(Edge*)f0->get(1,i);
							w[i]=(Vertex*)f0->get(0,i);
					}
					double d[5];
					double maxedge=0;
					for(int i=0;i!=5;i++){
						d[i]=w[i]->distance(w[(i+1)%5]);
						if (d[i]>maxedge)
							maxedge=d[i];
					}
					for(int i=0;i!=5;i++){
						if(e[i]->getClassification()->getId()==40000||e[i]->getClassification()->getId()==20000){
							double ratio=d[i]/maxedge;
							if(ratio<EPS){
							  Face* nbleft=findnb(f0,e[(i+4)%5]);							
						    Face* nbright=findnb(f0,e[(i+1)%5]);
								Edge* leftbd;
								Edge* rightbd;
								if(e[i]->getClassification()->getId()==20000){
							      leftbd=findbd(nbleft, 20000);
								    rightbd=findbd(nbright, 20000);
								}
								else{
							      leftbd=findbd(nbleft, 40000);
								    rightbd=findbd(nbright, 40000);
								}
							  Vertex* left=leftend(nbleft,leftbd);
							  Vertex* up=leftend(f0,e[(i+4)%5]);
							  mPoint pp=(w[i]->point()*((1+ratio)/2.0)+left->point()*((1-ratio)/2.0));
							//	Vertex* temp=createVertex(pp(0),pp(1),pp(2),leftbd->getClassification());
               // temp->print();
								mPoint newleft=intersection(&pp, up);
								//allEntities.del(temp);
								//  mPoint newposit=newleft->point();
									w[i]->moveposition(&newleft);
								//	w[i]->print();

                Vertex* right=rightend(nbright,rightbd);
								Vertex* down=rightend(f0,e[(i+1)%5]);
								pp=(w[(i+1)%5]->point()*((1+ratio)/2.0)+right->point()*((1-ratio)/2.0));
								//temp=createVertex(pp(0),pp(1),pp(2),leftbd->getClassification());
								mPoint newright=intersection(&pp, down);
								//allEntities.del(temp);
								//newposit=newright->point();
								w[(i+1)%5]->moveposition(&newright);
                 
							}
						}	
					}
					
														}break;
				default:break;

				}
				//v1->print();
			}

		}
}

void mDGMesh::splitQuaddial(Face* self, Vertex* A){
    Edge* e[4];
		Vertex* w[4];
		int a;
		for(int i=0;i!=4;i++){
				e[i]=(Edge*)self->get(1,i);
				w[i]=(Vertex*)self->get(0,i);
				if (w[i]==A) a=i;
		}
		Edge* g11=e[a];
		Edge* g12=e[(a+1)%4];
		Edge* g21=e[(a+2)%4];
		Edge* g22=e[(a+3)%4];

		Edge* enew=new Edge(A,w[(a+2)%4],self->getClassification());

		Face* f1=new Face(w[a],w[(a+1)%4],w[(a+2)%4],self->getClassification());
		Face* f2=new Face(w[(a+2)%4],w[(a+3)%4],w[a],self->getClassification());

		f1->add(g11); f1->add(g12); f1->add(enew);
		f2->add(g21); f2->add(g22); f2->add(enew);

		enew->add(f1); enew->add(f2);

		g11->del(self); g11->add(f1);
		g12->del(self); g12->add(f1);
		g21->del(self); g21->add(f2);
		g22->del(self); g22->add(f2);

		allEntities.add(enew);
		allEntities.add(f1);
		allEntities.add(f2);

		allEntities.del(self);
}

void mDGMesh::splitPant(Face* self, Vertex *A, Vertex *B)
{//using AB to split face self
    Edge* e[5];
		Vertex* w[5];
		int a, b;
		for(int i=0;i!=5;i++){
				e[i]=(Edge*)self->get(1,i);
				w[i]=(Vertex*)self->get(0,i);
				if (w[i]==A) a=i;
				if (w[i]==B) b=i;
		}
		if ((a+2)%5==b){
		    Edge* g11=e[a];
				Edge* g12=e[(a+1)%5];
				Edge* g21=e[b];
				Edge* g22=e[(b+1)%5];
				Edge* g23=e[(b+2)%5];
				Edge* edgenew=new Edge(A, B, self->getClassification());

				Face *f1=new Face (w[a],w[(a+1)%5],w[b],self->getClassification());
				Face *f2=new Face (w[b],w[(b+1)%5],w[(b+2)%5],w[a],self->getClassification());
		    
				f1->add(g11); f1->add(g12); f1->add(edgenew);
				f2->add(g21); f2->add(g22); f2->add(g23); f2->add(edgenew);

				edgenew->add(f1); edgenew->add(f2);
				g11->del(self); g11->add(f1);
				g12->del(self); g12->add(f1);
				g21->del(self); g21->add(f2);
				g22->del(self); g22->add(f2);
				g23->del(self); g23->add(f2);

				allEntities.add(edgenew);
				allEntities.add(f1);
				allEntities.add(f2);

				allEntities.del(self);

		}

		else if ((a+3)%5==b){
		    Edge* g11=e[b];
				Edge* g12=e[(b+1)%5];
				Edge* g21=e[a];
				Edge* g22=e[(a+1)%5];
				Edge* g23=e[(a+2)%5];
				Edge* edgenew=new Edge(A, B, self->getClassification());

				Face *f1=new Face (w[b],w[(b+1)%5],w[a],self->getClassification());
				Face *f2=new Face (w[a],w[(a+1)%5],w[(a+2)%5],w[b],self->getClassification());
		    
				f1->add(g11); f1->add(g12); f1->add(edgenew);
				f2->add(g21); f2->add(g22); f2->add(g23); f2->add(edgenew);

				edgenew->add(f1); edgenew->add(f2);
				g11->del(self); g11->add(f1);
				g12->del(self); g12->add(f1);
				g21->del(self); g21->add(f2);
				g22->del(self); g22->add(f2);
				g23->del(self); g23->add(f2);

				allEntities.add(edgenew);
				allEntities.add(f1);
				allEntities.add(f2);

				allEntities.del(self);

		}
		else ;
		
}

void mDGMesh::combineTQ(Face* q1, Face* t1, Edge* e1,Face* q2, Face* t2, Edge* e2){
	//q2 must be quad	
	Vertex* w1[3];
		Edge* ee1[3];
		int a1;
		for(int i=0;i!=3;i++){
				w1[i]=(Vertex*)t1->get(0,i);
				ee1[i]=(Edge*)t1->get(1,i);
				if(w1[i]!=(Vertex*)e1->get(0,0)&& w1[i]!=(Vertex*)e1->get(0,1))
										  a1=i;
					               }
		if(q1->getType()==mEntity::QUAD){
		    Vertex* w2[4];
				Edge* ee2[4];
				//Vertex* spec=rightend(q1,e1);
				int a2;
				for(int i=0;i!=4;i++){
					  w2[i]=(Vertex*)q1->get(0,i);
						ee2[i]=(Edge*)q1->get(1,i);
						if(ee2[i]==e1)a2=i;
				}
        Vertex* w3[3];
		    Edge* ee3[3];
		    int a3;
		    for(int i=0;i!=3;i++){
				    w3[i]=(Vertex*)t2->get(0,i);
				    ee3[i]=(Edge*)t2->get(1,i);
				    if(w3[i]!=(Vertex*)e2->get(0,0)&& w3[i]!=(Vertex*)e2->get(0,1))
										  a3=i;
					                   }
				Vertex* w4[4];
				Edge* ee4[4];
				//Vertex* spec=rightend(q1,e1);
				int a4;
				for(int i=0;i!=4;i++){
					  w4[i]=(Vertex*)q2->get(0,i);
						ee4[i]=(Edge*)q2->get(1,i);
						if(ee4[i]==e2)a4=i;
				}
				if(w1[(a1+1)%3]==w3[(a3+2)%3]){
				Vertex* n1w[4];
				n1w[0]=w1[a1];
				n1w[1]=w2[(a2+2)%4];
				n1w[2]=w2[(a2+3)%4];
				n1w[3]=w2[(a2)%4];

				Vertex* n2w[4];
				n2w[0]=w3[a3];
				n2w[1]=w4[(a4+1)%4];
				n2w[2]=w4[(a4+2)%4];
				n2w[3]=w4[(a4+3)%4];

				Edge* enew= new Edge (n1w[0],n1w[1], q1->getClassification());

				Face* f1=new Face(n1w[0],n1w[1],n1w[2],n1w[3],q1->getClassification());
				Face* f2=new Face(n2w[0],n2w[1],n2w[2],n2w[3],q2->getClassification());

				f1->add(enew); f1->add(ee2[(a2+2)%4]); f1->add(ee2[(a2+3)%4]); f1->add(ee1[(a1+2)%3]);
				f2->add(ee3[a3]); f2->add(ee4[(a4+1)%4]); f2->add(ee4[(a4+2)%4]); f2->add(enew);

				enew->add(f1); enew->add(f2);
				ee1[(a1+2)%3]->del(t1); ee1[(a1+2)%3]->add(f1);
				ee2[(a2+2)%4]->del(q1); ee2[(a2+2)%4]->add(f1);
				ee2[(a2+3)%4]->del(q1); ee2[(a2+3)%4]->add(f1);
				ee3[(a3)%3]->del(t2); ee3[(a3)%3]->add(f2);
				ee4[(a4+1)%4]->del(q2); ee4[(a4+1)%4]->add(f2);
				ee4[(a4+2)%4]->del(q2); ee4[(a4+2)%4]->add(f2);

				allEntities.add(enew);
				allEntities.add(f1);
				allEntities.add(f2);

				allEntities.del(t1);
				allEntities.del(t2);
				allEntities.del(q1);
				allEntities.del(q2);
				allEntities.del(ee2[(a2+1)%4]);
				allEntities.del(ee1[a1]);
				//allEntities.del(w1[(a1+1)%3]);
				allEntities.del(ee1[(a1+1)%3]);
				allEntities.del(ee3[(a3+1)%3]);

				}

				else if(w1[(a1+2)%3]==w3[(a3+1)%3]){
				Vertex* n1w[4];
				n1w[0]=w1[a1];
				n1w[1]=w2[(a2+1)%4];
				n1w[2]=w2[(a2+2)%4];
				n1w[3]=w2[(a2+3)%4];

				Vertex* n2w[4];
				n2w[0]=w3[a3];
				n2w[1]=w4[(a4+2)%4];
				n2w[2]=w4[(a4+3)%4];
				n2w[3]=w4[(a4)%4];

				Edge* enew= new Edge (n2w[0],n2w[1], q2->getClassification());

				Face* f1=new Face(n1w[0],n1w[1],n1w[2],n1w[3],q1->getClassification());
				Face* f2=new Face(n2w[0],n2w[1],n2w[2],n2w[3],q2->getClassification());

				f2->add(enew); f2->add(ee4[(a4+2)%4]); f2->add(ee4[(a4+3)%4]); f2->add(ee3[(a3+2)%3]);
				f1->add(ee1[a1]); f1->add(ee2[(a2+1)%4]); f1->add(ee2[(a2+2)%4]); f1->add(enew);

				enew->add(f1); enew->add(f2);
				ee1[(a1)%3]->del(t1); ee1[(a1)%3]->add(f1);
				ee2[(a2+1)%4]->del(q1); ee2[(a2+1)%4]->add(f1);
				ee2[(a2+2)%4]->del(q1); ee2[(a2+2)%4]->add(f1);
				ee3[(a3+2)%3]->del(t2); ee3[(a3+2)%3]->add(f2);
				ee4[(a4+2)%4]->del(q2); ee4[(a4+2)%4]->add(f2);
				ee4[(a4+3)%4]->del(q2); ee4[(a4+3)%4]->add(f2);

				allEntities.add(enew);
				allEntities.add(f1);
				allEntities.add(f2);

				allEntities.del(t1);
				allEntities.del(t2);
				allEntities.del(q1);
				allEntities.del(q2);
				allEntities.del(ee2[(a2+3)%4]);
				allEntities.del(ee3[a3]);
			//	allEntities.del(w3[(a3+1)%3]);
				allEntities.del(ee1[(a1+1)%3]);
				allEntities.del(ee3[(a3+1)%3]);

				}
				else{printf("situation not considered\n");}



		}
		else if(q1->getType()==mEntity::PANT){
			 Vertex* w2[5];
				Edge* ee2[5];
				int a2;
				for(int i=0;i!=5;i++){
					  w2[i]=(Vertex*)q1->get(0,i);
						ee2[i]=(Edge*)q1->get(1,i);
						if(ee2[i]==e1)a2=i;
				}
        Vertex* w3[3];
		    Edge* ee3[3];
		    int a3;
		    for(int i=0;i!=3;i++){
				    w3[i]=(Vertex*)t2->get(0,i);
				    ee3[i]=(Edge*)t2->get(1,i);
				    if(w3[i]!=(Vertex*)e2->get(0,0)&& w3[i]!=(Vertex*)e2->get(0,1))
										  a3=i;
					                   }
				Vertex* w4[4];
				Edge* ee4[4];
				int a4;
				for(int i=0;i!=4;i++){
					  w4[i]=(Vertex*)q2->get(0,i);
						ee4[i]=(Edge*)q2->get(1,i);
						if(ee4[i]==e2)a4=i;
				}
				
				if(w1[(a1+1)%3]==w3[(a3+2)%3]){
				Vertex* n1w[4];
				n1w[0]=w1[a1];
				n1w[1]=w2[(a2+2)%5];
				n1w[2]=w2[(a2+3)%5];
				n1w[3]=w2[(a2+4)%5];

				Vertex* n2w[4];
				n2w[0]=w3[a3];
				n2w[1]=w4[(a4+1)%4];
				n2w[2]=w4[(a4+2)%4];
				n2w[3]=w4[(a4+3)%4];

				Edge* enew= new Edge (n1w[0],n1w[1], q1->getClassification());
				Edge* enew1=new Edge (n1w[0],n1w[3], q1->getClassification());

				Face* f1=new Face(n1w[0],n1w[1],n1w[2],n1w[3],q1->getClassification());
				Face* f2=new Face(n2w[0],n2w[1],n2w[2],n2w[3],q2->getClassification());

				f1->add(enew); f1->add(ee2[(a2+2)%5]); f1->add(ee2[(a2+3)%5]); f1->add(enew1);
				f2->add(ee3[a3]); f2->add(ee4[(a4+1)%4]); f2->add(ee4[(a4+2)%4]); f2->add(enew);

				enew->add(f1); enew->add(f2);
				enew1->add(f1);
				//ee1[(a1+2)%3]->del(t1); ee1[(a1+2)%3]->add(f1);
				ee2[(a2+2)%5]->del(q1); ee2[(a2+2)%5]->add(f1);
				ee2[(a2+3)%5]->del(q1); ee2[(a2+3)%5]->add(f1);
				ee3[(a3)%3]->del(t2); ee3[(a3)%3]->add(f2);
				ee4[(a4+1)%4]->del(q2); ee4[(a4+1)%4]->add(f2);
				ee4[(a4+2)%4]->del(q2); ee4[(a4+2)%4]->add(f2);

				allEntities.add(enew);
				allEntities.add(enew1);
				allEntities.add(f1);
				allEntities.add(f2);

				allEntities.del(t1);
				allEntities.del(t2);
				allEntities.del(q1);
				allEntities.del(q2);
				allEntities.del(ee2[(a2+1)%5]);
				allEntities.del(ee1[a1]);
				//allEntities.del(w1[(a1+1)%3]);
				//allEntities.del(w1[(a1+2)%3]);
				allEntities.del(ee1[(a1+1)%3]);
				allEntities.del(ee1[(a1+2)%3]);
				allEntities.del(ee2[(a2+4)%5]);
				allEntities.del(ee4[(a4)%4]);

				}

				else if(w1[(a1+2)%3]==w3[(a3+1)%3]){
				Vertex* n1w[4];
				n1w[0]=w1[a1];
				n1w[1]=w2[(a2+2)%5];
				n1w[2]=w2[(a2+3)%5];
				n1w[3]=w2[(a2+4)%5];

				Vertex* n2w[4];
				n2w[0]=w3[a3];
				n2w[1]=w4[(a4+2)%4];
				n2w[2]=w4[(a4+3)%4];
				n2w[3]=w4[(a4)%4];

				Edge* enew= new Edge (n2w[0],n2w[1], q2->getClassification());
				Edge* enew1= new Edge (n1w[0],n1w[1], q1->getClassification());

				Face* f1=new Face(n1w[0],n1w[1],n1w[2],n1w[3],q1->getClassification());
				Face* f2=new Face(n2w[0],n2w[1],n2w[2],n2w[3],q2->getClassification());

				f2->add(enew); f2->add(ee4[(a4+2)%4]); f2->add(ee4[(a4+3)%4]); f2->add(ee3[(a3+2)%3]);
				f1->add(enew1); f1->add(ee2[(a2+2)%5]); f1->add(ee2[(a2+3)%5]); f1->add(enew);

				enew->add(f1); enew->add(f2);
				enew1->add(f1);
				//ee1[(a1)%3]->del(t1); ee1[(a1)%3]->add(f1);
				ee2[(a2+2)%5]->del(q1); ee2[(a2+2)%5]->add(f1);
				ee2[(a2+3)%5]->del(q1); ee2[(a2+3)%5]->add(f1);
				ee3[(a3+2)%3]->del(t2); ee3[(a3+2)%3]->add(f2);
				ee4[(a4+2)%4]->del(q2); ee4[(a4+2)%4]->add(f2);
				ee4[(a4+3)%4]->del(q2); ee4[(a4+3)%4]->add(f2);

				allEntities.add(enew);
				allEntities.add(enew1);
				allEntities.add(f1);
				allEntities.add(f2);

				allEntities.del(t1);
				allEntities.del(t2);
				allEntities.del(q1);
				allEntities.del(q2);
				allEntities.del(ee2[(a2+3)%4]);
				allEntities.del(ee3[a3]);
				//allEntities.del(w3[(a3+1)%3]);
				allEntities.del(ee1[(a1+1)%3]);
				allEntities.del(ee3[(a3+1)%3]);
				//allEntities.del(w1[(a1+1)%3]);
				allEntities.del(ee1[(a1)%3]);
				allEntities.del(ee2[(a2+1)%5]);
				allEntities.del(ee2[(a2)%5]);
				allEntities.del(ee4[(a4)%4]);

				}
				else{printf("situation not considered\n");}





		}
		else printf("situation not considered\n");

}

void mDGMesh::dealQuad(int dim){
    Edge *tempstore[SIZEN];
	int count=0;
	iter it;
	iter end_m1 =end (dim-1);
	for (it= begin(dim-1);it !=end_m1;++it)
		{
			mEntity *e=*it;
			if(e->getClassification()->getId()==40000||e->getClassification()->getId()==20000)
			{ 
				Edge* f0=(Edge*) e;
				tempstore[count]=f0;
				f0->print();
				count++;
			}
	}
	for(int count1=0;count1<count;++count1){
		Face *f0=(Face*)tempstore[count1]->get(2,0);
		//f0->print();
		if( f0->getType()==mEntity::QUAD){
			Edge* e[4];
		    Vertex* w[4];
		    int a;
		    for(int i=0;i!=4;i++){
				e[i]=(Edge*)f0->get(1,i);
				w[i]=(Vertex*)f0->get(0,i);
				if (e[i]->getClassification()->getId()==40000||e[i]->getClassification()->getId()==20000) a=i;
		    }
			double length01=w[(a+1)%4]->distance(w[(a+2)%4]);
			double lg=longest(f0);
			double ratio=length01/lg;
			if (ratio<EPS){
			//    Face *q1=findnb(f0, e[(a+2)%4]);
				splitQuaddial(f0,w[(a+1)%4]);
				Face *t1, *q1;
				if(e[(a+2)%4]->get(2,0)->getType()==mEntity::TRI){
				    t1=(Face*)e[(a+2)%4]->get(2,0);
				    q1=(Face*)e[(a+2)%4]->get(2,1);
				}
				else{
					t1=(Face*)e[(a+2)%4]->get(2,1);
					q1=(Face*)e[(a+2)%4]->get(2,0);
				}
				Face *t2=findnb(t1,e[(a+1)%4]);
				if(t2->getType()==mEntity::TRI){
					Edge* sube[3];
					    Vertex* subw[3];
							int bd;
					    for(int i=0;i!=3;i++){
					        sube[i]=(Edge*)t2->get(1,i);
							subw[i]=(Vertex*)t2->get(0,i);
							if (sube[i]==e[(a+1)%4])bd=i;
						}
						Face *q2=findnb(t2, sube[(bd+2)%3]);
						combineTQ(q2, t2, sube[(bd+2)%3], q1,t1,e[(a+2)%4]);
				}
				else if(t2->getType()==mEntity::QUAD){
					Edge* sube[4];
					    Vertex* subw[4];
							int bd;
					    for(int i=0;i!=4;i++){
					        sube[i]=(Edge*)t2->get(1,i);
							subw[i]=(Vertex*)t2->get(0,i);
							if (sube[i]==e[(a+1)%4])bd=i;
						}
						Face *q2=findnb(t2, sube[(bd+3)%4]);
						splitQuaddial(t2,subw[(bd+1)%4]);
						t2=findnb(q2,sube[(bd+3)%4]);
						combineTQ(q2, t2, sube[(bd+3)%4], q1,t1,e[(a+2)%4]);
				}
				else if(t2->getType()==mEntity::PANT){
					Edge* sube[5];
					    Vertex* subw[5];
							int bd;
					    for(int i=0;i!=5;i++){
					        sube[i]=(Edge*)t2->get(1,i);
							subw[i]=(Vertex*)t2->get(0,i);
							if (sube[i]==e[(a+1)%4])bd=i;
						}
						Face *q2=findnb(t2, sube[(bd+4)%5]);
						splitPant(t2,subw[(bd+1)%5],subw[(bd+4)%5]);
						t2=findnb(q2,sube[(bd+4)%5]);
						combineTQ(q2, t2, sube[(bd+4)%5], q1,t1,e[(a+2)%4]);


				}
				else {
					printf("not considered\n");
				}


			
			}
			else{
			//	w[a]->print();
				//w[(a+3)%4]->print();
			    double length02=w[a]->distance(w[(a+3)%4]);
				ratio=length02/lg;
				if (ratio<EPS){
					//Face *q1=findnb(f0, e[(a+2)%4]);
				    splitQuaddial(f0,w[(a)%4]);
					Face *t1, *q1;
				    if(e[(a+2)%4]->get(2,0)->getType()==mEntity::TRI){
				        t1=(Face*)e[(a+2)%4]->get(2,0);
				        q1=(Face*)e[(a+2)%4]->get(2,1);
				    }
				    else{
					    t1=(Face*)e[(a+2)%4]->get(2,1);
					    q1=(Face*)e[(a+2)%4]->get(2,0);
				    }
				    Face *t2=findnb(t1,e[(a+3)%4]);
				    if(t2->getType()==mEntity::TRI){
					    Edge* sube[3];
					    Vertex* subw[3];
						int bd;
					    for(int i=0;i!=3;i++){
					        sube[i]=(Edge*)t2->get(1,i);
							subw[i]=(Vertex*)t2->get(0,i);
							if (sube[i]==e[(a+3)%4])bd=i;
						}
						Face *q2=findnb(t2, sube[(bd+1)%3]);
						combineTQ(q2, t2, sube[(bd+1)%3],q1,t1,e[(a+2)%4]);
				    }
				    else if(t2->getType()==mEntity::QUAD){
					    Edge* sube[4];
					    Vertex* subw[4];
						int bd;
					    for(int i=0;i!=4;i++){
					        sube[i]=(Edge*)t2->get(1,i);
							subw[i]=(Vertex*)t2->get(0,i);
							if (sube[i]==e[(a+3)%4])bd=i;
						}
						Face *q2=findnb(t2, sube[(bd+1)%4]);
						splitQuaddial(t2,subw[(bd)%4]);
						t2=findnb(q2,sube[(bd+1)%4]);
						combineTQ(q2, t2, sube[(bd+1)%4], q1,t1,e[(a+2)%4]);
				    }
				    else if(t2->getType()==mEntity::PANT){
					    Edge* sube[5];
					    Vertex* subw[5];
						int bd;
					    for(int i=0;i!=5;i++){
					        sube[i]=(Edge*)t2->get(1,i);
							subw[i]=(Vertex*)t2->get(0,i);
							if (sube[i]==e[(a+3)%4])bd=i;
						}
						Face *q2=findnb(t2, sube[(bd+1)%5]);
						splitPant(t2,subw[(bd)%5],subw[(bd+2)%5]);
						t2=findnb(q2,sube[(bd+1)%5]);
						combineTQ(q2, t2, sube[(bd+1)%5], q1,t1,e[(a+2)%4]);


				    }
				    else {
					    printf("not considered\n");
				    }

				}
				else{}
			}
		}

}
}
void mDGMesh::merge(int dim){
	mImportExport io;
	Edge *tempstore[SIZEN];
	int count=0;
	iter it;
	iter end_m1 =end (dim-1);
	for (it= begin(dim-1);it !=end_m1;++it)
		{
			mEntity *e=*it;
			if(e->getClassification()->getId()==40000||e->getClassification()->getId()==20000)
			{ 
				Edge* f0=(Edge*) e;
				tempstore[count]=f0;
				count++;
			}
	}
	for(int count1=0;count1<count;++count1){
		Face *f0=(Face*)tempstore[count1]->get(2,0);
		//f0->print();
				switch (f0->getType()){
					case mEntity::TRI:{
							Edge* e[3];
					    Vertex* w[3];
							int bd;
					    for(int i=0;i!=3;i++){
					        e[i]=(Edge*)f0->get(1,i);
							    w[i]=(Vertex*)f0->get(0,i);
									if(e[i]->getClassification()->getId()==40000||e[i]->getClassification()->getId()==20000)
										  bd=i;
					                         }
              Face* nb01=findnb(f0,e[(bd+1)%3]);
							Face* nb02=findnb(f0,e[(bd+2)%3]);

							double length01=w[(bd+1)%3]->distance(w[(bd+2)%3]);
							double length02=w[(bd)%3]->distance(w[(bd+2)%3]);
							double lengthnb01=longest(nb01);
							double lengthnb02=longest(nb02);
							double ratio01=length01/lengthnb01;
							double ratio02=length02/lengthnb02;

							if (ratio01<EPS||ratio02<EPS){
								if(length01<length02){
									if(nb01->getType()==mEntity::PANT){
								    Vertex* subw[5];
										Edge* sube[5];
										int jd;
										for (int i=0;i!=5;i++){
									 	    subw[i]=(Vertex*)nb01->get(0,i);
												sube[i]=(Edge*)nb01->get(1,i);
										    if(subw[i]==w[(bd+2)%3])jd=i;
										}
										Face* nb03=findnb(nb01,sube[(jd+4)%5]);
										Edge* tempee=sube[(jd+4)%5];
										splitPant(nb01,subw[(jd+1)%5],subw[(jd+4)%5]);
										nb01=findnb(f0,e[(bd+1)%3]);
										combineTQ(nb02,f0,e[(bd+2)%3],nb03,nb01,tempee);
									}
									else if(nb01->getType()==mEntity::QUAD){
										Vertex* subw[4];
										Edge* sube[4];
										int jd;
										for (int i=0;i!=4;i++){
									 	    subw[i]=(Vertex*)nb01->get(0,i);
												sube[i]=(Edge*)nb01->get(1,i);
										    if(subw[i]==w[(bd+2)%3])jd=i;
										}
										Face* nb03=findnb(nb01,sube[(jd+3)%4]);
										Edge* tempee=sube[(jd+3)%4];
										splitQuaddial(nb01,subw[(jd+1)%4]);
										//io.exportGmshFile("check02.msh",this);
										nb01=findnb(f0,e[(bd+1)%3]);
										combineTQ(nb02,f0,e[(bd+2)%3],nb03,nb01,tempee);}
									else if(nb01->getType()==mEntity::TRI){
										Vertex* subw[3];
										Edge* sube[3];
										int jd;
										for (int i=0;i!=3;i++){
									 	    subw[i]=(Vertex*)nb01->get(0,i);
												sube[i]=(Edge*)nb01->get(1,i);
										    if(subw[i]==w[(bd+2)%3])jd=i;
										}
										Face* nb03=findnb(nb01,sube[(jd+2)%3]);
										Edge* tempee=sube[(jd+2)%3];
										//splitQuaddial(nb01,subw[(jd+1)%4]);
										//io.exportGmshFile("check02.msh",this);
										//nb01=findnb(f0,e[(bd+1)%3]);
										combineTQ(nb02,f0,e[(bd+2)%3],nb03,nb01,tempee);}
									else {
										printf("not considered\n");
									}
								  //  io.exportGmshFile("check02.msh",this);
								}
								else{
									if(nb02->getType()==mEntity::PANT){
								    Vertex* subw[5];
										Edge* sube[5];
										int jd;
										for (int i=0;i!=5;i++){
									 	    subw[i]=(Vertex*)nb02->get(0,i);
												sube[i]=(Edge*)nb02->get(1,i);
										    if(subw[i]==w[(bd+2)%3])jd=i;
										}
										Face* nb03=findnb(nb02,sube[(jd)%5]);
										Edge* tempee=sube[(jd)%5];
										splitPant(nb02,subw[(jd+1)%5],subw[(jd+4)%5]);
										nb02=findnb(f0,e[(bd+2)%3]);
										combineTQ(nb01,f0,e[(bd+1)%3],nb03,nb02,tempee);
									}
									else if(nb02->getType()==mEntity::QUAD){
										Vertex* subw[4];
										Edge* sube[4];
										int jd;
										for (int i=0;i!=4;i++){
									 	    subw[i]=(Vertex*)nb02->get(0,i);
												sube[i]=(Edge*)nb02->get(1,i);
										    if(subw[i]==w[(bd+2)%3])jd=i;
										}
										Face* nb03=findnb(nb02,sube[(jd)%4]);
										Edge* tempee=sube[(jd)%4];
										splitQuaddial(nb02,subw[(jd+1)%4]);
										nb02=findnb(f0,e[(bd+2)%3]);
										combineTQ(nb01,f0,e[(bd+1)%3],nb03,nb02,tempee);}
									else if(nb02->getType()==mEntity::TRI){
										Vertex* subw[3];
										Edge* sube[3];
										int jd;
										for (int i=0;i!=3;i++){
									 	    subw[i]=(Vertex*)nb02->get(0,i);
												sube[i]=(Edge*)nb02->get(1,i);
										    if(subw[i]==w[(bd+2)%3])jd=i;
										}
										Face* nb03=findnb(nb02,sube[(jd)%3]);
										Edge* tempee=sube[(jd)%3];
										//splitQuaddial(nb01,subw[(jd+1)%4]);
										//io.exportGmshFile("check02.msh",this);
										//nb01=findnb(f0,e[(bd+1)%3]);
										combineTQ(nb01,f0,e[(bd+1)%3],nb03,nb02,tempee);}
									else {
										printf("not considered\n");
									}
								}
							}
							else;

					    						
														}break;
					default:break;
				}
				
			
		}
}

void mDGMesh::resplitPant(int dim){
	Face *tempstore[SIZEN];
	int count=0;
	iter it;
	iter end_m1 =end (dim);
	for (it= begin(dim);it !=end_m1;++it)
		{
			Face *e=(Face*) *it;
			if(e->getType()==mEntity::PANT)
			{ 
				tempstore[count]=e;
				count++;
			}
	}
    /*iter it;
		iter end_m1 = end (dim);
		for (it= begin(dim);it !=end_m1;++it)
		{*/
	for (int count1=0; count1<count; count1++){
			Face* e=tempstore[count1];
			//mEntity *e=*it;
			//if(e->getClassification()->getId()==40000){
			if(e->getType()==mEntity::PANT){
				  //Face* f0=(Face*)e->get(2,0);
			  Face* f0=e;  
				Vertex* ee[5];
				Edge* ww[5];
				int jd;
					for (int i=0;i!=5;i++){
					    ee[i]=(Vertex*)f0->get(0,i);
							ww[i]=(Edge*)f0->get(1,i);
							if(ww[i]->getClassification()->getId()==40000||ww[i]->getClassification()->getId()==20000)
								jd=i;
					}
					splitPant(f0,ee[(jd+2)%5],ee[(jd+4)%5]);
			}
		}
}

void mDGMesh::resplitQuad(int dim){
	Face *tempstore[SIZEN];
	int count=0;
	iter it;
	iter end_m1 =end (dim);
	for (it= begin(dim);it !=end_m1;++it)
		{
			Face *e=(Face*) *it;
			if(e->getType()==mEntity::QUAD)
			{ 
				tempstore[count]=e;
				count++;
			}
	}
    /*iter it;
		iter end_m1 = end (dim);
		for (it= begin(dim);it !=end_m1;++it)
		{*/
	for (int count1=0; count1<count; count1++){
			Face* e=tempstore[count1];
			//mEntity *e=*it;
			//if(e->getClassification()->getId()==40000){
			if(e->getType()==mEntity::QUAD){
				  //Face* f0=(Face*)e->get(2,0);
			  Face* f0=e;  
				f0->print();
				if (f0->get(0,0)->getId()==463){
					printf("stop");
				}
				/*Vertex* ee[4];
				Edge* ww[4];
				int jd;
					for (int i=0;i!=5;i++){
					    ee[i]=(Vertex*)f0->get(0,i);
							ww[i]=(Edge*)f0->get(1,i);
							if(ww[i]->getClassification()->getId()==40000||ww[i]->getClassification()->getId()==20000)
								jd=i;
					}*/
					splitQuaddial(f0,(Vertex*)f0->get(0,0));
			}
		}
}

void mDGMesh::createConnections2(int to)
{
	//adding a pointer from edge to parent edge)
	  vector<mEntity *> tmpvec;
	  iter end_j=end(to);
      for(iter it=begin(to);it!=end_j;it++)
	{
	  mEntity *e = (*it);  
	  iter end_i=end(to);
	  for(iter im=it;im!=end_i;im++)
		  if(im!=it){
			  mEntity *n = (*im);
			  
			  Vertex *vn[2];
			  vn[0]=(Vertex*)n->get(0,0);
              vn[1]=(Vertex*)n->get(0,1);
			  
			  Vertex *ve[2];
			  ve[0]=(Vertex*)e->get(0,0);
              ve[1]=(Vertex*)e->get(0,1);

			  mPoint vnp0=vn[0]->point();
			  mPoint vnp1=vn[1]->point();
			  mPoint vep0=ve[0]->point();
			  mPoint vep1=ve[1]->point();

			  if(lineinside(vep0, vep1, vnp0)&&lineinside(vep0, vep1, vnp1)){
			      e->add(n);
				  n->setParent(e);
					Face *f0 = (Face*)e->get(2,0);
					n->add(f0);
					int r =  getRefinementLevel(e) + 1;
          attachRefinementLevel( n, r);
		        int vlength=tmpvec.size();
				int tmp=0;
				for (int ii=0;ii<vlength;ii++){
                    if (tmpvec[ii]==e)
						tmp=1;

				}
				if (tmp==0)
					tmpvec.push_back(e);
				//	if(e->isAdjacencyCreated(1));
				
				
           // allSplittedEntities.add(e);
					//	allEntities.del(e);
				}
			  
			  else if(lineinside(vnp0, vnp1, vep0)&&lineinside(vnp0, vnp1, vep1)){
				  n->add(e);
				  e->setParent(n);
					Face *f0 = (Face*)n->get(2,0);
					e->add(f0);
					int r =  getRefinementLevel(n) + 1;
                    attachRefinementLevel( e, r);
		        int vlength=tmpvec.size();
				int tmp=0;
				for (int ii=0;ii<vlength;ii++){
                    if (tmpvec[ii]==n)
						tmp=1;

				}
				if (tmp==0)
					tmpvec.push_back(n);
				//	if(n->isAdjacencyCreated(1));
				
				
						
          //  allSplittedEntities.add(n);
					//	allEntities.del(n);
			  }
			  else;
	      
		  }
	}
	int vlength=tmpvec.size();
				
				for (int ii=0;ii<vlength;ii++){
                    allEntities.del(tmpvec[ii]);

				}
				

}

int mDGMesh::lineinside(mPoint a, mPoint b, mPoint c){
	if((abs(c(0)-a(0))<exp(-10.)&&abs(c(1)-a(1))<exp(-10.))||(abs(c(0)-b(0))<exp(-10.)&&abs(c(1)-b(1))<exp(-10.))){
		return 1;
	}
	else if(abs(a(0)-b(0))<exp(-10.)){
		if (abs(c(0)-a(0))<exp(-10.)){
			double k2=(c(1)-a(1))/(b(1)-a(1));
			if(-exp(-10.)<k2&&k2<1+exp(-10.))return 1;
			else return 0;
		}
		else return 0;
	}
	else if(abs(a(1)-b(1))<exp(-10.)){
     if (abs(c(1)-a(1))<exp(-10.)){
			double k1=(c(0)-a(0))/(b(0)-a(0));
			if(-exp(-10.)<k1&&k1<1+exp(-10.))return 1;
			else return 0;
		}
		else return 0;
	}
	else{
    double k1=(c(0)-a(0))/(b(0)-a(0));
	double k2=(c(1)-a(1))/(b(1)-a(1));
	if (abs(k1-k2)<exp(-10.)&&-exp(-10.)<k1&&k1<1+exp(-10.))return 1;
	else return 0;
	}
}


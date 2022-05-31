#include "mImportExport.h"
#include "mMesh.h"
#include "Vertex.h"
#include "Edge.h"
#include "Face.h"
#include "mTet.h"
#include "mException.h"
#include <stdio.h>
#include <string.h>

void classifyUnclassifiedVerices (mMesh *m)
{
  for(int i=1;i<4;i++)
    {
      //for(mMesh::iter it = m->beginall(i);it != m->endall(i) ; ++it)
      for(mMesh::iter it = m->begin(i);it != m->end(i) ; ++it)
	{
	  const mEntity *e = *it;
	  int size = e->size(0);
	  for(int j=0;j<size;j++)
	    {
	      mEntity *v = e->get(0,j);
	      if(!v->getClassification())v->classify(e->getClassification());
	    }
	}
    }
}

void mImportExport::import (char *fName, mMesh *theMesh)
{
  char ext[6];
  strcpy(ext,fName+(strlen(fName)-4));  //it is assumed here that the name of the mesh file is in a format .***
  if(!strcmp(ext,".sms"))
    {
      importSmsFile(fName, theMesh);      
    }
  else if(!strcmp(ext,".msh"))
    {
      importGmshFile(fName, theMesh);
    }
  else
    {
      char text[256];
      sprintf(text,"unknown extension %s in file %s",ext,fName);
      throw new mExceptionFileNotFound (__LINE__,__FILE__,text);
    }
}

void mImportExport::importGmshFile (char *fName, mMesh *theMesh)
{
  FILE *f = fopen (fName,"r");
  if(!f)throw new mExceptionFileNotFound (__LINE__,__FILE__,"impossible to import the gmsh file");

  char line[256];
  while(!feof(f))
    {
      fscanf(f,"%s",line);
      if(!strcmp(line,"$NOE") ||!strcmp(line,"$NOD"))  //reading vertices 
	{
	  int NbNod;
	  fscanf(f,"%d",&NbNod);             //number of vertices in the mesh
	  printf("%d nodes\n",NbNod);
	  theMesh->resize(0,NbNod);
	  for(int i=0;i<NbNod;i++)
	    {
	      int iNod;
	      double x,y,z;
	      fscanf(f,"%d %lf %lf %lf",&iNod,&x,&y,&z);
	      theMesh->createVertex(iNod,x,y,z,0);
	    }
	}
      if(!strcmp(line,"$ELM"))
	{
	  int NbElm;
	  fscanf(f,"%d",&NbElm);
	  printf("%d elms\n",NbElm);   //gmsh writes in $ELM the number of cells + the number of edges marked for some purpose
	  for(int i=0;i<NbElm;i++)
	    { 
			

	      int iNbNod,iTyp,iGrp,iElm,iZon;
	      fscanf(f,"%d %d %d %d %d",&iElm,&iTyp,&iGrp,&iZon,&iNbNod);
	      switch(iTyp)
		{
		  case 1 : // edge
		  {
		    int iNod1,iNod2;
		    fscanf(f,"%d %d",&iNod1,&iNod2);
		    theMesh->add(theMesh->createEdge(iNod1,iNod2,theMesh->getGEntity(iGrp,1)));
		  }
		  break;		
		case 2 : // triangle
		  {
		    int iNod1,iNod2,iNod3;
		    fscanf(f,"%d %d %d",&iNod1,&iNod2,&iNod3);
		    theMesh->add(theMesh->createFaceWithVertices(iNod1, iNod2,iNod3,theMesh->getGEntity(iGrp,2)));
		  }
		  break;		
		case 3 : // quad
		  {
		    int iNod1,iNod2,iNod3,iNod4;
		    fscanf(f,"%d %d %d %d",&iNod1,&iNod2,&iNod3,&iNod4);
		    theMesh->add(theMesh->createFaceWithVertices(iNod1, iNod2,iNod3, iNod4, theMesh->getGEntity(iGrp,2)));

		  }
		  break;
		case 4 : // tetraedron
		  {
		    int iNod1,iNod2,iNod3,iNod4;
		    fscanf(f,"%d %d %d %d",&iNod1,&iNod2,&iNod3,&iNod4);
		    theMesh->add(theMesh->createTetWithVertices(iNod1, iNod2,iNod3,iNod4,theMesh->getGEntity(iGrp,3)));
		  }
		  break;
		case 5 : // Hex
		  {
		    int iNod1,iNod2,iNod3,iNod4,iNod5,iNod6,iNod7,iNod8;
		    fscanf(f,"%d %d %d %d %d %d %d %d",&iNod1,&iNod2,&iNod3,&iNod4,&iNod5,&iNod6,&iNod7,&iNod8);
		    theMesh->createHexWithVertices(iNod1, iNod2,iNod3,iNod4,
						   iNod5, iNod6,iNod7,iNod8,
						   theMesh->getGEntity(iGrp,3));
		  }
		  break;
		  case 6 : // pantagon
		  {
		    int iNod1,iNod2,iNod3,iNod4,iNod5;
		    fscanf(f,"%d %d %d %d %d",&iNod1,&iNod2,&iNod3,&iNod4,&iNod5);
		    theMesh->add(theMesh->createFaceWithVertices(iNod1, iNod2,iNod3, iNod4, iNod5, theMesh->getGEntity(iGrp,2)));

		  }
		  break;
		default:
			printf("This type of element is not recognized\n");
			exit(0);
		}
	    }
	}
    }
  fclose(f);
  classifyUnclassifiedVerices (theMesh);
}

// Same with streams
void mImportExport::importGmshFile (istream &f, mMesh *theMesh)
{
  char line[255];
  while(f.getline(line,255))  
  {
	  if(!strcmp(line,"$NOE") ||!strcmp(line,"$NOD"))
	{
	  int NbNod;
	  f >> NbNod;
	  	printf("%d nodes\n",NbNod);
	  for(int i=0;i<NbNod;i++)
	    {
	      int iNod;
	      double x,y,z;
	      f >> iNod >> x >> y >> z;
	      theMesh->createVertex(iNod,x,y,z,0);
	    }
	}
      if(!strcmp(line,"$ELM"))
	{
	  int NbElm;
	  f >> NbElm; 
	  printf("%d elms\n",NbElm);
	  for(int i=0;i<NbElm;i++)
	    {
	      int iNbNod,iTyp,iGrp,iElm,iZon;
	      f >> iElm >> iTyp >> iGrp >> iZon >> iNbNod;
	      switch(iTyp)
		{
		// tetraedron
		case 4 : 
		  {
		    int iNod1,iNod2,iNod3,iNod4;
		    f >> iNod1 >> iNod2 >> iNod3 >> iNod4;
		    theMesh->createTetWithVertices(iNod1, iNod2,iNod3,iNod4,theMesh->getGEntity(iGrp,3));
		  }
		  break;
		// Hex
		case 5 : 
		  {
		    int iNod1,iNod2,iNod3,iNod4,iNod5,iNod6,iNod7,iNod8;
		    f >> iNod1 >> iNod2 >> iNod3 >> iNod4 >> iNod5 >> iNod6 >> iNod7 >> iNod8;
		    theMesh->createHexWithVertices(iNod1, iNod2,iNod3,iNod4,
						   iNod5, iNod6,iNod7,iNod8,
						   theMesh->getGEntity(iGrp,3));
		  }
		  break;
				// triangle
		case 2 : 
		  {
		    int iNod1,iNod2,iNod3;
		    f >> iNod1 >> iNod2 >> iNod3;
		    theMesh->createFaceWithVertices(iNod1, iNod2,iNod3,theMesh->getGEntity(iGrp,2));
		  }
		  break;
				// quad
		case 3 : 
		  {
		    int iNod1,iNod2,iNod3,iNod4;
		    f >> iNod1 >> iNod2 >> iNod3 >> iNod4;
		    theMesh->createFaceWithVertices(iNod1, iNod2,iNod3, iNod4, theMesh->getGEntity(iGrp,2));
		  }
		  break;
				// line
		case 1 : 
		  {
		    int iNod1,iNod2;
		    f >> iNod1 >> iNod2;
		    theMesh->createEdge(iNod1,iNod2,theMesh->getGEntity(iGrp,1));
		  }
		  break;
		case 15 : 
		  {
		    int iNod;
		    Vertex *v;
		    f >> iNod;
		    v = theMesh->getVertex(iNod); 
		    if(iGrp > 0)
		      {
			v->classify(theMesh->getGEntity(iGrp,0));
		      }		    
		  }
		  break;
		}
	    }
	}
    }
  classifyUnclassifiedVerices (theMesh);
}



void mImportExport::importSmsFile(char *fName, mMesh *theMesh)
{
  /*
  FILE *in = fopen (fName,"r");
  if(!in)throw new mExceptionFileNotFound (__LINE__,__FILE__,"unable to open file in loadSmsFile");

  char line[256];
  int i;
  
  int NbRegions,NbFaces,NbEdges,NbVertices,NbPoints,
    GEntityType,GEntityId,EntityNbConnections,Dummy,
    Edge1,Edge2,Edge3,Edge4,Face[4];
  int VertexId1,VertexId2,NbEdgesOnFace,NbFacesOnRegion;
  double x,y,z;

  fgets(line,255,in);
  fgets(line,255,in);
  sscanf(line,"%d %d %d %d %d",&NbRegions,&NbFaces,&NbEdges,&NbVertices,&NbPoints);
  printf("reading %d vertices ...\n",NbVertices);
  for(i=0;i<NbVertices;i++)
    {
      fgets(line,255,in);
      sscanf(line,"%d",&GEntityId); 
      if(GEntityId)
	{
	  sscanf(line,"%d %d %d",&GEntityId,&GEntityType,&EntityNbConnections); 
	  fgets(line,255,in);
	  sscanf(line,"%le %le %le",&x,&y,&z);
	  theMesh->createVertex(i+1,x,y,z,0);
	}
    }

  Edge **allEd = new Edge*[NbEdges];

  printf("reading %d edges ...\n",NbEdges);
  for(i=0;i<NbEdges;i++)
    {
      fgets(line,255,in);
      sscanf(line,"%d %d %d %d %d %d",&GEntityId,&GEntityType, &VertexId1,&VertexId2,&EntityNbConnections,&Dummy); 
      allEd[i] = theMesh->createEdge(VertexId1,VertexId2,theMesh->getGEntity(GEntityId,GEntityType));
    }
  
  Face **allFac = new 
    Face* [NbFaces];
  printf("reading %d faces ...\n",NbFaces);
  for(i=0;i<NbFaces;i++)
    {
      fgets(line,255,in);
      sscanf(line,"%d %d %d",&GEntityId,&GEntityType, &NbEdgesOnFace);
      if(NbEdgesOnFace == 3)
	{
	  sscanf(line,"%d %d %d %d %d %d",&GEntityId,&GEntityType, &NbEdgesOnFace,&Edge1,&Edge2,&Edge3);
	  allFac[i] = theMesh->createFaceWithEdges(allEd[abs(Edge1)-1],
						   allEd[abs(Edge2)-1],allEd[abs(Edge3)-1],					    	                                                  theMesh->getGEntity(GEntityId,GEntityType));		
	}
      else if(NbEdgesOnFace == 4)
	{
	  sscanf(line,"%d %d %d %d %d %d %d",&GEntityId,&GEntityType, &NbEdgesOnFace,&Edge1,&Edge2,&Edge3,&Edge4);

	  allFac[i] = theMesh->createFaceWithEdges(allEd[abs(Edge1)-1],allEd[abs(Edge2)-1],
						   allEd[abs(Edge3)-1],
						   allEd[abs(Edge4)-1],
						   theMesh->getGEntity(GEntityId,GEntityType));	
	}
    }
  delete [] allEd;
	
  printf("reading %d regions ...\n",NbRegions);
  for(i=0;i<NbRegions;i++)
    {
      fgets(line,255,in);
      sscanf(line,"%d %d",&GEntityId,&NbFacesOnRegion);
      if(NbFacesOnRegion == 4)
	{
	  sscanf(line,"%d %d %d %d %d %d",&GEntityId, &NbFacesOnRegion,&Face[0],&Face[1],&Face[2],&Face[3]);
	  theMesh->createTetWithFaces(allFac[abs(Face[0])-1],allFac[abs(Face[1])-1],
				      allFac[abs(Face[2])-1],
				      allFac[abs(Face[3])-1],
				      theMesh->getGEntity(GEntityId,GEntityType));
	  	  
	}
    }
  delete [] allFac;
  fclose (in);
  // classifyUnclassifiedVerices (theMesh);
  */
}

/*
void mImportExport::printSmsFile(char *fName)
{
  FILE *out = fopen (fName,"w");
  if(!out)throw new mExceptionFileNotFound (__LINE__,__FILE__);

  fprintf(out,"msh2sms 0.1\n");
  fprintf(out,"%d %d %d %d %d\n",allMRegions.size(),allMFaces.size(),
	  allMEdges.size(),allMVertices.size(),allMVertices.size());

  mEntityContainer<CMP>::const_iterator it    = allMVertices.begin();
  mEntityContainer<CMP>::const_iterator itEnd = allMVertices.end();
  int counter = 1;
  for(;it != itEnd; ++it)
    {
      Vertex *v = (Vertex*)(*it);
      fprintf(out,"%d %d %d\n",v->getClassification()->getId(),getClassificationNb(v),7);
      if(getClassificationNb(v) == 0 ||
	 getClassificationNb(v) == 3)
	fprintf(out,"%12.5E %12.5E %12.5E\n",v->point()(0),
		v->point()(1),v->point()(2));
      else if(getClassificationNb(v) == 1)
	fprintf(out,"%12.5E %12.5E %12.5E 0\n",v->point()(0),
		v->point()(1),v->point()(2));
      else if(getClassificationNb(v) == 2)
	fprintf(out,"%12.5E %12.5E %12.5E 0 0\n",v->point()(0),
		v->point()(1),v->point()(2));
      int *id = new int(counter);
      v->attachData("id",id);
      counter ++;
    }

  it    = allMEdges.begin();
  itEnd = allMEdges.end();
  counter = 1;
  for(;it != itEnd; ++it)
    {

      Edge *e = (Edge*)(*it);
      Vertex *v1 = (Vertex*)(e->get(mEntity::VERTEX,0));
      Vertex *v2 = (Vertex*)(e->get(mEntity::VERTEX,1));
      int *id1 = (int*)v1->getData("id");
      int *id2 = (int*)v2->getData("id");
      fprintf(out,"%d %d %d %d %d\n",e->getClassification()->getId(),getClassificationNb(e),*id1,*id2);
      int *id = new int(counter);
      e->attachData("id",id);
      counter ++;
    }


  it    = allMFaces.begin();
  itEnd = allMFaces.end();
  counter = 1;
  for(;it != itEnd; ++it)
    {
      Face *f = (Face*)(*it);
      Edge *e1 = (Edge*)(f->get(mEntity::EDGE,0));
      Edge *e2 = (Edge*)(f->get(mEntity::EDGE,1));
      Edge *e3 = (Edge*)(f->get(mEntity::EDGE,2));
      int id1 = (*((int*)e1->getData("id"))) * ((f->getUse(e1) == 1)?1:-1); 
      int id2 = (*((int*)e2->getData("id"))) * ((f->getUse(e2) == 1)?1:-1); 
      int id3 = (*((int*)e3->getData("id"))) * ((f->getUse(e3) == 1)?1:-1); 
      fprintf(out,"%d %d 3 %d %d %d 0\n",f->getClassification()->getId(),getClassificationNb(f),id1,id2,id3);
      int *id = new int(counter);
      f->attachData("id",id);
      counter ++;
    }
  fclose (out);
}
*/

void printGmshFace (FILE *out, Face *f, int &k)
{
  /*int typ = (f->size(0) == 3)?2:3;
  if(typ == 2)
    fprintf(out,"%d %d %d %d %d %d %d %d\n",
	    k++,typ,f->getClassification()->getId(),
	    f->getClassification()->getId(),3,
	    f->get(0,0)->getId(),
	    f->get(0,1)->getId(),
	    f->get(0,2)->getId());
  else
    fprintf(out,"%d %d %d %d %d %d %d %d %d\n",
	    k++,typ,f->getClassification()->getId(),
	    f->getClassification()->getId(),4,
	    f->get(0,0)->getId(),
	    f->get(0,1)->getId(),
	    f->get(0,2)->getId(),
	    f->get(0,3)->getId());*/
	int typ;
	if (f->size(0)==3){
		typ=2;
		fprintf(out,"%d %d %d %d %d %d %d %d\n",
	    k++,typ,f->getClassification()->getId(),
	    f->getClassification()->getId(),3,
	    f->get(0,0)->getId(),
	    f->get(0,1)->getId(),
	    f->get(0,2)->getId());
	}
	else if (f->size(0)==4){
		typ=3;
		fprintf(out,"%d %d %d %d %d %d %d %d %d\n",
	    k++,typ,f->getClassification()->getId(),
	    f->getClassification()->getId(),4,
	    f->get(0,0)->getId(),
	    f->get(0,1)->getId(),
	    f->get(0,2)->getId(),
	    f->get(0,3)->getId());
	}
	else {
		typ=6;
		fprintf(out,"%d %d %d %d %d %d %d %d %d %d\n",
	    k++,typ,f->getClassification()->getId(),
	    f->getClassification()->getId(),5,
	    f->get(0,0)->getId(),
	    f->get(0,1)->getId(),
	    f->get(0,2)->getId(),
	    f->get(0,3)->getId(),
		  f->get(0,4)->getId());
	}

}

void printGmshFace (ostream &fs, Face *f, int &k, int recur)
{
  if(recur == 2)return;
	int i;

  for(i=0;i<recur;i++) fs << " ";
  int nbrec = (f->isAdjacencyCreated(2))?f->size(2):0;
  int typ = (f->size(0) == 3)?2:3;
  if(typ == 2)
    fs << k++ << " " << typ << " " << f->getClassification()->getId() << " "
       << -nbrec << " " <<
      3 << " " <<  f->get(0,0)->getId() << " " << f->get(0,1)->getId() << " " << 
      f->get(0,2)->getId() << "\n";
  else
    fs << k++ << " " << typ << " " << f->getClassification()->getId() << " "  
       << -nbrec << " " <<
      4 << " " <<  f->get(0,0)->getId() << " " << f->get(0,1)->getId() << " " << 
      f->get(0,2)->getId() << " " << f->get(0,3)->getId() << "\n";
  for(i=0;i<nbrec;i++)printGmshFace(fs,(Face*)f->get(2,i),k,recur+1);
}

void printGmshRegion (ostream &fs, mRegion *f, int &k, int recur)
{
  if(recur == 2)return;
	int i;
  for(i=0;i<recur;i++) fs << " ";
  int nbrec = (f->isAdjacencyCreated(3))?f->size(3):0;
  switch(f->getType())
    {
    case mEntity::HEX :
      fs << 
	k++ << " " << 
	5 << " " << 
	f->getClassification()->getId() << " "<< 
	-nbrec << " " <<
	8 << " " <<  
	f->get(0,0)->getId() << " " << 
	f->get(0,1)->getId() << " " << 
	f->get(0,2)->getId() << " " << 
	f->get(0,3)->getId() << " " << 
	f->get(0,4)->getId() << " " << 
	f->get(0,5)->getId() << " " << 
	f->get(0,6)->getId() << " " << 
	f->get(0,7)->getId() << "\n";
      break;
    case mEntity::TET :
      fs << 
	k++ << " " << 
	4 << " " << 
	f->getClassification()->getId() << " "<< 
	-nbrec << " " <<
	4 << " " <<  
	f->get(0,0)->getId() << " " << 
	f->get(0,1)->getId() << " " << 
	f->get(0,2)->getId() << " " << 
	f->get(0,3)->getId() << "\n";
      break;
    }
  for(i=0;i<nbrec;i++)printGmshRegion(fs,(mRegion*)f->get(3,i),k,recur+1);
}

void printGmshEdge (ostream &fs, Edge *f, int &k, int recur)
{
  if(recur == 2)return;
  int i;
  for(i=0;i<recur;i++) fs << " ";
  int nbrec = (f->isAdjacencyCreated(1))?f->size(1):0;
  fs << k++ << " 1 "  << f->getClassification()->getId() << " "
     << -nbrec << " " <<
    2 << " " <<  f->get(0,0)->getId() << " " << f->get(0,1)->getId() << "\n";
  for(i=0;i<nbrec;i++)printGmshEdge(fs,(Edge*)f->get(1,i),k,recur+1);  
}


void mImportExport::exportGmshFile(char *fName,mMesh *theMesh)
{
  FILE *out = fopen (fName,"w");
  if(!out)throw new mExceptionFileNotFound (__LINE__,__FILE__,"");
// I change NOE to NOD
  fprintf(out,"$NOD\n");
  fprintf(out,"%d\n",theMesh->size(0));
  
  mMesh::iter it;
  for( it= theMesh->begin(0);it != theMesh->end(0); ++it)
    {
      Vertex *v = (Vertex*)(*it);
      fprintf(out,"%d %12.5E %12.5E %12.5E\n",
	      v->getId(),
	      v->point()(0),
	      v->point()(1),
	      v->point()(2));
    }
  fprintf(out,"$ENDNOD\n");
  fprintf(out,"$ELM\n");
	int count=theMesh->size(2);
	for (it = theMesh->begin(1);it != theMesh->end(1); ++it)
		{ 
      Edge *e = (Edge*)(*it);
			if (e->getClassification()->getId()==40000||e->getClassification()->getId()==20000)
				count++;
	}
  //fprintf(out,"%d\n",theMesh->size(2));
	fprintf(out,"%d\n",count);
  int k = 1;
  for(it = theMesh->begin(2);it != theMesh->end(2); ++it)
    {
      Face *f = (Face*)(*it);
      printGmshFace (out, f, k);
    }
	int pk=theMesh->size(2)+1;
	for (it = theMesh->begin(1);it != theMesh->end(1); ++it)
		{ 
			
      Edge *e = (Edge*)(*it);
			if (e->getClassification()->getId()==40000){
				fprintf(out,"%d %d %d %d %d %d %d\n",
				    pk,
						1,
						40000,
						10,
						2,
						e->get(0,0)->getId(),
						e->get(0,1)->getId());
				pk++;
			}
			if (e->getClassification()->getId()==20000){
				fprintf(out,"%d %d %d %d %d %d %d\n",
				    pk,
						1,
						20000,
						10,
						2,
						e->get(0,0)->getId(),
						e->get(0,1)->getId());
				pk++;
			}
    }
  fprintf(out,"$ENDELM\n");
  fclose (out);
}
/*
void mImportExport::exportGmshFile(ostream &f ,mMesh *theMesh) const
{
  f << "$NOE\n";
  f << theMesh->size(0) << "\n";
  for(mMesh::iter it = theMesh->begin(0);it != theMesh->end(0); ++it)
    {
      Vertex *v = (Vertex*)(*it);
      f << v->getId() << " " << v->point()(0) << " " << v->point()(1) << " " << v->point()(2) << "\n";
    }
  f << "$ENDNOE\n";
  int k = 0;

  int lev = 3;
  if(!theMesh->size(3))lev = 2;
  if(!theMesh->size(2) && !theMesh->size(3))lev = 1;

  for(mMesh::iter it = theMesh->begin(1);it != theMesh->end(1); ++it)
    {
      Edge *edge = (Edge*)(*it);
      if(lev==1 || edge->getClassification()->getLevel() == 1)
	k++;
    }
  for(mMesh::iterall it = theMesh->beginall(2);it != theMesh->endall(2); ++it)
    {
      Face *face = (Face*)(*it);
      if(lev == 2 || face->getClassification()->getLevel() == 2)
	k++;
    }

  f << "$ELM\n";
  f << k+theMesh->size(3) << "\n";

  k = 1;
  for(mMesh::iter it = theMesh->begin(1);it != theMesh->end(1); ++it)
    {
      Edge *edge = (Edge*)(*it);
      if(lev == 1 || edge->getClassification()->getLevel() == 1)
	printGmshEdge (f,edge, k,0);
    }
  for(mMesh::iterall it = theMesh->beginall(2);it != theMesh->endall(2); ++it)
    {
      Face *face = (Face*)(*it);
      if(lev == 2 || face->getClassification()->getLevel() == 2)
	printGmshFace (f,face, k,1);
    }
  for(mMesh::iter it = theMesh->begin(3);it != theMesh->end(3); ++it)
    {
      printGmshRegion (f,(mRegion*)(*it), k,0);
    }
  f << "$ENDELM\n";
}
*/

void mImportExport::exportDGFile(ostream &f ,mMesh *theMesh)
{
  f << "$NOE\n";
  f << theMesh->size(0) << "\n";
  mMesh::iter it;
  mMesh::iter end0 = theMesh->end(0);
  f.precision(16);
  for( it= theMesh->begin(0);it != end0; ++it)
    {
      const Vertex *v = (Vertex*)(*it);
      f << v->getId() << " " << v->point()(0) << " " << v->point()(1) << " " << v->point()(2) << "\n";
    }
  f << "$ENDNOE\n";

  f << "$ELM\n";
  f << theMesh->size(1) + theMesh->size(2) + theMesh->size(3) << "\n";

  int k = 1;
  mMesh::iter end1 = theMesh->end(1);
  for(it = theMesh->begin(1);it != end1; ++it)
    {
      Edge *edge = (Edge*)(*it);
      printGmshEdge (f,edge, k,0);
    }
  mMesh::iter end2 = theMesh->end(2);
  for(it = theMesh->begin(2);it != end2; ++it)
    {
      Face *face = (Face*)(*it);
      printGmshFace (f,face, k,0);
    }
  f << "$ENDELM\n";
}

mEntity * readGmshELM(istream &f, mMesh *theMesh, int recur)
{
  int i,iNbNod,iTyp,iGrp,iElm,iNbSub,nod[100];
  f >> iElm >> iTyp >> iGrp >> iNbSub >> iNbNod;
  for(i=0;i<iNbNod;i++)f >> nod[i];
  mEntity *theEntity = 0;

  switch(iTyp)
    {
    case 4 : 
      theEntity = theMesh->createTetWithVertices(nod[0],nod[1],nod[2],nod[3],theMesh->getGEntity(iGrp,3));
      break;
    case 5:
      theEntity = (mEntity*)theMesh->createHexWithVertices(nod[0],nod[1],nod[2],nod[3],nod[4],nod[5],nod[6],nod[7],
						 theMesh->getGEntity(iGrp,3));
      break;
    case 2 :
      theEntity = theMesh->createFaceWithVertices(nod[0],nod[1],nod[2],theMesh->getGEntity(iGrp,2));
      break;
    case 3 : 
      theMesh->createFaceWithVertices(nod[0],nod[1],nod[2],nod[3],theMesh->getGEntity(iGrp,2));
      break;
    case 1 : 
      theMesh->createEdge(nod[0],nod[1],theMesh->getGEntity(iGrp,1));
      break;
    case 15 : 
      Vertex *v = theMesh->getVertex(nod[0]); 
      if(iGrp > 0)
	{
	  v->classify(theMesh->getGEntity(iGrp,0));
	}		    
      break;
    }
  if(recur == 1 || iNbSub >= 0)return theEntity;
  for(i=0;i<-iNbSub;i++)
    {
      mEntity *sub = readGmshELM(f,theMesh,recur+1);
      theEntity->add(sub);
    }
  return theEntity;
}


void mImportExport::importDGFile(istream &f ,mMesh *theMesh)
{
  int NbNod,NbElm;
  char line[256];

  while(1)
    {
      if(!strcmp(line,"$NOE") ||!strcmp(line,"$NOD"))
	{
	  int NbNod;
	  f >> NbNod;
	  //	  printf("%d nodes\n",NbNod);
	  for(int i=0;i<NbNod;i++)
	    {
	      int iNod;
	      double x,y,z;
	      f >> iNod >> x >> y >> z;
	      theMesh->createVertex(iNod,x,y,z,0);
	    }
	}
      if(!strcmp(line,"$ELM"))
	{
	  int NbElm;
	  f >> NbElm; 
	  //	  printf("%d elms\n",NbElm);
	  for(int i=0;i<NbElm;i++)
	    {
	      readGmshELM(f,theMesh,0);
	    }
	}
    }
  //classifyUnclassifiedVerices (theMesh);
}





function [A,E,C,D]=createmesh(Nx,Ny,flname)
shapenum=1;
dom=[-1,1,-1,1]*10;
%dom=[-1,1,-1,1]*10+[-0.1,-0.1,0.005,0.005];
%dom=[0,1.2,-0.2,0.2]+[-0.1,-0.1,0.025,0.025];
dx=(dom(2)-dom(1))/Nx;
dy=(dom(4)-dom(3))/Ny;
[A,B]=coarsegen(Nx,Ny,dom);
% [C,D]=coarsegen(Nx,Ny,dom);
  [A,B]=refgen(A,B,2,[-1,1,-1,1]*5);
  [A,B]=refgen(A,B,2,[-1,1,-1,1]*2.5);
 % [A,B]=refgen(A,B,2,[-1,1,-1,1]*1.9);
  [A,B]=refgen(A,B,2,[-1,1,-1,1]*1.25);
% [A,B]=refgen(C,D,2,[-1,3,-1,1]*0.4*5+[-0.1,-0.1,0.005,0.005]);
% [A,B]=refgen(A,B,2,[-1,4,-1,1]*0.2*5+[-0.1,-0.1,0.005,0.005]);
% [A,B]=refgen(A,B,2,[-1,4,-1,1]*0.1*4+[-0.1,-0.1,0.005,0.005]);
% [A,B]=refgen(A,B,2,[0,1.2,-0.2,0.2]+[-0.1,-0.1,0.005,0.005]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A stores vertices, B stores faces
E=[];
C=[];
D=[];
F=[];
indA=size(A,1)+1;
for vv=1:size(B,1)
    BDindices=[];
    [inside_num,Loop, Loop2,boundary,boundary2]=intersection2(B(vv,1),B(vv,2),B(vv,3),B(vv,4),shapenum);

switch size(Loop,1)
            
            case 0
            case 1
            case 3
                indices=0;
                for ii=1:size(Loop,1)
                    indices(ii)=findindex(Loop(ii,1),Loop(ii,2),A);
                    if indices(ii)==0
                        A=[A;indA,Loop(ii,:)];
                        indA=indA+1;
                        indices(ii)=indA-1;
                    end
                end
               E=[E;indices,0,0];
               if (~isempty(Loop2))
                   indices=0;
                   for ii=1:size(Loop2,1)
                    indices(ii)=findindex(Loop2(ii,1),Loop2(ii,2),A);
                    if indices(ii)==0
                        A=[A;indA,Loop2(ii,:)];
                        indA=indA+1;
                        indices(ii)=indA-1;
                    end
                   end
%                   E=[E;indices,0,0];
              indicessize=size(indices,2);
              Esize=size(E,1);
              E(Esize+1,1:indicessize)=indices;
               end
            case 5
                indices=0;
                jj=1;
                for ii=1:5
                    indices(jj)=findindex(Loop(ii,1),Loop(ii,2),A);
                    if indices(jj)==0
                        A=[A;indA,Loop(ii,:)];
                        indA=indA+1;
                   
                    indices(jj)=indA-1; 
                    end
                    jj=jj+1;
                end
                E=[E;indices];
                
       
          case 4
                indices=0;
                for ii=1:size(Loop,1)
                    indices(ii)=findindex(Loop(ii,1),Loop(ii,2),A);
                    if indices(ii)==0
                        A=[A;indA,Loop(ii,:)];
                        indA=indA+1;
                        indices(ii)=indA-1;
                    end
                end
                E=[E;indices,0];
                if ~isempty(Loop2)
                    indices=0;
                    for ii=1:size(Loop2,1)
                    indices(ii)=findindex(Loop2(ii,1),Loop2(ii,2),A);
                    if indices(ii)==0
                        A=[A;indA,Loop2(ii,:)];
                        indA=indA+1;
                        indices(ii)=indA-1;
                    end
                    end
                    %E=[E;indices,0];
                    indicessize=size(indices,2);
              Esize=size(E,1);
              E(Esize+1,1:indicessize)=indices;
                end
            
    otherwise
            indices=0;
                for ii=1:size(Loop,1)
                    indices(ii)=findindex(Loop(ii,1),Loop(ii,2),A);
                    if indices(ii)==0
                        A=[A;indA,Loop(ii,:)];
                        indA=indA+1;
                        indices(ii)=indA-1;
                    end
                end
                Fsize=size(F,1);%for elements with many edges
                Lindices=size(indices,2);
               F(Fsize+1,1:Lindices)=indices;
        end
        if isempty(boundary)
        else
            if size(boundary, 1)==1
                fprintf('%f,%f', boundary(1,1),boundary(1,2));
            end
            if norm(boundary(1,:)-boundary(2,:),inf)>10^(-10)
            for ii=1:size(boundary,1)
                BDindices(ii)=findindex(boundary(ii,1),boundary(ii,2),A);
            end
            C=[C;BDindices];
            end
        end
        if isempty(boundary2)
        else
             if size(boundary2, 1)==1
                fprintf('%f,%f', boundary2(1,1),boundary2(1,2));
            end
            if norm(boundary2(1,:)-boundary2(2,:),inf)>10^(-10)
            for ii=1:size(boundary2,1)
                BDindices(ii)=findindex(boundary2(ii,1),boundary2(ii,2),A);
            end
            C=[C;BDindices];
            end
        end
    end
Esize=size(E,1);
Fsize=size(F);
E(Esize+1:Esize+Fsize(1),1:Fsize(2))=F;


for ii=1:4
    switch ii
        case 1
            for jj=1:Nx
                x1=dom(1)+(jj-1)*dx;
                y1=dom(3);
                x2=x1+dx;
                y2=y1;
                D=[D;findindex(x1,y1,A),findindex(x2,y2,A)];
            end
        case 2
            for jj=1:Ny
                x1=dom(2);
                y1=dom(3)+(jj-1)*dy;
                x2=x1;
                y2=y1+dy;
                D=[D;findindex(x1,y1,A),findindex(x2,y2,A)];
            end
        case 3
            for jj=1:Nx
                x1=dom(1)+(jj-1)*dx;
                y1=dom(4);
                x2=x1+dx;
                y2=y1;
                D=[D;findindex(x1,y1,A),findindex(x2,y2,A)];
            end
        case 4
            for jj=1:Ny
                x1=dom(1);
                y1=dom(3)+(jj-1)*dy;
                x2=x1;
                y2=y1+dy;
                D=[D;findindex(x1,y1,A),findindex(x2,y2,A)];
            end
    end
end
writetofile(A,E,C,D,flname);
end
function state=writetofile(A,B,C,D,flname)
state=0;
fid=fopen(flname,'wt');
%fprintf(fid, '$MeshFormat\n');
%fprintf(fid,'%d %d %d\n',[2 0 8]);
%fprintf(fid,'$EndMeshFormat\n$Nodes\n');
fprintf(fid,'$NOD\n');
fprintf(fid,'%d\n',size(A,1));
for ii=1:size(A,1)
    fprintf(fid,'%d %22.15E %22.15E %d\n',[A(ii,:),0]);
end
%fprintf(fid,'$EndNodes\n$Elements\n');
fprintf(fid,'$ENDNOD\n$ELM\n');
%fprintf(fid,'%d\n',size(B,1)+size(C,1)+size(D,1));
fprintf(fid,'%d\n',size(B,1)+size(C,1)+size(D,1));
%fprintf(fid,'%d\n',size(B,1));

%n1=size(D,1);
n1=0;
for ii=1:size(B,1)
    if (B(ii,4)==0)
        fprintf(fid,'%d %d %d %d %d %d %d %d\n',[ii+n1, 2, 601, 600, 3, B(ii,1:3)]);
    else if B(ii,5)==0
        fprintf(fid,'%d %d %d %d %d %d %d %d %d\n',[ii+n1, 3, 601, 600, 4, B(ii,1:4)]);
      %   else if B(ii,6)==0
        else
             fprintf(fid,'%d %d %d %d %d %d %d %d %d %d\n',[ii+n1, 6, 601, 600, 5, B(ii,1:5)]);
%             else if B(ii,7)==0
%                     fprintf(fid,'%d %d %d %d %d %d %d %d %d %d %d\n',[ii+n1, 7, 601, 600, 6, B(ii,1:6)]);
%                  else 
%                         fprintf(fid,'%d %d %d %d %d %d %d %d %d %d %d %d\n',[ii+n1, 8, 601, 600, 7, B(ii,1:7)]);
%                     
%                 end
 %           end
         end
    end
end

n1=n1+size(B,1);
for ii=1:size(C,1)
    fprintf(fid,'%d %d %d %d %d %d %d\n',[ii+n1, 1, 20000, 10, 2, C(ii,:)]);
end
n1=n1+size(D,1);
for ii=1:size(D,1)
    fprintf(fid,'%d %d %d %d %d %d %d\n',[ii+n1, 1, 30000, 1, 2, D(ii,:)]);
end



%fprintf(fid,'$EndElements');
fprintf(fid,'$ENDELM\n');
fclose(fid);
       
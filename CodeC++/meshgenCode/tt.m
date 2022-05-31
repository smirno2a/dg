function [a b]=tt(finame)
fid=fopen(finame);
a=[];
b=[];
st=0;
ind=0;
while ~feof(fid)
  tline=fgetl(fid);
  if tline(1)=='$'
    if strcmp(tline(1:4),'$NOD')
      st=1;
    else if strcmp(tline(1:4),'$ELM')
          st=2;
        end
    end
  else
      if st==1
          temp=sscanf(tline,'%d %f %f %f');
          if length(temp)>1
               a=[a;temp(1:3)']; 
          end
      else
%          temp=sscanf(tline,'%d %d %d %d %d %d %d %d %d %d');
          
%           if length(temp)==10
%               b=[b;temp(6:10)'];
%           else if length(temp)==9
%                   b=[b;temp(6:9)',0];
%               else if length(temp)==8
%                       b=[b;temp(6:8)',0,0];
%                   end
%               end
%           end
           temp=sscanf(tline,'%d');
           %if(size(temp,1)>1&&temp(3)==40000)
           if(size(temp,1)>1)
             b(ind+1,1:size(temp,1)-5)=temp(6:end)';
             ind=ind+1;
           end
      end
  end
end
fclose(fid);
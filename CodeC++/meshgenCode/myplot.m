function myplot(a,b)
figure;
plot(0,0);
hold on;
% for ii=1:size(b,1)
%    
%    if b(ii,4)==0
%        for jj=1:3
%            hold on;
%            plot([a(b(ii,jj),2),a(b(ii,(jj+1-(jj==3)*3)),2)],[a(b(ii,jj),3),a(b(ii,(jj+1-(jj==3)*3)),3)]); 
%        end
%    else if b(ii,5)==0
%            for jj=1:4
%                hold on;
%                plot([a(b(ii,jj),2),a(b(ii,(jj+1-(jj==4)*4)),2)],[a(b(ii,jj),3),a(b(ii,(jj+1-(jj==4)*4)),3)]); 
%            end
%        else
%            for jj=1:5
%                hold on;
%                plot([a(b(ii,jj),2),a(b(ii,(jj+1-(jj==5)*5)),2)],[a(b(ii,jj),3),a(b(ii,(jj+1-(jj==5)*5)),3)]); 
%            end
%        end
%    end
% end
for ii=1:size(b,1)
  for jj=1:size(b,2)
      if jj<size(b,2)
         if b(ii,jj)==0
         else if b(ii,jj+1)==0
                 plot([a(b(ii,jj),2),a(b(ii,1),2)],[a(b(ii,jj),3),a(b(ii,1),3)]);
             else
                 plot([a(b(ii,jj),2),a(b(ii,jj+1),2)],[a(b(ii,jj),3),a(b(ii,jj+1),3)]);
             end
         end
      else
          if b(ii,jj)==0
          else
              plot([a(b(ii,jj),2),a(b(ii,1),2)],[a(b(ii,jj),3),a(b(ii,1),3)]);
          end
      end
      
      
  end
end
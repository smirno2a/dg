function exst=findindex(x,y,A)
exst=0;
for ii=1:size(A,1)
    if norm([x,y]-A(ii,2:3),inf)<10^(-10)
        exst=A(ii,1);
        break;
    end
end
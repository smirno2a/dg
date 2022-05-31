function bool = sq_inside( x, y, dx, dy )
%SQ_INSIDE Summary of this function goes here
%  return 1 if at least one vertex is inside the domain
bool=0;    
for ii=1:4
        x1=x+(ii==2||ii==3)*dx;
        y1=y+(ii==3||ii==4)*dy;
        if is_inside(x1,y1,2)
            bool=1;
        end
   %     break;
end


function [inside_num,Loop,Loop2,boundary,boundary2]=intersection2(xl,yl,dx,dy,shapenum)
inside_num=0;
Loop2=[];
boundary2=[];
for ii=1:4
    x=xl+(ii==2||ii==3)*dx;
    y=yl+(ii==3||ii==4)*dy;
    if is_inside(x,y,shapenum)
        inside_num=inside_num+1;
    end
end
switch inside_num
    case 0
        Loop=[];
        boundary=[];
    case 4
        Loop=[];
        boundary=[];
        for ii=1:4
            x=xl+(ii==2||ii==3)*dx;
            y=yl+(ii==3||ii==4)*dy;
            Loop=[Loop;x,y];
        end
    otherwise
        Loop=[];
        boundary=[];
        for ii=1:4
            x(ii)=xl+(ii==2||ii==3)*dx;
            y(ii)=yl+(ii==3||ii==4)*dy;
        end
        for ii=1:4
            if is_inside(x(ii),y(ii),shapenum)
                Loop=[Loop;x(ii),y(ii)];
                if is_inside(x(mod(ii+1,4)+(ii==3)*4),y(mod(ii+1,4)+(ii==3)*4),shapenum)
                else
                    [xi,yi]=find_intersection(x(ii),y(ii),x(mod(ii+1,4)+(ii==3)*4),y(mod(ii+1,4)+(ii==3)*4),shapenum);
                    boundary=[boundary;xi,yi];
                    if norm([xi,yi]-[x(ii),y(ii)],inf)>10^(-10)
                    Loop=[Loop;xi,yi];
                    end
                end
                    
            else if is_inside(x(mod(ii+1,4)+(ii==3)*4),y(mod(ii+1,4)+(ii==3)*4),shapenum)
                    [xi,yi]=find_intersection(x(ii),y(ii),x(mod(ii+1,4)+(ii==3)*4),y(mod(ii+1,4)+(ii==3)*4),shapenum);
                     boundary=[boundary;xi,yi];
                    if norm([xi,yi]-[x(mod(ii+1,4)+(ii==3)*4),y(mod(ii+1,4)+(ii==3)*4)],inf)>10^(-10)
                    Loop=[Loop;xi,yi];
                   
                    end
                else
                end
                
            end
        end
    
end
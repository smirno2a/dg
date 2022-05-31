function [inside_num,Loop, Loop2, boundary, boundary2]=intersection(xl,yl,dx,dy,shapenum)
corner=[1.009, 0];
inside_num=0;
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
        Loop2=[];
        boundary=[];
        boundary2=[];
    case 2
        Loop=[];
        Loop2=[];
        boundary=[];
        boundary2=[];
        [xi,yi,wi,zi]=find_intersection(xl+dx,yl,xl+dx,yl+dy,2);
        [xi2,yi2,wi2,zi2]=find_intersection(xl,yl+dy,xl,yl,2);
        if (isempty(wi))&&(isempty(wi2))
              for ii=1:4
            x(ii)=xl+(ii==2||ii==3)*dx;
            y(ii)=yl+(ii==3||ii==4)*dy;
        end
        for ii=1:4
            if is_inside(x(ii),y(ii),shapenum)
                Loop=[Loop;x(ii),y(ii)];
                if is_inside(x(mod(ii+1,4)+(ii==3)*4),y(mod(ii+1,4)+(ii==3)*4),shapenum)
                else
                    [xi,yi,wi,zi]=find_intersection(x(ii),y(ii),x(mod(ii+1,4)+(ii==3)*4),y(mod(ii+1,4)+(ii==3)*4),shapenum);
                    boundary=[boundary;xi,yi];
                    if norm([xi,yi]-[x(ii),y(ii)],inf)>10^(-10)
                    Loop=[Loop;xi,yi];                    
                    end
                end
                    
            else if is_inside(x(mod(ii+1,4)+(ii==3)*4),y(mod(ii+1,4)+(ii==3)*4),shapenum)
                    [xi,yi,wi,zi]=find_intersection(x(ii),y(ii),x(mod(ii+1,4)+(ii==3)*4),y(mod(ii+1,4)+(ii==3)*4),shapenum);
                     boundary=[boundary;xi,yi];
                    if norm([xi,yi]-[x(mod(ii+1,4)+(ii==3)*4),y(mod(ii+1,4)+(ii==3)*4)],inf)>10^(-10)
                    Loop=[Loop;xi,yi];                   
                    end

                end
                
            end
        end
        else
          if is_inside(xl,yl,shapenum)
            Loop=[xl,yl];
            [xi,yi,wi,zi]=find_intersection(xl,yl,xl+dx,yl,2);
            Loop=[Loop;xi,yi];
            boundary=[xi,yi];
            [xi,yi,wi,zi]=find_intersection(xl+dx,yl+dy,xl,yl+dy,2);
            Loop2=[xi,yi;xl,yl+dy];
            boundary2=[xi,yi];
            [xi,yi,wi,zi]=find_intersection(xl,yl+dy,xl,yl,2);
            Loop2=[Loop2;xi,yi];
            boundary2=[boundary2;xi,yi];
            Loop=[Loop;wi,zi];
            boundary=[boundary;wi,zi];
          else
            [xi,yi,wi,zi]=find_intersection(xl,yl,xl+dx,yl,2);
            Loop=[xi,yi;xl+dx,yl];
            boundary=[xi,yi];
            [xi,yi,wi,zi]=find_intersection(xl+dx,yl,xl+dx,yl+dy,2);
            Loop=[Loop;xi,yi];
            boundary=[boundary;xi,yi];
            Loop2=[wi,zi;xl+dx,yl+dy];
            boundary2=[wi,zi];
            [xi,yi,wi,zi]=find_intersection(xl+dx,yl+dy,xl,yl+dy,2);
            Loop2=[Loop2;xi,yi];
            boundary2=[boundary2;xi,yi];
          end
        end
    case 4
        Loop=[];
        Loop2=[];
        boundary=[];
        boundary2=[];
        if cornerinside(xl,yl,dx,dy)
            [x1,y1,x2,y2]=find_spec_interception(xl);
            for ii=1:4
                x=xl+(ii==2||ii==3)*dx;
                y=yl+(ii==3||ii==4)*dy;
                Loop=[Loop;x,y];
            end
            Loop=[Loop;x1,y1;corner;x2,y2];
            boundary=[x1,y1;corner];
            boundary2=[corner;x2,y2];
        else
            Loop=[xl,yl];
            Loop=[Loop;xl+dx,yl;];
            [xi,yi,wi,zi]=find_intersection(xl+dx,yl,xl+dx,yl+dy,2);
            [xi2,yi2,wi2,zi2]=find_intersection(xl,yl+dy,xl,yl,2);
            if isempty(wi)||isempty(wi2)
                Loop=[Loop;xl+dx,yl+dy;xl,yl+dy];
            else
                Loop=[Loop;xi,yi;];
                Loop2=[wi,zi;xl+dx,yl+dy;xl,yl+dy];
                [xi2,yi2,wi2,zi2]=find_intersection(xl,yl+dy,xl,yl,2);
                Loop2=[Loop2;xi2,yi2];
                Loop=[Loop;wi2,zi2;];
                boundary=[xi,yi;wi2,zi2];
                boundary2=[xi2,yi2;wi,zi];
            
                
            end
%             for ii=2:1:4
%                 x=xl+(ii==2||ii==3)*dx;
%                 y=yl+(ii==3||ii==4)*dy;
%                 
%                 Loop=[Loop;x,y];
%             end
        end

       
       
        
            
                
           
    otherwise
        Loop=[];
        Loop2=[];
        boundary=[];
        boundary2=[];
        check=0;
        [xi,yi,wi,zi]=find_intersection(xl+dx,yl,xl+dx,yl+dy,2);
        if (~isempty(wi))
            Loop=[xl,yl];
            Loop=[Loop;xl+dx,yl;];
            [xi,yi,wi,zi]=find_intersection(xl+dx,yl,xl+dx,yl+dy,2);
            [xi2,yi2,wi2,zi2]=find_intersection(xl+dx,yl+dy,xl,yl+dy,2);
            [xi3,yi3,wi3,zi3]=find_intersection(xl,yl+dy,xl,yl,2);
           
                Loop=[Loop;xi,yi;];
                Loop2=[wi,zi;xl+dx,yl+dy;xi2,yi2];
             %   [xi2,yi2,wi2,zi2]=find_intersection(xl,yl+dy,xl,yl,2);
              %  Loop2=[Loop2;xi2,yi2];
                Loop=[Loop;xi3,yi3;];
                boundary=[xi,yi;xi3,yi3];
                boundary2=[xi2,yi2;wi,zi];
        else
        for ii=1:4
            x(ii)=xl+(ii==2||ii==3)*dx;
            y(ii)=yl+(ii==3||ii==4)*dy;
        end
        for ii=1:4
            if is_inside(x(ii),y(ii),shapenum)
                Loop=[Loop;x(ii),y(ii)];
                if is_inside(x(mod(ii+1,4)+(ii==3)*4),y(mod(ii+1,4)+(ii==3)*4),shapenum)
                else
                    [xi,yi,wi,zi]=find_intersection(x(ii),y(ii),x(mod(ii+1,4)+(ii==3)*4),y(mod(ii+1,4)+(ii==3)*4),shapenum);
                    boundary=[boundary;xi,yi];
                    if norm([xi,yi]-[x(ii),y(ii)],inf)>10^(-10)
                    Loop=[Loop;xi,yi];                    
                    end
%                     if check==0
%                         Loop=[Loop;corner];
%                         boundary=[boundary;corner];
%                         check=1;
%                     end
                end
                    
            else if is_inside(x(mod(ii+1,4)+(ii==3)*4),y(mod(ii+1,4)+(ii==3)*4),shapenum)
                    [xi,yi,wi,zi]=find_intersection(x(ii),y(ii),x(mod(ii+1,4)+(ii==3)*4),y(mod(ii+1,4)+(ii==3)*4),shapenum);
                     boundary=[boundary;xi,yi];
                    if norm([xi,yi]-[x(mod(ii+1,4)+(ii==3)*4),y(mod(ii+1,4)+(ii==3)*4)],inf)>10^(-10)
                    Loop=[Loop;xi,yi];                   
                    end
%                     if check==0
%                         Loop=[Loop;corner];
%                         boundary=[boundary;corner];
%                         check=1;
%                     end
                else
                end
                
            end
        end
        end
    
end
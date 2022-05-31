function [x,y,w,z]=find_intersection(x1,y1,x2,y2,shapenum)
R=1.0;
epsl=10^(-10);
x=[];
y=[];
w=[];
z=[];
switch shapenum
    case 1
        if abs(x1-x2)<10^(-10)
            x=x1;
            y=sqrt(R^2-x1^2);
            if (y>(y1+10^(-10))&&y>(y2+10^(-10)))
            y=-y;
            end
        else
            y=y1;
            x=sqrt(R^2-y1^2);
                if (x>(x1+10^(-10))&&x>(x2+10^(-10)))
                x=-x;
                end
        end
    case 2
        if abs(x1-x2)<10^(-10)
            x=x1;
            if (x<=1.009+epsl&&x>=-epsl)
                y=airf(x1);
                if (y>(min(y1,y2)-10^(-10))&&y<(max(y1,y2)+10^(-10)))
                    if(-y>(min(y1,y2)-10^(-10))&&-y<(max(y1,y2)+10^(-10)))
                        w=x1;
                        z=-y;
                        if (y1<y2)
                            z=-z;
                            y=-y;
                        end
                        if(abs(y-z)<1e-8)
                            x=0;
                            y=0;
                            z=0;
                            w=0;
                        end
                    end
                else if (-y>(min(y1,y2)-10^(-10))&&-y<(max(y1,y2)+10^(-10)))
                        y=-y;
                    else
                        x=[];
                        y=[];
                    end
                    
                end
            end
        else
            y=y1;
            if y1>0
%                 if x1<0
%                     x=fzero(@(x)(airf(x)-y1),0);
%                 else if x2<0
%                         x=fzero(@(x)(airf(x)-y1),0);
%                     else
                        x=fzero(@(x)(airf(x)-y1),x1);
%                     end
                    
%                 end
                  if (x>min(x1,x2)&&x<max(x1,x2))
                  else
                      x=[];
                      y=[];
                  end
            else if y1<0
%                     if x1<0
%                         x=fzero(@(x)(-airf(x)-y1),0);
%                     else if x2<0
%                             x=fzero(@(x)(-airf(x)-y1),0);
%                         else
                            x=fzero(@(x)(-airf(x)-y1),x1);
%                         end
%                     end
%                       
                    if (x>min(x1,x2)&&x<max(x1,x2))
                    else
                      x=[];
                      y=[];
                    end
                 else x=0;
                end
            end
        end
    otherwise
end
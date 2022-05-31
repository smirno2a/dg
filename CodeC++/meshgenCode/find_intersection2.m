function [x,y]=find_intersection2(x1,y1,x2,y2,shapenum)
R=1.0;
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
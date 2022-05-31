function bool=is_inside(x,y,shapenum)
R=1.0;
epsl=10^(-10);
bool=0;
switch shapenum
    case 1
        if(x^2+y^2>=R^2-epsl)
       % if (x^2<0.01-10^(-15)&&y^2<0.01-10^(-12))      
            bool=1;
        end
    case 2
        if(x<epsl||x>1.009-epsl)
             bool=1;
        else
            if (abs(y)>airf(x)-epsl)
                bool=1;
            end
        end
    otherwise
end
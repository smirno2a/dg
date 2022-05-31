function bol = needref( x, y, dx, dy, rang)
%NEEDREF Summary of this function goes here
%   Detailed explanation goes here
bol=0;    
for ii=1:4
        x1=x+(ii==2||ii==3)*dx;
        y1=y+(ii==3||ii==4)*dy;
        if x1>rang(1)-10^-10 && x1<rang(2)+10^-10 && y1>rang(3)-10^-10 && y1<rang(4)+10^-10
            bol=bol+0.25;
        end
end
if bol<1-10^-3
    bol=0;
else bol=1;
end

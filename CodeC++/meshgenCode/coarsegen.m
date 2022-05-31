function [ A, B ] = coarsegen( Nx, Ny, dom)
%COARSEGEN Summary of this function goes here
%   dom=[xl,xr,yl,yr]
shapenum=2;
dx=(dom(2)-dom(1))/Nx;
dy=(dom(4)-dom(3))/Ny;
A=[];
B=[];
indA=1;
for m=1:(Ny+1)
    for n=1:(Nx+1)
        x=dom(1)+(n-1)*dx;
        y=dom(3)+(m-1)*dy;
        if is_inside(x,y,shapenum)
            A=[A;indA, x, y;];%index and coordinates
            indA=indA+1;
        end
    end
end

for m=1:Ny
    for n=1:Nx
        x=dom(1)+(n-1)*dx;
        y=dom(3)+(m-1)*dy;
        if sq_inside(x,y,dx,dy)
            B=[B; x, y, dx, dy];%left bottom corner and size
        end
    end
end

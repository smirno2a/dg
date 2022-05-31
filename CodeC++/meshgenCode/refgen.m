function [ A, B] = refgen( C, D, refnum, rang )
%REFGEN Summary of this function goes here
%   Detailed explanation goes here
shapenum=2;
length=size(D,1);
ind=1;
A=C;
B=D;
while ind<=length
    if needref(B(ind,1),B(ind,2),B(ind,3),B(ind,4),rang)
        dxr=B(ind,3)/refnum;
        dyr=B(ind,4)/refnum;%sizes for small elements
        for m=1:(refnum-1)
            for n=1:(refnum+1)
                xr=B(ind,1)+m*dxr;
                yr=B(ind,2)+(n-1)*dyr;
                if is_inside(xr, yr, shapenum)
                    indA=size(A,1);
                    A=[A;indA+1, xr, yr];
                end
            end
            xr=B(ind,1);
            yr=B(ind,2)+m*dyr;
            if is_inside(xr, yr, shapenum)
                    indA=size(A,1);
                    A=[A;indA+1, xr, yr];
            end
            xr=B(ind,1)+B(ind,3);
            if is_inside(xr, yr, shapenum)
                    indA=size(A,1);
                    A=[A;indA+1, xr, yr];
            end
        end
        Btemp=B(ind,:);
        B(ind,:)=[];
        ind=ind-1;
        length=length-1;
        for m=1:refnum
            for n=1:refnum
                xr=Btemp(1)+(m-1)*dxr;
                yr=Btemp(2)+(n-1)*dyr;
                if sq_inside(xr, yr, dxr, dyr);
                    B=[B; xr, yr, dxr, dyr];
                end
            end
        end
    end
    ind=ind+1;
end

end


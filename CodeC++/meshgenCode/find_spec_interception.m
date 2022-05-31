function [x1,y1,x2,y2]=find_spec_interception(xl)
x1=xl;
x2=xl;
y1=airf(xl);
y2=-y1;

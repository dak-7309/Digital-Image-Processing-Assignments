%input matrix
matrix=[10,10,13,17;11,12,7,11;10,13,5,19;23,17,9,8]

%input coordinates at which pixel value is to be interpolated
x=input('Enter Q- x coordinate ')
y=input('Enter Q- y coordinate ')

%finding 4 neighbours to (x,y), handling corner cases
if ceil(x)==x
x1=x-1;
x2=x+1;
else
x1=floor(x);
x2=ceil(x);
end

if ceil(y)==y
y1=y-1;
y2=y+1;
else
y1=floor(y);
y2=ceil(y);
end

%x matrix for bilinear interpolation
X_matrix=[x1,y1,x1*y1,1; x2,y1,x2*y1,1; x2,y2,x2*y2,1; x1,y2,x1*y2,1]
SS=size(matrix);
rows=SS(1);
columns=SS(2);

%finding pixel values at the 4 neighbours to (x,y), handling corner cases
if (y1+1)<1 || (y1+1)>rows || (x1+1)<1 || (x1+1)>columns
pixel1=0;
else
pixel1=matrix(y1+1,x1+1);
end

if (y1+1)<1 || (y1+1)>rows || (x2+1)<1 || (x2+1)>columns
pixel2=0;
else
pixel2=matrix(y1+1,x2+1);
end

if (y2+1)<1 || (y2+1)>rows || (x2+1)<1 || (x2+1)>columns
pixel3=0;
else
pixel3=matrix(y2+1,x2+1);
end

if (y2+1)<1 || (y2+1)>rows || (x1+1)<1 || (x1+1)>columns
pixel4=0;
else
pixel4=matrix(y2+1,x1+1);
end


pixel_matrix=[pixel1;pixel2;pixel3;pixel4]
%using matrix operation to find out a from xa=p, a=inv(x)*p
coeff=(X_matrix)\pixel_matrix

%now we have a,b,c,d, we can find out pixel values at (x,y) by multiplying
%coefficients to x matrix for (x,y)
Q_matrix=[x,y,x*y,1]

%pixel value obtained
Q_value=Q_matrix*coeff
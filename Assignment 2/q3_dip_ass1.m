%point selection tool offered by Matlab, helps us identify 3 non colinear
%points, and save these points as 

fixedPoints=[177.833333333333,140.250000000000;183.750000000000,235.750000000000;140.916666666667,221.083333333333];
movingPoints=[153.333333333333,112.583333333333;176.750000000000,204.583333333333;131.750000000000,199.750000000000];


%uncomment the below line for manual input
%cpselect('cameraman.tif', 'referenceIm.jpg');


im11=imread('referenceIm.jpg');
figure
imshow(im11)
%displays refernce image
h=gca;
h.Visible='On';

z=ones(3,1);

%fixedPoints and movingPoints are 3 set of points obtained from the
%aforementioned tool
reference=[fixedPoints z];
unregistered=[movingPoints z];

%obtaining the transformation matrix to move from refernce image to
%unregistered image, then taking its inverse
T_dash=pinv(reference)*unregistered;
T=inv(T_dash);


%now we have T, we need to apply this on our unregistered image using
%inverse mapping resuing code from q2, which will give us the registered
%image and we can compare this to the reference image

image=imread('cameraman.tif');
image=imcrop(image,[1 1 255 255]);

figure
%unregistered image
imshow(image)
h=gca;
h.Visible='On';

[row,col]=size(image);
%reusing code from q2 for inverse mapping
[X,Y]=meshgrid(-2*col:2*col,-2*row:2*row);
Z=ones(4*col+1,4*row+1);
sourceCoor=[X(:) Y(:) Z(:)]*inv(T);
Xs=(sourceCoor(:,1));
Ys=(sourceCoor(:,2));

[X,Y]=meshgrid(0:col-1,0:row-1);
%Vq=griddata(X(:),Y(:),double(image(:)),Xs(:),Ys(:));

X_vector=X(:);
Y_vector=Y(:);

OUTPUT=ones(1050625,1);

for i=1:size(Xs,1)
    x=Xs(i);
    y=Ys(i);
    
    %reusing code from q1 to do bilinear interpolation
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



    X_matrix=[x1,y1,x1*y1,1; x2,y1,x2*y1,1; x2,y2,x2*y2,1; x1,y2,x1*y2,1];

    if (y1+1)<1 || (y1+1)>255 || (x1+1)<1 || (x1+1)>255
        pixel1=double(0);
    else
        pixel1=double(image(y1+1,x1+1));
    end

    if (y1+1)<1 || (y1+1)>255 || (x2+1)<1 || (x2+1)>255
        pixel2=double(0);
    else
        pixel2=double(image(y1+1,x2+1));
    end

    if (y2+1)<1 || (y2+1)>255 || (x2+1)<1 || (x2+1)>255
        pixel3=double(0);
    else
        pixel3=double(image(y2+1,x2+1));
    end

    if (y2+1)<1 || (y2+1)>255 || (x1+1)<1 || (x1+1)>255
        pixel4=double(0);
    else
        pixel4=double(image(y2+1,x1+1));
    end

    pixel_matrix=[pixel1;pixel2;pixel3;pixel4];
    coeff=inv(X_matrix)*pixel_matrix;

    Q_matrix=[x,y,x*y,1];
    OUTPUT(i)=Q_matrix*coeff;
 
    
end


figure
V=reshape(OUTPUT,4*col+1,4*row+1);
imshow(V,[])
%registered image obtained
h=gca;
h.Visible='On';










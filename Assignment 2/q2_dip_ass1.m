%as per clear instructions by Sir, I have written a code using cameraman.tif as input,
%any other input matrix can be used as well, q1 matrix can be used too but
%according to Sirs 1st clarification, Ive directly extended it to
%cameraman.tif

%given transformation matrix for scaling
T=[2 0 0;0 2 0; 0 0 2];

%reading image
image=imread('cameraman.tif');
image=imcrop(image,[1 1 255 255]);

%displaying initial image
figure
imshow(image)
h=gca;
h.Visible='On';


[row,col]=size(image);
%creating output grid of 4 times the height and width of hour grid
[X,Y]=meshgrid(-2*col:2*col,-2*row:2*row);
Z=ones(4*col+1,4*row+1);
sourceCoor=[X(:) Y(:) Z(:)]*inv(T);
%using inverse mapping to find coordinate points to be interpolated in the
%input image
Xs=(sourceCoor(:,1));
Ys=(sourceCoor(:,2));

[X,Y]=meshgrid(0:col-1,0:row-1);

%using the loop, we iterate over every pixel in output scaled image and
%assign it a pixel value using the pixel value at its corresponding point
%to be interpolated
X_vector=X(:);
Y_vector=Y(:);

OUTPUT=ones(1050625,1);

for i=1:size(Xs,1)
    x=Xs(i);
    y=Ys(i);
    
    %here we reuse the code from q1 for bilinear interpolation
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
    
    %assign the output pixel with the pixel value found at the interpolated
    %point
    OUTPUT(i)=Q_matrix*coeff;
 
    
end


figure
%reshaping our vector into a matrix of size of meshgrid
V=reshape(OUTPUT,4*col+1,4*row+1);

%the scaled image
imshow(V,[])
h=gca;
h.Visible='On';



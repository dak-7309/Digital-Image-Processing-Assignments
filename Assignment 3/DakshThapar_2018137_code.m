%%Question 1 Solution
disp("Solution to 1 goes here")

image1=imread('DakshThapar_2018137_cameraman.tif');
%figure
%imshow(image1)
%h=gca; 
%h.Visible='On';
%display input image


image2=imread('DakshThapar_2018137_specified image.jpg');
%figure
%imshow(image2)
%h=gca;
%h.Visible='On';
%display specified image

    

x1=double(imread('DakshThapar_2018137_cameraman.tif'));
h1=zeros(256,1);

%create histogram for input image
for i=0:255
    % for every pixel value from 0 to 255, count number of occurences for
    % the whole image
    h1(i+1)=sum(sum(x1==i));
end
r1=0:255;

%figure,stem(r1,h1);
%plot input image histogram

%finding the transfer function for input image
s=sum(h1);
Tr=zeros(256,1);
for i=0:255
    %(L-1)*CDF
    % CDF is calculated by summing over all pixel values uptil i in every
    % iteration divided by sum of total occurences of all pixels
    Tr(i+1)=round(255*(sum(h1(1:i))/s));
end
    

%_________________________________________________

x2=double(imread('DakshThapar_2018137_specified image.jpg'));
h2=zeros(256,1);
%create histogram for sepcified image
for i=0:255
    % for every pixel value from 0 to 255, count number of occurences for
    % the whole image
    h2(i+1)=sum(sum(x2==i));
end
r2=0:255;
%figure,stem(r2,h2);
%plot specified image histogram

s=sum(h2);
%finding the transfer function for input image
Gz=zeros(256,1);
for i=0:255
    %(L-1)*CDF
    % CDF is calculated by summing over all pixel values uptil i in every
    % iteration divided by sum of total occurences of all pixels
    Gz(i+1)=round(255*(sum(h2(1:i))/s));
end
 
%_________________________________________________


Match=zeros(256,1,'uint8');


for idx=1:256
    %for every value in Tr, finds the value closest to it and saves it
    %index, basically minimizes difference of values and return index
    val=abs(Tr(idx)-Gz);
    [extra,ind]=min(val);
    ppp=ind-1;
    Match(idx)=ppp;
end

%obtain the ouput image by referencing for every pixel in imsge matrix the
%mapping of our Match function
out=zeros(256,256);
for i=1:256
    for j=1:256
        out(i,j)=Match(double(image1(i,j))+1);
    end
end


%figure
%imshow(out,[])
%h=gca;
%h.Visible='On';
%display matched image

x3=double(out);
h3=zeros(256,1);
for i=0:255
    h3(i+1)=sum(sum(x3==i));
end
%display matched image histogram
r3=0:255;
%figure,stem(r3,h3);






subplot(2,3,1),imshow(image1),title('input image');
subplot(2,3,2),imshow(image2),title('specified image');
subplot(2,3,3),imshow(out,[]),title('matched image');
subplot(2,3,4),stem(r1,h1);
subplot(2,3,5),stem(r2,h2);
subplot(2,3,6),stem(r3,h3);


%% Question 3b Solution
disp('Solution to 3b goes here')

f=[6,7;8,9];
%f=[6,7,2;8,9,3;1,2,3];
H=[0,1,0;1,-4,1;0,1,0];

h=zeros(3,3);
for i=1:3
    for j=1:3
        h(i,j)=H(3-i+1,3-j+1);
    end
end

[r1,c1]=size(f);
[r2,c2]=size(h);

matrix=zeros(r1+r2*2-2,c1+c2*2-2);
%big matrix that will fit our convolution brackets (with padding)

for i=r2:r1+r2-1
    for j=c2:c1+c2-1
        %filling matrix center with our image
        matrix(i,j)=f(i-r2+1,j-c2+1);
    end
end    

%final result matrix that will have convolved values (without padding)
result=zeros(r1+r2-1,c1+c2-1);
[R2,C2]=size(result);

deltaX=r2-1;
deltaY=c2-1;

for i=1:R2
    for j=1:C2
        centerX=deltaX+i-1;
        centerY=deltaY+j-1;
        %convolution matrix operations, (element wise multiplication)
        result(i,j)=sum(sum(matrix(centerX-1:centerX+1,centerY-1:centerY+1).*h));
    end
end

disp(result)

%__________________________________________________________________________
%% Question 4 Solution
disp('Solution to 4 goes here')

f=imread("DakshThapar_2018137_chandrayaan.jpg");
H=[0,1,0;1,-4,1;0,1,0];

h=zeros(3,3);
for i=1:3
    for j=1:3
        % flip filter matrix by 180 degrees
        h(i,j)=H(3-i+1,3-j+1);
    end
end

[r1,c1]=size(f);
[r2,c2]=size(h);

%reusing code from Q3 a)
matrix=zeros(r1+r2*2-2,c1+c2*2-2); 
[R1,C1]=size(matrix);

for i=r2:r1+r2-1
    for j=c2:c1+c2-1
        matrix(i,j)=f(i-r2+1,j-c2+1);
    end
end    

result=zeros(r1+r2-1,c1+c2-1);
[R2,C2]=size(result);

deltaX=r2-1;
deltaY=c2-1;

for i=1:R2
    for j=1:C2
        centerX=deltaX+i-1;
        centerY=deltaY+j-1;
        result(i,j)=sum(sum(matrix(centerX-1:centerX+1,centerY-1:centerY+1).*h));
    end
end
[RR,CC]=size(result);
TEMP=result(2:RR-1,2:CC-1);
%cropping result matrix of same size of image

f=double(f);
g4=f-TEMP;
%image sharpening, c=-1 as center is negative, immage - grad function*c
figure
subplot(1,2,1),imshow(uint8(f)),title('original image');
subplot(1,2,2),imshow(uint8(g4)),title('sharpened image');
%__________________________________________________________________________

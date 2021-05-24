%% Question 1 Solution
disp('Solution to 1 goes here');

image=imread("DakshThapar_2018137_Noise-lines.jpg");
imshow(image)
title('Original Image')
%reading noisy image (noise is horizontal lines) 

hh=zeros(512,512);
for i=1:512
    for j=254:258
        hh(i,j)=1;
    end
end
for i=255:256
    for j=1:512
        hh(i,j)=0;
    end
end
%creating a manual bandpass filter that rejects the frequencies comprising
%of noise
    
noise=hh;
%bandreject filter to filter out the noise
H=1-hh;

S=log(1+abs((noise)));
figure,imshow(S,[])
title('Magnitude Spectrum of Noise')
%displaying magnitude spectrum of noise

S=log(1+abs((H)));
figure,imshow(S,[])
title('Magnitude Spectrum of Filter')
%displaying magnitude spectrum of the filter

figure,mesh(abs((H)))
title('Mesh Plot Magnitude Spectrum of Filter')
colormap('default')
%displaying mesh plot of the filter

%filter in spatial domain
H_time=real(ifft2(ifftshift(H)));

%standardizing spatial domain filter by subtracting min value and dividing
%by range
real_H_time=255*((H_time-min(min(H_time)))./(max(max(H_time)))-min(min(H_time)));
figure,imshow(fftshift(log(1+real_H_time)),[])
title('Spatial Domain Filter')
%displaying spatial domain filter

figure,mesh(fftshift(real_H_time))
title('Mesh Plot Spatial Domain Filter')
colormap('default')
figure
%displaying mesh plot time domain filter

IMG=fftshift(fft2(double(image),size(H,1),size(H,2)));
%applying DFT to image and shifting it to the size of the filter H
imshow(mat2gray(log(1+abs(IMG))))
title('Magnitude Spectrum of Original Image');
%displaying magnitude spectrum of original image

%applying hadamard product in fourier domain
fourier_filtered_image=(H).*IMG;
%taking real part of inverse fft (shifted)
output_image=real(ifft2(ifftshift((fourier_filtered_image))));
%cropping out portion from top left of image as same size as input image
output_image=output_image(1:size(image,1),1:size(image,2));

S=log(1+abs(fourier_filtered_image));
figure,imshow(S,[])
title('Magnitude Spectrum of Filtered Image')
%displaying magnitude spectrum of filtered image

figure,imshow(output_image,[])
title('Filtered Image via Fourier Domain Filtering Process')
%displaying filtered image via Fourier Domain Filtering Process

figure,imshow(double(image)-double(output_image),[])
title('Filtered Noise')
%subtracting output and input images to display noise

figure
hh=conv2(image,ifftshift(H_time),'same');
%performing convolution in time domain itself and converting to grayscale
gg=mat2gray(hh);
imshow(gg)
title('Filtered Image obtained via Spatial Domain Filtering Process')
%displaying filtered image obtained via Spatial Domain Filtering Process





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







%% Question 2 Solution
disp('Solution to 2 goes here');


img=double(imread("DakshThapar_2018137_Barbara.bmp"));
figure,imshow(img,[]);
title('Original image');

f=[0 1 0; 1 -4 1; 0 1 0];
%input filter

f_=f(2:3,2:3);
fil = padarray(f_, [512 512], 0, 'post');
%shifting centre and hence padding with zeros to give final size 514x514 and
%then using mod operation to put the pixels in the first row and first
%column into their updated spots

fil(mod(-1,514)+1,mod(-1,514)+1)=f(1,1);
fil(mod(0,514)+1,mod(-1,514)+1)=f(2,1);
fil(mod(1,514)+1,mod(-1,514)+1)=f(3,1);
fil(mod(-1,514)+1,mod(1,514)+1)=f(1,3);
fil(mod(-1,514)+1,mod(0,514)+1)=f(1,2);

mesh(fil);
image = padarray(img, [2 2], 0, 'post');
%padding the input image with zeros to give final size 514x514

H=fftshift(fft2(double(fil)));
IMG=fftshift(fft2(double(image)));
%applying fft to both filter and image

figure,imshow(H,[])
title('Filter in Fourier Domain');
%displaying filter in fourier domain

fourier_filtered_image = (H).*IMG;
%taking hadamard product of filter and image in fourier domain
output_image=real(ifft2(ifftshift((fourier_filtered_image)))); 
%applying inverse fft to bring it back to spatial domain
output_image=output_image(1:size(img,1), 1:size(img,2));
%cropping out the size, to be the same size as the input image

figure, imshow(output_image,[])
title('Filtered Image')
%displaying filtered image

figure,imshow(uint8(img-output_image))
%applying laplacian filter to actually sharpen our image
title('Output Sharpened Image');
%display sharpened image




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%% Question 3 Solution
disp('Solution to 3 goes here');

image=imread("DakshThapar_2018137_Barbara.bmp");
imshow(image)
title('Original Image')
[M,N]=size(image);
%input image

u=0:(M-1);
v=0:(N-1);
[U,V]=meshgrid(v,u);
%preparing meshgrid of filter of required size

D0=20;
D=sqrt((U-size(U,1)/2).^2+(V-size(V,2)/2).^2);
H=exp(-(D.^2)./(2*(D0^2)));
%gaussian lpf with D0=20

S=log(1+abs((H)));
figure,imshow(S,[])
title('Magnitude Spectrum of Filter')
%displaying magnitude spectrum of filter

figure,mesh(abs((H)))
title('Mesh Plot of Magnitude Spectrum of Filter')
colormap('default')
%displaying mesh plot for magnitude spectrum of filter

IMG=fftshift(fft2(double(image),size(H,1),size(H,2)));
imshow(mat2gray(log(1+abs(IMG))))
title('Magnitude Spectrum of Original Image');
%performing fft on image and displaying magnitude spectrum

%performing hadamard product in fourier domain
fourier_filtered_image=(H).*IMG;
output_image=real(ifft2(ifftshift((fourier_filtered_image))));
output_image=output_image(1:size(image,1),1:size(image,2));
%performing inverse fft on fourier domain product and resizing to bring it
%to same size as input image

Fcf=(fourier_filtered_image);
S=log(1+abs(Fcf));
figure,imshow(S,[])
title('Magnitude Spectrum of Filtered Image')
%displaying magnitude spectrum of filtered image

figure,imshow(output_image,[])
title('Filtered Image via Fourier Domain Filtering Process')
%displaying filtered image via Fourier Domain Filtering Process
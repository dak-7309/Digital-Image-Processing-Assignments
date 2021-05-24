%% Question 1b Solution
disp("Solution to 1b goes here")

% reading input image
image=imread("DakshThapar_2018137_Barbara.bmp");
figure,imshow(image,[]);
title('Original Image')

% degradation filter h
h=[1.6,2.9,0;1.3,1,0;0,0,0];
h=h./sum(sum(h));

% convolution of h and input image
g=conv2(h,image,"full");

%AWGN component added to g
g=uint8(g)+uint8(30*randn(514));
figure,imshow(g,[]);
title('Noisy Image')

%fft to convert into fourier domain
G=fft2(g);

% 0 padding h filter for fft2, and then fft2 to fourier domain
padded_h=padarray(h,[511,511],0,'post');
H=fft2(padded_h);

% laplacian filter
laplacian=[-1,-1,-1; -1,8,-1; -1, -1 ,-1 ];

% 0 padding h filter for fft2, and then fft2 to fourier domain
padded_laplacian=padarray(h,[511,511],0,'post');
L=fft2(padded_laplacian);

%values of gamma to be trialled out 
gamma_array=0.003:0.003:2;

%|H|^2
mod_H_2=abs(H.*H);
%|L|^2
mod_L_2=abs(L.*L);
%H*
conj_H=conj(H);
% array of 1s
ones_2=ones(514,514);

%keeping some large initial values for comparison in loop
j=99999999999999999;
IND=-1;

% grid search: iterating over gamma values to find which value yields least mse 
for i=1:length(gamma_array)
    % the term W is found as is proven in 1a
    denom=mod_H_2+ones(514,514) + (gamma_array(i)*(mod_L_2.*(mod_H_2+ones_2)));
    W=conj_H./denom;
    
    % hadamard product in fourier domain
    output=W.*G;
    % real part of ifft to bring it to spatial domain, and slicing to crop out image of original size 
    output_ifft=real(ifft2(output));
    output_ifft=output_ifft(1:512,1:512);
    % mse value
    mse_val=mean(mean(double(output_ifft)-double(image)).^2);
    
    % if mse value is lesser than the previous one, update the mse and keep
    % track of value of gamma
    if (mse_val<j)
        j=mse_val;
        IND=i;
    end
    
end

best_gamma_value=gamma_array(IND)

% the final expression is independent of K, so for all values, the same minimum 
% mse exists

% best gamma, finding and displaying mse
denom1=mod_H_2+ones(514,514) + (gamma_array(IND)*(mod_L_2.*(mod_H_2+ones_2)));
W1=conj_H./denom1;
output1=W1.*G;
output_ifft1=real(ifft2(output1));
output_ifft1=output_ifft1(1:512,1:512);
mse_bestgamma=mean(mean(double(output_ifft1)-double(image)).^2)

% display denoised image
figure,imshow(output_ifft,[])
title('Denoised Image with best gamma')


% gamma =0, finding and displaying mse
denom2=mod_H_2+ones(514,514);
W2=conj_H./denom2;
output2=W2.*G;
output_ifft2=real(ifft2(output2));
output_ifft2=output_ifft2(1:512,1:512);
mse_0gamma=mean(mean(double(output_ifft2)-double(image)).^2)

% display denoised image
figure,imshow(output_ifft,[])
title('Denoised Image with 0 gamma')

% PSNR values for both cases mentioned above
PSNR_bestgamma=10*log(255*255/mse_bestgamma)
PSNR_0gamma=10*log(255*255/mse_0gamma)

% here we can observe PSNR_0gamma>PSNR_bestgamma, i.e. wiener filter has
% higher PSNR value by a very small margin and hence works better

%% Question 3 Solution
disp("Solution to 3 goes here")

% given matrix
matrix=[1,1,1,1; 0,10,10,1; 0,2,3,1; 0,5,15,8];

% sobel operator filter matrices
sobel_op_x=[-1,-2,-1;0,0,0;1,2,1];
sobel_op_y=[-1,0,1;-2,0,2;-1,0,1];

% flipping matrices by 180 degrees
sobel_op_x=rot90(sobel_op_x,2);
sobel_op_y=rot90(sobel_op_y,2);

% convolution for given matrix and sobel filter matrices
matrix_x=conv2(matrix,sobel_op_x,"same");
matrix_y=conv2(matrix,sobel_op_y,"same");

% magnitude matrix obtained by summing the absolute value of above matrices
Magnitude_matrix=abs(matrix_x)+abs(matrix_y)

% angle matrix obtained by taking tan inverse of matrix_uy over matrix_x
Angle_matrix=atand(matrix_y./matrix_x)

% non max suppression matrix initialised by matrix of 0s
non_max_suppression=zeros(4,4);
% initialising 2nd and 3rd rows and columns by the magnitude matrix, as
% know that the border elements are 0
non_max_suppression(2:3,2:3)=Magnitude_matrix(2:3,2:3);

% determining direction for every coordinate by comparing angles as below
% for angles in -90 degrees and 90 degrees
string="";
for i=2:3
    for j=2:3
        if Angle_matrix(i,j)>=-22.5 && Angle_matrix(i,j)<=22.5
            string="horizontal";
        elseif (Angle_matrix(i,j)>=22.5 && Angle_matrix(i,j)<=67.5) 
            string="right_diag";
        elseif (Angle_matrix(i,j)<=-22.5 && Angle_matrix(i,j)>=-67.5) 
            string="left_diag";
        elseif (Angle_matrix(i,j)>=67.5 || Angle_matrix(i,j)<=-67.5) 
            string="vertical";
        end
        
        % comparing 2 immediate neighbours in corresponding directions and
        % assigning 0 value of magnitude value is less than any of the neighbours in that direction  
        if string=="horizontal"
            if Magnitude_matrix(i,j)<=Magnitude_matrix(i,j+1) || Magnitude_matrix(i,j)<=Magnitude_matrix(i,j-1)
                non_max_suppression(i,j)=0;
            end
            
        elseif string=="vertical"
            if Magnitude_matrix(i,j)<=Magnitude_matrix(i-1,j) || Magnitude_matrix(i,j)<=Magnitude_matrix(i+1,j)
                non_max_suppression(i,j)=0;
            end
        
        elseif string=="right_diag"
            if Magnitude_matrix(i,j)<=Magnitude_matrix(i-1,j+1) || Magnitude_matrix(i,j)<=Magnitude_matrix(i+1,j-1)
                non_max_suppression(i,j)=0;
            end
        
        elseif string=="left_diag"
            if Magnitude_matrix(i,j)<=Magnitude_matrix(i-1,j-1) || Magnitude_matrix(i,j)<=Magnitude_matrix(i+1,j+1)
                non_max_suppression(i,j)=0;
            end 
        end
    end
end

%displayng non max suppression
non_max_suppression


% hysterisis thresholding for 90% and 30% of maximum values of Magnitude
% matrices and comparing values at every index with thresholds and
% assigning 0 if value < threshold
thresh_h=max(max(Magnitude_matrix))*0.9;
thresh_l=max(max(Magnitude_matrix))*0.3;

hysterisis_thresh_h=non_max_suppression;
hysterisis_thresh_l=non_max_suppression;

for i=1:4
    for j=1:4
        if hysterisis_thresh_h(i,j)<thresh_h
            hysterisis_thresh_h(i,j)=0;
        end
        if hysterisis_thresh_l(i,j)<thresh_l
            hysterisis_thresh_l(i,j)=0;
        end
    end
end

hysterisis_thresh_h
hysterisis_thresh_l






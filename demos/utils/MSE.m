
% load image1
data1 = load('demos/images/img1_down.mat');
rgbImage1 = data1.img1_down;

% load image2
data2 = load('demos/images/img2_down.mat');
rgbImage2 = data2.img2_down;


image1Size = size(rgbImage1);  % size of image1
nrow1 = image1Size(1);         % number of rows of image1
ncol1 = image1Size(2);         % number of columns of image1
img1 = rgbImage1(1:nrow1-mod(nrow1,8),1:ncol1-mod(ncol1,8),:);  % resize image1 so that its new dimensions are multiples of 8


image2Size = size(rgbImage2);  % size of image2
nrow2 = image2Size(1);         % number of rows of image2
ncol2 = image2Size(2);         % number of columns of image2
img2 = rgbImage2(1:nrow2-mod(nrow2,8),1:ncol2-mod(ncol2,8),:);  % resize image2 so that its new dimensions are multiples of 8


qScale = [0.1; 0.3; 0.6; 1; 2; 5; 10];  % vector of the quantization scaling factors

n = length(qScale);  % number of factors

MSE1 = zeros(n,1);   % initialize MSE vector of image1
MSE2 = zeros(n,1);   % initialize MSE vector of image2

encBitSize1 = zeros(n,1);   % initialize a vector that will contain the total bit size of image1
encBitSize2 = zeros(n,1);   % initialize a vector that will contain the total bit size of image2
% the bit size of each image is measured by the sum length of the image's huffman encoded streams

for i=1:n   % for each scaling factor

    % encode image1 
    JPEGenc1 = JPEGencode(img1,[4 2 2],qScale(i));

    nblocks1 = size(JPEGenc1,1);   % number of blocks in image1
    
    % calculate the sum huffman stream length of all image1 blocks
    for k = 2:nblocks1
        encBitSize1(i) = encBitSize1(i) + strlength(JPEGenc1{k}.huffStream);
    end 

    % reconstruct image1 by decoding
    imgRec1 = JPEGdecode(JPEGenc1);
    
    s1=size(img1);   % size of image1
    M1=s1(1);        % number of image1 rows
    N1=s1(2);        % number of image1 columns 
    mse = sum(sum((img1-imgRec1).^2))/(M1*N1);  % compute MSE of each image component
    MSE1(i) = mse(1) + mse(2) + mse(3);         % the total MSE is the sum of the above 3 MSEs
    
    % print result in command window
    x1 = sprintf('MSE of image1, qScale=%g : %f',qScale(i),MSE1(i));
    disp(x1)
        
end

% plot the MSE of image1 according to the qScale factor
figure()
h1 = plot(qScale,MSE1,'-o','MarkerSize',3);
grid on
ax=gca; 
ax.GridLineStyle = ':';
ax.GridAlpha = 0.3;
set(h1, 'markerfacecolor', get(h1, 'color'));
xlabel('qScale')
ylabel('MSE')
title('Image 1')

% plot the bit size of image1 according to its MSE
figure()
h1 = plot(MSE1,encBitSize1,'-o','MarkerSize',3);
grid on
ax=gca; 
ax.GridLineStyle = ':';
ax.GridAlpha = 0.3;
set(h1, 'markerfacecolor', get(h1, 'color'));
xlabel('MSE')
ylabel('bit size')
title('Image 1')


for i=1:n   % for each scaling factor
    
    % encode image2
    JPEGenc2 = JPEGencode(img2,[4 4 4],qScale(i));
    
    nblocks2 = size(JPEGenc2,1);   % number of blocks in image2

    % calculate the sum huffman stream length of all image2 blocks
    for k = 2:nblocks2
        encBitSize2(i) = encBitSize2(i) + strlength(JPEGenc2{k}.huffStream);
    end
    
    % reconstruct image2 by decoding
    imgRec2 = JPEGdecode(JPEGenc2);
    
    s2=size(img2);   % size of image2
    M2=s2(1);        % number of image2 rows
    N2=s2(2);        % number of image2 columns 
    mse = sum(sum((img2-imgRec2).^2))/(M2*N2);  % compute MSE of each image component
    MSE2(i) = mse(1) + mse(2) + mse(3);         % the total MSE is the sum of the above 3 MSEs

    % print result in command window
    x2 = sprintf('MSE of image2, qScale=%g : %f',qScale(i),MSE2(i));
    disp(x2)
    
end


% plot the MSE of image2 according to the qScale factor
figure()
h2 = plot(qScale,MSE2,'-o','MarkerSize',3);
grid on
ax=gca; 
ax.GridLineStyle = ':';
ax.GridAlpha = 0.3;
set(h2, 'markerfacecolor', get(h2, 'color'));
xlabel('qScale')
ylabel('MSE')
title('Image 2')

% plot the bit size of image2 according to its MSE
figure()
h2 = plot(MSE2,encBitSize2,'-o','MarkerSize',3);
grid on
ax=gca; 
ax.GridLineStyle = ':';
ax.GridAlpha = 0.3;
set(h2, 'markerfacecolor', get(h2, 'color'));
xlabel('MSE')
ylabel('bit size')
title('Image 2')



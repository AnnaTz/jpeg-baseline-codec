clear all

%load image 1
data1 = load('demos/images/img1_down.mat');
rgbImage1 = data1.img1_down;

%load image 2
data2 = load('demos/images/img2_down.mat');
rgbImage2 = data2.img2_down;


%%%%%%%%% a %%%%%%%%%

%image1

%show image1
figure()
imshow(rgbImage1)
title('a: Original Image1')

%convert image1 to YCbCr format
[imageY1,imageCb1,imageCr1] = convert2ycbcr(rgbImage1,[4 2 2]);

%reconvert the YCbCr format of image1 to RGB format
recRGBimage1 = convert2rgb(imageY1,imageCr1,imageCb1,[4 2 2]);

%show the reconstructed image1
figure()
imshow(recRGBimage1)
title('a: Reconstructed Image1')


%image2

%show image2
figure()
imshow(rgbImage2)
title('a: Original Image2')

%convert image2 to YCbCr format
[imageY2,imageCb2,imageCr2] = convert2ycbcr(rgbImage2,[4 4 4]);

%reconvert the YCbCr format of image2 to RGB format
recRGBimage2 = convert2rgb(imageY2,imageCr2,imageCb2,[4 4 4]);

%show the reconstructed image2
figure()
imshow(recRGBimage2)
title('a: Reconstructed Image2')


%%%%%%%% b %%%%%%%%%

%save the luma quantization table
qTableL = [16 11 10 16 24 40 51 61; ...
           12 12 14 19 26 58 60 55; ...
           14 13 16 24 40 57 69 56; ...
           14 17 22 29 51 87 80 62; ...
           18 22 37 56 68 109 103 77; ...
           24 35 55 64 81 104 113 92; ...
           49 64 78 87 103 121 120 101; ...
           72 92 95 98 112 100 103 99];
     
%save the chroma quantization table
qTableC = [17 18 24 47 99 99 99 99; ...
           18 21 26 66 99 99 99 99; ...
           24 26 56 99 99 99 99 99; ...
           47 66 99 99 99 99 99 99; ...
           99 99 99 99 99 99 99 99; ...
           99 99 99 99 99 99 99 99; ...
           99 99 99 99 99 99 99 99; ...
           99 99 99 99 99 99 99 99];

%image1

% show image1
figure()
imshow(rgbImage1)
title('b: Original Image1')

% convert image1 from RGB to YCbCr format
[imageY1,imageCb1,imageCr1] = convert2ycbcr(rgbImage1,[4 2 2]);

% resize the Y, Cb and Cr components
[imageY1,imageCb1,imageCr1] = resize(imageY1,imageCb1,imageCr1);


qBlocksY = zeros(size(imageY1));    % matrix that will contain the quantized version of the Y component
imageSizeY1 = size(imageY1);    % size of the Y component of image1
nrowY1 = imageSizeY1(1);        % number of Y rows
ncolY1 = imageSizeY1(2);        % number of Y columns
for i=1:8:nrowY1-7                                      % iterate rows and columns of the Y component by 8
    for j=1:8:ncolY1-7                                  % and each time
        block = imageY1(i:i+7,j:j+7,:);                 % extract the next 8x8 block 
        dctBlock = blockDCT(block);                     % apply DCT to the block
        qBlock = quantizeJPEG(dctBlock,qTableL,0.6);    % quantize the dc transformed block (using the luma quantization table)
        qBlocksY(i:i+7,j:j+7,:) = qBlock;               % save the quantized block 
    end
end

% exactly the same procedure for the Cb component
% the only difference is that the chroma quantization table is used 
qBlocksCb = zeros(size(imageCb1));
imageSizeCb1 = size(imageCb1);
nrowCb1 = imageSizeCb1(1);
ncolCb1 = imageSizeCb1(2);
for i=1:8:nrowCb1-7
    for j=1:8:ncolCb1-7
        block = imageCb1(i:i+7,j:j+7,:);
        dctBlock = blockDCT(block);
        qBlock = quantizeJPEG(dctBlock,qTableC,0.6);
        qBlocksCb(i:i+7,j:j+7,:) = qBlock;
    end
end

% exactly the same procedure for the Cr component
% using the chroma quantization table like above
qBlocksCr = zeros(size(imageCr1));
imageSizeCr1 = size(imageCr1);
nrowCr1 = imageSizeCr1(1);
ncolCr1 = imageSizeCr1(2);
for i=1:8:nrowCr1-7
    for j=1:8:ncolCr1-7
        block = imageCr1(i:i+7,j:j+7,:);
        dctBlock = blockDCT(block);
        qBlock = quantizeJPEG(dctBlock,qTableC,0.6);
        qBlocksCr(i:i+7,j:j+7,:) = qBlock;
    end
end

% inverse the above procedure for the Y component  
% the reconstructed blocks will be saved directly to the Y image component
for i=1:8:nrowY1-7                                          % iterate rows and columns of the quantized Y component by 8
    for j=1:8:ncolY1-7                                      % and each time
        qBlock = qBlocksY(i:i+7,j:j+7,:);                   % extract the next 8x8 quantized block 
        dctBlock = dequantizeJPEG(qBlock,qTableL,0.6);      % dequantize the quantized block (using the luma quantization table)
        block = iBlockDCT(dctBlock);                        % apply inverse DCT to the dequantized block
        imageY1(i:i+7,j:j+7,:) = block;                     % save the reconstructed block 
    end
end

% exactly the same inverting procedure for the Cb component
% the only difference is that the chroma quantization table is used 
for i=1:8:nrowCb1-7
    for j=1:8:ncolCb1-7
        qBlock = qBlocksCb(i:i+7,j:j+7,:);
        dctBlock = dequantizeJPEG(qBlock,qTableC,0.6);
        block = iBlockDCT(dctBlock);
        imageCb1(i:i+7,j:j+7,:) = block;
    end
end

% exactly the same inverting procedure for the Cr component
% using the chroma quantization table like above
for i=1:8:nrowCr1-7
    for j=1:8:ncolCr1-7
        qBlock = qBlocksCr(i:i+7,j:j+7,:);
        dctBlock = dequantizeJPEG(qBlock,qTableC,0.6);
        block = iBlockDCT(dctBlock);
        imageCr1(i:i+7,j:j+7,:) = block;
    end
end


% convert the reconstructed image1 from YCbCr to RGB format
recRGBimage1 = convert2rgb(imageY1,imageCr1,imageCb1,[4 2 2]);

% show the reconstructed image1
figure()
imshow(recRGBimage1)
title('b: Reconstructed Image1')


%image2

% show image2
figure()
imshow(rgbImage2)
title('b: Original Image2')


% convert image2 from RGB to YCbCr format
[imageY2,imageCb2,imageCr2] = convert2ycbcr(rgbImage2,[4 4 4]);

% resize the Y, Cb and Cr components
[imageY2,imageCb2,imageCr2] = resize(imageY2,imageCb2,imageCr2);


qBlocksY = zeros(size(imageY2));    % matrix that will contain the quantized version of the Y component
imageSizeY2 = size(imageY2);    % size of the Y component of image2
nrowY2 = imageSizeY2(1);        % number of Y rows
ncolY2 = imageSizeY2(2);        % number of Y columns
for i=1:8:nrowY2-7                                    % iterate rows and columns of the Y component by 8
    for j=1:8:ncolY2-7                                % and each time
        block = imageY2(i:i+7,j:j+7,:);               % extract the next 8x8 block 
        dctBlock = blockDCT(block);                   % apply DCT to the block
        qBlock = quantizeJPEG(dctBlock,qTableL,5);    % quantize the dc transformed block (using the luma quantization table)
        qBlocksY(i:i+7,j:j+7,:) = qBlock;             % save the quantized block 
    end
end

% exactly the same procedure for the Cb component
% the only difference is that the chroma quantization table is used 
qBlocksCb = zeros(size(imageCb2));
imageSizeCb2 = size(imageCb2);
nrowCb2 = imageSizeCb2(1);
ncolCb2 = imageSizeCb2(2);
for i=1:8:nrowCb2-7
    for j=1:8:ncolCb2-7
        block = imageCb2(i:i+7,j:j+7,:);
        dctBlock = blockDCT(block);
        qBlock = quantizeJPEG(dctBlock,qTableC,5);
        qBlocksCb(i:i+7,j:j+7,:) = qBlock;
    end
end

% exactly the same procedure for the Cr component
% using the chroma quantization table like above
qBlocksCr = zeros(size(imageCr2));
imageSizeCr2 = size(imageCr2);
nrowCr2 = imageSizeCr2(1);
ncolCr2 = imageSizeCr2(2);
for i=1:8:nrowCr2-7
    for j=1:8:ncolCr2-7
        block = imageCr2(i:i+7,j:j+7,:);
        dctBlock = blockDCT(block);
        qBlock = quantizeJPEG(dctBlock,qTableC,5);
        qBlocksCr(i:i+7,j:j+7,:) = qBlock;
    end
end

% inverse the above procedure for the Y component  
% the reconstructed blocks will be saved directly to the Y image component
for i=1:8:nrowY2-7                                         % iterate rows and columns of the quantized Y component by 8
    for j=1:8:ncolY2-7                                     % and each time 
        qBlock = qBlocksY(i:i+7,j:j+7,:);                  % extract the next 8x8 quantized block 
        dctBlock = dequantizeJPEG(qBlock,qTableL,5);       % dequantize the quantized block (using the luma quantization table)
        block = iBlockDCT(dctBlock);                       % apply inverse DCT to the dequantized block
        imageY2(i:i+7,j:j+7,:) = block;                    % save the reconstructed block 
    end
end

% exactly the same inverting procedure for the Cb component
% the only difference is that the chroma quantization table is used 
for i=1:8:nrowCb2-7
    for j=1:8:ncolCb2-7
        qBlock = qBlocksCb(i:i+7,j:j+7,:);
        dctBlock = dequantizeJPEG(qBlock,qTableC,5);
        block = iBlockDCT(dctBlock);
        imageCb2(i:i+7,j:j+7,:) = block;
    end
end

% exactly the same inverting procedure for the Cr component
% using the chroma quantization table like above
for i=1:8:nrowCr2-7       
    for j=1:8:ncolCr2-7
        qBlock = qBlocksCr(i:i+7,j:j+7,:);
        dctBlock = dequantizeJPEG(qBlock,qTableC,5);
        block = iBlockDCT(dctBlock);
        imageCr2(i:i+7,j:j+7,:) = block;
    end
end


% convert the reconstructed image2 from YCbCr to RGB format
recRGBimage2 = convert2rgb(imageY2,imageCr2,imageCb2,[4 4 4]);

% show the reconstructed image2
figure()
imshow(recRGBimage2)
title('b: Reconstructed Image2')



function [imageY,imageCb,imageCr] = resize(imgY,imgCb,imgCr)

imageSizeY = size(imgY);    %size of the Y component of the image
nrowY = imageSizeY(1);      %number of Y rows
ncolY = imageSizeY(2);      %number of Y columns
imageY = imgY(1:nrowY-mod(nrowY,8),1:ncolY-mod(ncolY,8),:);  %resize Y so that its new dimensions are multiples of 8

imageSizeCb = size(imgCb);      %size of the Cb component of the image
nrowCb = imageSizeCb(1);        %number of Cb rows
ncolCb = imageSizeCb(2);        %number of Cb columns
imageCb = imgCb(1:nrowCb-mod(nrowCb,8),1:ncolCb-mod(ncolCb,8),:);   %resize Cb so that its new dimensions are multiples of 8

imageSizeCr = size(imgCr);      %size of the Cr component of the image
nrowCr = imageSizeCr(1);        %number of Cr rows
ncolCr = imageSizeCr(2);        %number of Cr columns
imageCr = imgCr(1:nrowCr-mod(nrowCr,8),1:ncolCr-mod(ncolCr,8),:);   %resize Cr so that its new dimensions are multiples of 8

end

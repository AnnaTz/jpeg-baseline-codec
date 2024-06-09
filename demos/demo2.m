
% load image1
data1 = load('demos/images/img1_down.mat');
rgbImage1 = data1.img1_down;

size1 = size(rgbImage1);    % size of image1
M1 = size1(1);              % number of image1 rows
N1 = size1(2);              % number of image1 columns

% load image2
data2 = load('demos/images/img2_down.mat');
rgbImage2 = data2.img2_down;

size2 = size(rgbImage2);    % size of image2
M2 = size2(1);              % number of image2 rows
N2 = size2(2);              % number of image2 columns


rgbEntropy1 = 0;   % initialize RGB entropy of image1
for i=1:3                                                          % for each component of the RGB image1 
    rgbEntropy1 = rgbEntropy1 + matrixEntropy(rgbImage1(:,:,i));   % compute its entropy and add it to the total entropy 
end

% compute the total entropy of the RGB image
% multiply the mean entropy by the number of symbols it contains
rgbEntropy1 = rgbEntropy1 * M1 * N1

% convert image1 from RGB to YCbCr format
[imgY1,imgCb1,imgCr1] = convert2ycbcr(rgbImage1,[4 2 2]);

% compute entropy for the quantized dct coefficients and for the run-lengths
[dctEntropy1,rleEntropy1] = encEntropy(imgY1,imgCb1,imgCr1,0.6)


rgbEntropy2 = 0;   % initialize RGB entropy of image2
for i=1:3                                                          % for each component of the RGB image2 
    rgbEntropy2 = rgbEntropy2 + matrixEntropy(rgbImage1(:,:,i));   % compute its entropy and add it to the total entropy 
end

% compute the total entropy of the RGB image
% multiply the mean entropy by the number of symbols it contains
rgbEntropy2 = rgbEntropy2 * M2 * N2

% convert image2 from RGB to YCbCr format
[imgY2,imgCb2,imgCr2] = convert2ycbcr(rgbImage2,[4 4 4]);

% compute entropy for the quantized dct coefficients and for the run-lengths
[dctEntropy2,rleEntropy2] = encEntropy(imgY2,imgCb2,imgCr2,5)


function [dctEntropy,rleEntropy] = encEntropy(imgY,imgCb,imgCr,qScale)

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

  
% the dcts and the run lengths of all the Y blocks are saved in two vectors
% so that their total frequences can be computed afterwards
% likewise for all Cb blocks and for all Cr blocks

% initialize vectors that will contain the dct coefficients of Y, Cb, Cr components
dctY = [];
dctCb = [];
dctCr = [];

% initialize vectors that will contain the run length symbols of Y, Cb, Cr components
rlY = [];
rlCb = [];
rlCr = [];


imageSizeY = size(imgY);   % size of the Y component 
nrowY = imageSizeY(1);     % number of Y rows
ncolY = imageSizeY(2);     % number of Y columns
    
for i=1:8:nrowY-7      % iterate rows and columns of the Y component by 8
    for j=1:8:ncolY-7   
        block = imgY(i:i+7,j:j+7,:);   % extract the next 8x8 block 
          
        dctBlock = blockDCT(block);    % apply DCT to the block
        
        dctY = [dctY; uint8(dctBlock(:))];    % save the dct coefficients of the block in the total vector
        
        % quantize the dc transformed block 
        % (using the luma quantization table and the appropriate qScale factor according to the image)
        qBlock = quantizeJPEG(dctBlock,qTableL,qScale);
        
        % compute the run lengths of the quantized dc transformed block
        if i==1 && j==1
            runSymbols = runLength(qBlock,0);        % for the first block use 0 as the predictor for the DC differential encoding
            DCpred = qBlock(1,1);                    % save the block's DC element as the predictor for the next block
        else
            runSymbols = runLength(qBlock,DCpred);   % for the rest of the blocks use the saved predictor
            DCpred = qBlock(1,1);                    % save the block's DC element as the predictor for the next block
        end  
        
        rlY = [rlY; runSymbols];    % save the run lengths of the block in the total vector
        
    end
end

% exactly the same procedure for the Cb component
% the only difference is that the chroma quantization table is used 
imageSizeCb = size(imgCb);
nrowCb = imageSizeCb(1);
ncolCb = imageSizeCb(2);
for i=1:8:nrowCb-7
    for j=1:8:ncolCb-7
        block = imgCb(i:i+7,j:j+7,:);
        
        dctBlock = blockDCT(block);
        
        dctCb = [dctCb; uint8(dctBlock(:))];
        
        qBlock = quantizeJPEG(dctBlock,qTableC,qScale);
        
        if i==1 && j==1
            runSymbols = runLength(qBlock,0);
            DCpred = qBlock(1,1);
        else
            runSymbols = runLength(qBlock,DCpred);
            DCpred = qBlock(1,1);
        end      
        
        rlCb = [rlCb; runSymbols];
                
    end
end

% exactly the same procedure for the Cr component
% using the chroma quantization table like above
imageSizeCr = size(imgCr);
nrowCr = imageSizeCr(1);
ncolCr = imageSizeCr(2);
for i=1:8:nrowCr-7
    for j=1:8:ncolCr-7
        block = imgCr(i:i+7,j:j+7,:);
        
        dctBlock = blockDCT(block);
        
        dctCr = [dctCr; uint8(dctBlock(:))];
        
        qBlock = quantizeJPEG(dctBlock,qTableC,qScale);

        if i==1 && j==1
            runSymbols = runLength(qBlock,0);
            DCpred = qBlock(1,1);
        else
            runSymbols = runLength(qBlock,DCpred);
            DCpred = qBlock(1,1);
        end   
        
        rlCr = [rlCr; runSymbols];
                
    end
end

% compute the total entropy of the dct coefficients 
% multiply each mean entropy by the number of symbols it contains and then sum the 3 entropies
dctEntropy = matrixEntropy(dctY) * size(dctY,1) + matrixEntropy(dctCb) * size(dctCb,1) + matrixEntropy(dctCr) * size(dctCr,1);

% compute the total entropy of the run lengths 
% multiply each mean entropy by the number of symbols it contains and then sum the 3 entropies
rleEntropy = runLengthEntropy(rlY) * size(rlY,1) + runLengthEntropy(rlCb) * size(rlCb,1) + runLengthEntropy(rlCr) * size(rlCr,1);

end


function H = matrixEntropy(X)

% assemble observed alphabet
Alphabet = unique(X);

% initialize vector of frequencies
Frequency = zeros(size(Alphabet));

% calculate sample frequencies
for symbol = 1:length(Alphabet)
    Frequency(symbol) = sum(sum(X == Alphabet(symbol)),2);
end

% calculate sample class probabilities
P = Frequency / sum(Frequency);

% calculate Shannon entropy
H = -sum(P .* log2(P));

end


function H = runLengthEntropy(X)

% assemble observed alphabet
Alphabet = unique(X,'rows');

sizeA = size(Alphabet);  % size of the alphabet

% initialize vector of frequencies
Frequency = zeros(sizeA(1));

% calculate sample frequencies
for symbol = 1:sizeA(1)
    Frequency(symbol) = sum(sum(ismember(X,Alphabet(symbol,:)),2) == 2);
end

% calculate sample class probabilities
P = Frequency / sum(Frequency);

% calculate Shannon entropy
H = -sum(P .* log2(P));

end





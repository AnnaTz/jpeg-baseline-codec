% function that plots the entropy of the run length symbols according to qScale

% load image1
data1 = load('demos/images/img1_down.mat');
rgbImage1 = data1.img1_down;


% load image2
data2 = load('demos/images/img2_down.mat');
rgbImage2 = data2.img2_down;


qScale = [0.1; 0.3; 0.6; 1; 2; 5; 10];  % vector of the quantization scaling factors

n = length(qScale);  % number of factors

rlEntropy1 = zeros(n,1);     % initialize vector that will contain the entropy of the run length symbols of image1
rlEntropy2 = zeros(n,1);     % initialize vector that will contain the entropy of the run length symbols of image2

for i=1:n   % for each scaling factor


    % convert image1 from RGB to YCbCr format
    [imgY1,imgCb1,imgCr1] = convert2ycbcr(rgbImage1,[4 4 4]);

    % compute entropy for the run-lengths
    rleEnt1 = encEntropy(imgY1,imgCb1,imgCr1,qScale(i));

    % save the computed entropy
    rlEntropy1(i) = rleEnt1;
    

    % convert image2 from RGB to YCbCr format
    [imgY2,imgCb2,imgCr2] = convert2ycbcr(rgbImage2,[4 4 4]);

    % compute entropy for the run-lengths
    rleEnt2 = encEntropy(imgY2,imgCb2,imgCr2,qScale(i));
  
    % save the computed entropy
    rlEntropy2(i) = rleEnt2;
    
end

% plot the entropy of the run length symbols of image1 according to qScale
figure()
h = plot(qScale,rlEntropy1,'-o','MarkerSize',3);
grid on
ax=gca; 
ax.GridLineStyle = ':';
ax.GridAlpha = 0.3;
set(h, 'markerfacecolor', get(h, 'color'));
xlabel('qScale')
ylabel('Entropy of run-length symbols')
title('Image 1')

% plot the entropy of the run length symbols of image2 according to qScale
figure()
h = plot(qScale,rlEntropy2,'-o','MarkerSize',3);
grid on
ax=gca; 
ax.GridLineStyle = ':';
ax.GridAlpha = 0.3;
set(h, 'markerfacecolor', get(h, 'color'));
xlabel('qScale')
ylabel('Entropy of run-length symbols')
title('Image 2')


function rleEntropy = encEntropy(imgY,imgCb,imgCr,qScale)

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

  
% the run lengths of all the Y blocks are saved in two vectors
% so that their total frequences can be computed afterwards
% likewise for all Cb blocks and for all Cr blocks


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

% compute the total entropy of the run lengths 
% multiply each mean entropy by the number of symbols it contains and then sum the 3 entropies
rleEntropy = runLengthEntropy(rlY) * size(rlY,1) + runLengthEntropy(rlCb) * size(rlCb,1) + runLengthEntropy(rlCr) * size(rlCr,1);

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





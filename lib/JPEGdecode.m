function imgRec = JPEGdecode(JPEGenc)

% DC and AC tables will be global variables and will contain either luma or chroma values as needed
global DCtable ACtable

% extract luma and chroma quantization tables from the cell's first structure 
qTableL = JPEGenc{1,1}.qTableL;
qTableC = JPEGenc{1,1}.qTableC;

% extract luma and chroma DC and AC tables from the cell's first structure 
DCL = {JPEGenc{1,1}.DCL};
DCC = {JPEGenc{1,1}.DCC};
ACL = {JPEGenc{1,1}.ACL};
ACC = {JPEGenc{1,1}.ACC};


sizeJPEGenc = size(JPEGenc);   % number of structures that the cell contains

% initialize the maximum horizontal and vertical indexes for Y, Cb and Cr 
maxiHorY = 0;
maxiVerY = 0;
maxiHorCb = 0;
maxiVerCb = 0;
maxiHorCr = 0;
maxiVerCr = 0;

for i=2:sizeJPEGenc(1)              % for every block struct in the cell
    blkType = JPEGenc{i,1}.blkType; % extract the block's type
    indHor = JPEGenc{i,1}.indHor;   % extract the block's horizontal index 
    indVer = JPEGenc{i,1}.indVer;   % extract the block's vertical index
    if strcmp(blkType,'Y')          % if it's a Y block
        if indHor > maxiHorY        % update the maximum Y horizontal index
            maxiHorY = indHor;
        end
        if indVer > maxiVerY        % update the maximum Y vertical index
            maxiVerY = indVer;
        end
    elseif strcmp(blkType,'Cb')     % if it's a Cb block
        if indHor > maxiHorCb       % update the maximum Cb horizontal index
            maxiHorCb = indHor;
        end
        if indVer > maxiVerCb       % update the maximum Cb vertical index
            maxiVerCb = indVer;
        end
    elseif strcmp(blkType,'Cr')     % if it's a Cr block
        if indHor > maxiHorCr       % update the maximum Cr horizontal index
            maxiHorCr = indHor;
        end
        if indVer > maxiVerCr       % update the maximum Cr vertical index
            maxiVerCr = indVer;
        end
    else
        error('Invalid block type')
    end   
end


% initialize the Y, Cb, Cr image components, of type uint8 
% and size based on the maximum indexes that were calculated above
imgY = zeros(maxiHorY*8,maxiVerY*8,'uint8');
imgCb = zeros(maxiHorCb*8,maxiVerCb*8,'uint8');
imgCr = zeros(maxiHorCr*8,maxiVerCr*8,'uint8');

% initialize flags that will indicate if it's the first Y,Cb or Cr block being decoded
fY = 0;
fCb = 0;
fCr = 0;

for i=2:sizeJPEGenc(1)                     % for every block struct in the cell
    blkType = JPEGenc{i,1}.blkType;        % extract the block's type
    indHor = JPEGenc{i,1}.indHor;          % extract the block's horizontal index
    indVer = JPEGenc{i,1}.indVer;          % extract the block's vertical index
    huffStream = JPEGenc{i,1}.huffStream;  % extract the block's huffman encoded stream
    
    if strcmp(blkType,'Y')      % if it's a Y block,
        DCtable = DCL;          % update the DC and AC tables with luma values
        ACtable = ACL;   
    else                        % else,
        DCtable = DCC;          % update the DC and AC tables with chroma values
        ACtable = ACC;
    end  
    
    runSymbols = huffDec(huffStream);    % decode the huffman stream to run length sumbols
    
    % compute the quantized block from the run length symbols
    if strcmp(blkType,'Y')                               % if it's a Y block,
        if fY == 0                                       % check if it's the first one
            qBlock = irunLength(runSymbols, 0);          % for the first block use 0 as the predictor for the DC differential encoding
            DCpredY = qBlock(1,1);                       % save the quantized block's DC element as the predictor for the next block
        else
            qBlock = irunLength(runSymbols, DCpredY);    % for the rest of the blocks use the saved predictor
            DCpredY = qBlock(1,1);                       % save the quantized block's DC element as the predictor for the next block
        end
    elseif strcmp(blkType,'Cb')                          % likewise if it's a Cb block
        if fCb == 0
            qBlock = irunLength(runSymbols, 0);
            DCpredCb = qBlock(1,1);
        else    
            qBlock = irunLength(runSymbols, DCpredCb);
            DCpredCb = qBlock(1,1);
        end
    elseif strcmp(blkType,'Cr')                          % likewise if it's a Cr block
        if fCr == 0
            qBlock = irunLength(runSymbols, 0);
            DCpredCr = qBlock(1,1);
        else
            qBlock = irunLength(runSymbols, DCpredCr);
            DCpredCr = qBlock(1,1);
        end
    end

    % dequantize the quantized block 
    if strcmp(blkType,'Y')                              % if it's a Y block, 
        dctBlock = dequantizeJPEG(qBlock, qTableL, 1);  % use the luma quantization table 
    else                                                % else,
        dctBlock = dequantizeJPEG(qBlock, qTableC, 1);  % use the chroma quantization table
    end
    
    block = iBlockDCT(dctBlock);   % apply inverse DCT to the dequantized block
    
    if strcmp(blkType,'Y')                                              % if it's a Y block, 
        imgY((indHor-1)*8+1:indHor*8,(indVer-1)*8+1:indVer*8) = block;  % save the reconstructed block in the Y image component 
                                                                        % in the position corresponding to the  blocks' indexes
        fY = 1;                                                         % and update the Y flag
    elseif strcmp(blkType,'Cb')                                         % likewise if it's a Cb block
        imgCb((indHor-1)*8+1:indHor*8,(indVer-1)*8+1:indVer*8) = block;
        fCb = 1;
    elseif strcmp(blkType,'Cr')                                         % likewise if it's a Cr block
        imgCr((indHor-1)*8+1:indHor*8,(indVer-1)*8+1:indVer*8) = block;
        fCr = 1;
    end        
end


if maxiVerY == maxiVerCb        % if the Y and Cb have the same number of rows,
    subimg = [4 4 4];           % [4 4 4] subsampling has been used
elseif maxiHorY == maxiHorCb    % else if the Y and Cb have the same number of columns,
    subimg = [4 2 2];           % [4 2 2] subsampling has been used
else                            % else,
    subimg = [4 2 0];           % [4 2 0] subsampling has been used
end


% convert the Y,Cb,Cr components to the reconstructed RGB image
imgRec = convert2rgb(imgY,imgCr,imgCb,subimg);


end


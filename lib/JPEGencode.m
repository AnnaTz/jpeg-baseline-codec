function JPEGenc = JPEGencode(img,subimg,qScale)


imageSize = size(img);   % size of the image
nrow = imageSize(1);     % number of image rows
ncol = imageSize(2);     % number of image columns
img = img(1:nrow-mod(nrow,8),1:ncol-mod(ncol,8),:);  % resize the image so that its new dimensions are multiples of 8

% convert the image from RGB to YCbCr format
[imgY,imgCb,imgCr] = convert2ycbcr(img,subimg);

imageSizeY = size(imgY);   % size of the Y image component
nrowY = imageSizeY(1);     % number of Y rows
ncolY = imageSizeY(2);     % number of Y columns

% same for the Cb component
imageSizeCb = size(imgCb); 
nrowCb = imageSizeCb(1);
ncolCb = imageSizeCb(2);

% and same for the Cr component
imageSizeCr = size(imgCr);
nrowCr = imageSizeCr(1);
ncolCr = imageSizeCr(2);

% the number of encoded structs is equal to the total number of 8x8 blocks of the Y, Cb and Cr components
N = nrowY/8 * ncolY/8 + nrowCb/8 * ncolCb/8 + nrowCr/8 * ncolCr/8;

JPEGenc = cell(N+1,1);  % initialize the cell that will contain the encoded structures

% DC and AC tables will be global variables and will contain either luma or chroma values as needed
global DCtable ACtable

% luma quantization table 
qTableL = [16 11 10 16 24 40 51 61; ...
           12 12 14 19 26 58 60 55; ...
           14 13 16 24 40 57 69 56; ...
           14 17 22 29 51 87 80 62; ...
           18 22 37 56 68 109 103 77; ...
           24 35 55 64 81 104 113 92; ...
           49 64 78 87 103 121 120 101; ...
           72 92 95 98 112 100 103 99];
       
% chroma quantization table
qTableC = [17 18 24 47 99 99 99 99; ...
           18 21 26 66 99 99 99 99; ...
           24 26 56 99 99 99 99 99; ...
           47 66 99 99 99 99 99 99; ...
           99 99 99 99 99 99 99 99; ...
           99 99 99 99 99 99 99 99; ...
           99 99 99 99 99 99 99 99; ...
           99 99 99 99 99 99 99 99];
 
% specification of typical tables for luma DC difference coding
K3 = [00 01 05 01 01 01 01 01 01 00 00 00 00 00 00 00 ...
      00 01 02 03 04 05 06 07 08 09 0x0A 0x0B];

% specification of typical tables for chroma DC difference coding
K4 = [00 03 01 01 01 01 01 01 01 01 01 00 00 00 00 00 ...
      00 01 02 03 04 05 06 07 08 09 0x0A 0x0B];

% specification of typical tables for luma AC difference coding
K5 = [0x00 0x02 0x01 0x03 0x03 0x02 0x04 0x03 0x05 0x05 0x04 0x04 0x00 0x00 0x01 0x7D ...
      0x01 0x02 0x03 0x00 0x04 0x11 0x05 0x12 0x21 0x31 0x41 0x06 0x13 0x51 0x61 0x07 ...
      0x22 0x71 0x14 0x32 0x81 0x91 0xA1 0x08 0x23 0x42 0xB1 0xC1 0x15 0x52 0xD1 0xF0 ...
      0x24 0x33 0x62 0x72 0x82 0x09 0x0A 0x16 0x17 0x18 0x19 0x1A 0x25 0x26 0x27 0x28 ...
      0x29 0x2A 0x34 0x35 0x36 0x37 0x38 0x39 0x3A 0x43 0x44 0x45 0x46 0x47 0x48 0x49 ...
      0x4A 0x53 0x54 0x55 0x56 0x57 0x58 0x59 0x5A 0x63 0x64 0x65 0x66 0x67 0x68 0x69 ...
      0x6A 0x73 0x74 0x75 0x76 0x77 0x78 0x79 0x7A 0x83 0x84 0x85 0x86 0x87 0x88 0x89 ...
      0x8A 0x92 0x93 0x94 0x95 0x96 0x97 0x98 0x99 0x9A 0xA2 0xA3 0xA4 0xA5 0xA6 0xA7 ...
      0xA8 0xA9 0xAA 0xB2 0xB3 0xB4 0xB5 0xB6 0xB7 0xB8 0xB9 0xBA 0xC2 0xC3 0xC4 0xC5 ...
      0xC6 0xC7 0xC8 0xC9 0xCA 0xD2 0xD3 0xD4 0xD5 0xD6 0xD7 0xD8 0xD9 0xDA 0xE1 0xE2 ...
      0xE3 0xE4 0xE5 0xE6 0xE7 0xE8 0xE9 0xEA 0xF1 0xF2 0xF3 0xF4 0xF5 0xF6 0xF7 0xF8 ...
      0xF9 0xFA];
  
% specification of typical tables for chroma AC difference coding
K6 = [0x00 0x02 0x01 0x02 0x04 0x04 0x03 0x04 0x07 0x05 0x04 0x04 0x00 0x01 0x02 0x77 ...
      0x00 0x01 0x02 0x03 0x11 0x04 0x05 0x21 0x31 0x06 0x12 0x41 0x51 0x07 0x61 0x71 ...
      0x13 0x22 0x32 0x81 0x08 0x14 0x42 0x91 0xA1 0xB1 0xC1 0x09 0x23 0x33 0x52 0xF0 ...
      0x15 0x62 0x72 0xD1 0x0A 0x16 0x24 0x34 0xE1 0x25 0xF1 0x17 0x18 0x19 0x1A 0x26 ...
      0x27 0x28 0x29 0x2A 0x35 0x36 0x37 0x38 0x39 0x3A 0x43 0x44 0x45 0x46 0x47 0x48 ...
      0x49 0x4A 0x53 0x54 0x55 0x56 0x57 0x58 0x59 0x5A 0x63 0x64 0x65 0x66 0x67 0x68 ...
      0x69 0x6A 0x73 0x74 0x75 0x76 0x77 0x78 0x79 0x7A 0x82 0x83 0x84 0x85 0x86 0x87 ...
      0x88 0x89 0x8A 0x92 0x93 0x94 0x95 0x96 0x97 0x98 0x99 0x9A 0xA2 0xA3 0xA4 0xA5 ...
      0xA6 0xA7 0xA8 0xA9 0xAA 0xB2 0xB3 0xB4 0xB5 0xB6 0xB7 0xB8 0xB9 0xBA 0xC2 0xC3 ...
      0xC4 0xC5 0xC6 0xC7 0xC8 0xC9 0xCA 0xD2 0xD3 0xD4 0xD5 0xD6 0xD7 0xD8 0xD9 0xDA ...
      0xE2 0xE3 0xE4 0xE5 0xE6 0xE7 0xE8 0xE9 0xEA 0xF2 0xF3 0xF4 0xF5 0xF6 0xF7 0xF8 ...
      0xF9 0xFA];
  
  
% the AC and DC tables are constructed based on the specifications above
% more details are explained inside the local functions GenCodeTab, findHsize and ConstrHcode
  
% construct the luma DC table
Bits = K3(1:16);
HuffVal = K3(17:end);
[Huffsize,HuffVal,HuffCode] = ConstrHcode(Bits, HuffVal);
DCL = {dec2hex(HuffVal,2) Huffsize HuffCode'};

% construct the chroma DC table 
Bits = K4(1:16);
HuffVal = K4(17:end);
[Huffsize,HuffVal,HuffCode] = ConstrHcode(Bits, HuffVal);
DCC = {dec2hex(HuffVal,2) Huffsize HuffCode'};

% construct the luma AC table 
Bits = K5(1:16);
HuffVal = K5(17:end);
[Huffsize,HuffVal,HuffCode] = ConstrHcode(Bits, HuffVal);
ACL = {dec2hex(HuffVal,2) Huffsize HuffCode'};

% construct the chroma AC table 
Bits = K6(1:16);
HuffVal = K6(17:end);
[Huffsize,HuffVal,HuffCode] = ConstrHcode(Bits, HuffVal);
ACC = {dec2hex(HuffVal,2) Huffsize HuffCode'};


% the quantization tables are multiplied by the qScale factor before being encoded, 
% since the qScale will not be contained in the encoded cell
qTableL = qTableL.*qScale;
qTableC = qTableC.*qScale;

% the first struct contains the quantization tables and the DC and AC tables, for luma and chroma
JPEGenc{1,1} = struct('qTableL',qTableL,'qTableC',qTableC,'DCL',DCL,'DCC',DCC,'ACL',ACL,'ACC',ACC);


% update global tables to luma DC and AC values
DCtable = DCL;
ACtable = ACL;

h = 1;      % initialize the horizontal index of the encoded blocks
v = 1;      % initialize the vertical index of the encoded blocks
k = 2;      % initialize the index of the structs inside the cell
for i=1:8:nrowY-7        % iterate rows and columns of the Y component by 8
    for j=1:8:ncolY-7
        block = imgY(i:i+7,j:j+7,:);   % extract the next 8x8 block 
        
        dctBlock = blockDCT(block);    % apply DCT to the block
        
        qBlock = quantizeJPEG(dctBlock,qTableL,1);    % quantize the dc transformed block (using the luma quantization table)
        
        % compute the run lengths of the quantized dc transformed block
        if i==1 && j==1
            runSymbols = runLength(qBlock,0);        % for the first block use 0 as the predictor for the DC differential encoding
            DCpred = qBlock(1,1);                    % save the block's DC element as the predictor for the next block
        else
            runSymbols = runLength(qBlock,DCpred);   % for the rest of the blocks use the saved predictor
            DCpred = qBlock(1,1);                    % save the block's DC element as the predictor for the next block
        end
        
        huffStream = huffEnc(runSymbols);  % apply huffman encoding to the run length symbols
        
        % save the block type (which is Y), the horizontal and vertical
        % indexes of the block and its huffman encoded stream 
        JPEGenc{k,1} = struct('blkType','Y','indHor',h,'indVer',v,'huffStream',huffStream);
        
        k = k + 1;   % increment the structs' index 
        v = v + 1;   % increment the vertical blocks' index
    end
    h = h + 1;      % increment the horizontal blocks' index
    v = 1;          % re-initialize the vertical blocks' index
end

% update global tables to chroma DC and AC values
DCtable = DCC;
ACtable = ACC;

% exactly the same procedure for the Cb component
% the only difference is that the chroma quantization table is used 
% and that the saved block type is Cb
h = 1;
v = 1;
for i=1:8:nrowCb-7
    for j=1:8:ncolCb-7
        block = imgCb(i:i+7,j:j+7,:);
        
        dctBlock = blockDCT(block);
        qBlock = quantizeJPEG(dctBlock,qTableC,1);
        if i==1 && j==1
            runSymbols = runLength(qBlock,0);
            DCpred = qBlock(1,1);
        else
            runSymbols = runLength(qBlock,DCpred);
            DCpred = qBlock(1,1);
        end
        huffStream = huffEnc(runSymbols);
        
        JPEGenc{k,1} = struct('blkType','Cb','indHor',h,'indVer',v,'huffStream',huffStream);
        
        k = k + 1;       
        v = v + 1;
    end
    h = h + 1;
    v = 1;
end

% exactly the same procedure for the Cr component
% using the chroma quantization table like above
% and with the saved block type being Cr
h = 1;
v = 1;
for i=1:8:nrowCr-7
    for j=1:8:ncolCr-7
        block = imgCr(i:i+7,j:j+7,:);
        
        dctBlock = blockDCT(block);
        qBlock = quantizeJPEG(dctBlock,qTableC,1);
        if i==1 && j==1
            runSymbols = runLength(qBlock,0);
            DCpred = qBlock(1,1);
        else
            runSymbols = runLength(qBlock,DCpred);
            DCpred = qBlock(1,1);
        end
        huffStream = huffEnc(runSymbols);
        
        JPEGenc{k,1} = struct('blkType','Cr','indHor',h,'indVer',v,'huffStream',huffStream);
        
        k = k + 1;       
        v = v + 1;
    end
    h = h + 1;
    v = 1;
end

end


% HuffCode: The Huffman code in decimal format 
% HuffVal: The array of symbol values 
% Bits: The array containing the number of codes of each size (1-16)

function HuffCode = GenCodeTab(HuffSize)
HuffSize = [HuffSize; 0];                       % add a zero at the end of the HuffSize array
k = 1;                                          % initialize the index of HuffSize
code = 0;                                       % the first code is always 0
si = HuffSize(1);                               % retrieve the first element from HuffSize
while(1)
    HuffCode(k) = code;                         % assign code into HuffCode
    code = code + 1;                            % increase the code by 1
    k = k + 1;                                  % move index k to the next
    while HuffSize(k) == si                     % if Huffsize(k) is the same as the previous size
        HuffCode(k) = code;                     % assign the new code into HuffCode
        code = code + 1;                        % increase the code by 1
        k = k + 1;                              % move index k to the next
    end
    if HuffSize(k) == 0                         % if k has reached the end of Huffsize 
        return;                                 % done, so return
    end                                         % otherwise, 
    code = bitand(bitshift(code,1),2^16-1);     % shift code to the left
    si = si + 1;                                % increase the current size si by 1
    while(HuffSize(k)~=si)                      % until si reaches HuffSize(k)
        code = bitand(bitshift(code,1),2^16-1); % shift code to the left
        si = si + 1;                            % increase the current size si by 1
    end
end
end

function Huffsize = findHsize(Bits)
VCsize = find(Bits>0);                                  % find where Bits is not zero
VCcnt = Bits(VCsize);                                   % keep only those values
L = length(VCcnt);                                      % number of non zero Bits' values
Huffsize = [];                                          % initialize empty Huffsize, which will be the list of sizes
for k = 1:L                                             % for each non zero value of Bits
    Huffsize = [Huffsize ones(1,VCcnt(k))*VCsize(k)];   % save the size in the list 
end
Huffsize = Huffsize(:);                                 % turn Huffsize into column array
end

function [Huffsize,HuffVal,HuffCode] = ConstrHcode(Bits, HuffVal)
HuffVal = HuffVal(:);             % turn HuffVal into column array
Bits = Bits(:);                   % turn Bits into column array
Huffsize = findHsize(Bits);       % generate Huffsize from Bits
                                  % Huffsize is the list of the sizes corresponding to the codes in HuffVal array
[Huffsize, I] = sort(Huffsize);   % sort Huffsize
HuffVal = HuffVal(I);             % reorganize HuffVal to match the sorted Huffsize
HuffCode = GenCodeTab(Huffsize);  % Generate Huffman Code table based on the list Huffsize
end


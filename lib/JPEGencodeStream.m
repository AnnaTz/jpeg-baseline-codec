function JPEGencStream = JPEGencodeStream(img,subimg,qScale)

imageSize = size(img);      % size of image
nrow = imageSize(1);        % number of rows
ncol = imageSize(2);        % number of columns
img = img(1:nrow-mod(nrow,8),1:ncol-mod(ncol,8),:);   % resize the image so that its new dimensions are multiples of 8
imageSize = size(img);      % new size of image
nrow = imageSize(1);        % new number of rows
ncol = imageSize(2);        % new number of columns


% all markers and encoded tables will firstly be saved in binary format in a cell

% Start of Image (SOI) marker: FFD8=255,216
SOI = [255 216];

% append SOI to the binary cell
bin = dec2bin(SOI,8);
totbin{1} = bin;

% JFIF marker: FFE0=255,224
%        Length=000,016
%        Identifier: 4A46494600=074,070,073,070,000
%        Version=001,002
%        Units=000
%        Xdensity=000,001      
%        Ydensity=000,001      
%        Xthumbnail=000
%        Ythumbnail=000

APP0 = [255	224	000	016	074	070	073	070	000	001	002	000	000	001	000	001	000	000];

% append APP0 to the binary cell
bin = dec2bin(APP0,8);
totbin{end+1} = bin;


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

       
% the quantization table values are multiplied with a scaling factor before being encoded
qTableL = floor(qTableL * qScale);
qTableC = floor(qTableC * qScale);
% dec2bin will round down any float values, so the quantization tables must be floored before being encoded


maxQtL = max(qTableL(:));   % maximum value of luma quantization table
maxQtC = max(qTableC(:));   % maximum value of chroma quantization table


% According to the JPEG standard, quantization table elements should be 8-bit values for 8-bit sample precision (baseline)
% However for experimentation purposes this codec will support both 8-bit and 16-bit quantization table precision

       
% Define Quantization table marker (luma): FFDB=255,219
%         Length:two bytes that indicate the number of bytes that this header contains=000,067, if the quantization table contains 8-bit values
%                                                                                      000,131, if the quantization table contains 16-bit values
%         Precision=000, for 8-bit precision
%                   001, for 16-bit precision
%         Destination identifier=000
%         Quantization values=016,011,012,014,012,010,016,...,101,103,099 (zigzag)


% save luma DQT marker,length,precision and identifier 
if maxQtL > 255
    DQTluma = [255	219	000 131 001 000];
else
    DQTluma = [255	219	000 067 000 000];
end

% append them to the binary cell, using the right number of bits for every part
bin = dec2bin(DQTluma(1:4),8);
totbin{end+1} = bin; 
bin = dec2bin(DQTluma(5:6),4);
totbin{end+1} = bin; 


% construct the vector of indexes in zig zag order
ind = reshape(1:numel(qTableL), size(qTableL));  % indices of elements
ind = fliplr( spdiags( fliplr(ind) ) );          % get the anti-diagonals
ind(:,1:2:end) = flipud( ind(:,1:2:end) );       % reverse order of odd columns
ind(ind==0) = [];                                % keep non-zero indices
                          
% save the luma quantization table in zig zag order
ZZqTableL = qTableL(ind);


if maxQtL > 255     % if the maximum value is greater than 255, the luma table needs 16-bits/element
    bin = dec2bin(ZZqTableL,16);
else                % else, the table needs 8-bits/element
    bin = dec2bin(ZZqTableL,8);
end

% append the binary quantization values to the cell
totbin{end+1} = bin;


% Define Quantization table marker (chroma): FFDB=255,219
%         Length:two bytes that indicate the number of bytes that this header contains=000,067, if the quantization table contains 8-bit values
%                                                                                      000,131, if the quantization table contains 16-bit values
%         Precision=000, for 8-bit precision
%                   001, for 16-bit precision
%         Destination identifier=001
%         Quantization values=017,018,018,024,021,024,047,...,099,099,099 (zigzag)


% save luma DQT marker,length,precision and identifier 
if maxQtC > 255
    DQTchroma = [255	219	000 131 001 001];
else
    DQTchroma = [255	219	000 067 000 001];
end

% append them to the binary cell, using the right number of bits for every part
bin = dec2bin(DQTchroma(1:4),8);
totbin{end+1} = bin; 
bin = dec2bin(DQTchroma(5:6),4);
totbin{end+1} = bin; 


% construct the vector of indexes in zig zag order
ind = reshape(1:numel(qTableC), size(qTableC));  % indices of elements
ind = fliplr( spdiags( fliplr(ind) ) );          % get the anti-diagonals
ind(:,1:2:end) = flipud( ind(:,1:2:end) );       % reverse order of odd columns
ind(ind==0) = [];                                % keep non-zero indices
                          
% save the chroma quantization table in zig zag order
ZZqTableC = qTableC(ind);

if maxQtC > 255     % if the maximum value is greater than 255, the chroma table needs 16-bits/element
    bin = dec2bin(ZZqTableC,16);
else                % else, the table needs 8-bits/element
    bin = dec2bin(ZZqTableC,8);
end

% append the binary quantization values to the cell
totbin{end+1} = bin;


% Start of frame marker: FFC0=255,192 (Baseline DCT)
%         Length:two bytes that indicate the number of bytes that this header contains=000,017
%         Sample precision=008
%         Y=nrow (variable defined earlier)
%         X=ncol (variable defined earlier)
%         Number of components in the image=003 (color baseline)
%         1st Component parameters:
%            ID=001
%            H and V sampling factors=001,001 (for [4 4 4] subsampling)
%                                     002,001 (for [4 2 2] subsampling)
%                                     002,002 (for [4 2 0] subsampling)
%            Quantization table number=000 (luma)
%         2nd Component parameters:
%            ID=002
%            H and V sampling factors=001,001 (for [4 4 4],[4 2 2] and [4 2 0] subsampling)
%            Quantization table number=001 (chroma)
%         3rd Component parameters:
%            ID=003
%            H and V sampling factors=001,001 (for [4 4 4],[4 2 2] and [4 2 0] subsampling)
%            Quantization table number=001 (chroma)


% save SOF marker, length and sample precision 
SOF = [255 192 000 017 008]; 

% append them to the binary cell, using the right number of bits for every part
bin = dec2bin(SOF,8);
totbin{end+1} = bin; 

% append X and Y to the binary cell, using the right number of bits
bin = dec2bin(nrow,16);
totbin{end+1} = bin; 
bin = dec2bin(ncol,16);
totbin{end+1} = bin; 

% save number of components
SOF = [003];

% save ID, sampling factors and quantization table number for each of the 3 components (as described above)
if(isequal(subimg,[4 4 4]))
    SOF = [SOF 001 001 001 000];
    SOF = [SOF 002 001 001 001];
    SOF = [SOF 003 001 001 001];
elseif(isequal(subimg,[4 2 2]))
    SOF = [SOF 001 002 001 000];
    SOF = [SOF 002 001 001 001];
    SOF = [SOF 003 001 001 001];
elseif(isequal(subimg,[4 2 0]))
    SOF = [SOF 001 002 002 000];
    SOF = [SOF 002 001 001 001];
    SOF = [SOF 003 001 001 001];  
else
    error('Unsupported subsampling factors')
end

% append everything to the binary cell, using the right number of bits for every part

bin = dec2bin(SOF(1:2),8);
totbin{end+1} = bin; 
bin = dec2bin(SOF(3:4),4);
totbin{end+1} = bin; 
bin = dec2bin(SOF(5),8);
totbin{end+1} = bin; 

bin = dec2bin(SOF(6),8);
totbin{end+1} = bin; 
bin = dec2bin(SOF(7:8),4);
totbin{end+1} = bin; 
bin = dec2bin(SOF(9),8);
totbin{end+1} = bin; 

bin = dec2bin(SOF(10),8);
totbin{end+1} = bin; 
bin = dec2bin(SOF(11:12),4);
totbin{end+1} = bin; 
bin = dec2bin(SOF(13),8);
totbin{end+1} = bin; 


% specification of typical tables for luma DC coding
K3 = [00 01 05 01 01 01 01 01 01 00 00 00 00 00 00 00 ...
      00 01 02 03 04 05 06 07 08 09 0x0A 0x0B];

% specification of typical tables for chroma DC coding
K4 = [00 03 01 01 01 01 01 01 01 01 01 00 00 00 00 00 ...
      00 01 02 03 04 05 06 07 08 09 0x0A 0x0B];

% specification of typical tables for luma AC coding
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
  
% specification of typical tables for chroma AC coding
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
  
  
% Define Huffman table marker (DCL):FFC4=255,196
%        Length:two bytes that indicate the number of bytes that this header contains=000,031
%        Class=000
%        Identifier=000
%        Bits=The next 16 bytes from an array of unsigned 1-byte integers whose elements give the number of Huffman codes for each possible code length (1-16).
%        Huffman values=000,001,002,...,010,011
  
% save DHT marker, length, class and identifier
DHTdcl = [255 196 000 031 000 000];  

% append them to the binary cell, using the right number of bits for every part
bin = dec2bin(DHTdcl(1:4),8);
totbin{end+1} = bin;   
bin = dec2bin(DHTdcl(5:6),4);
totbin{end+1} = bin; 

% append bits and huffman values of DCL to the binary cell
bin = dec2bin(K3,8);
totbin{end+1} = bin; 

% Define Huffman table marker (ACL):FFC4=255,196
%        Length:two bytes that indicate the number of bytes that this header contains=000,181
%        Class=001
%        Identifier=000
%        Bits=The next 16 bytes from an array of unsigned 1-byte integers whose elements give the number of Huffman codes for each possible code length (1-16).
%        Huffman values=001,002,003,...,249,250

% save DHT marker, length, class and identifier 
DHTacl = [255 196 000 181 001 000];  

% append them to the binary cell, using the right number of bits for every part
bin = dec2bin(DHTacl(1:4),8);
totbin{end+1} = bin;   
bin = dec2bin(DHTacl(5:6),4);
totbin{end+1} = bin; 

% append bits and huffman values of ACL to the binary cell
bin = dec2bin(K5,8);
totbin{end+1} = bin; 

% Define Huffman table marker (DCC):FFC4=255,196
%        Length:two bytes that indicate the number of bytes that this header contains=000,031
%        Class=000
%        Identifier=001
%        Bits=The next 16 bytes from an array of unsigned 1-byte integers whose elements give the number of Huffman codes for each possible code length (1-16).
%        Huffman values=000,001,002,...,010,011

% save DHT marker, length, class and identifier
DHTdcc = [255 196 000 031 000 001];  

% append them to the binary cell, using the right number of bits for every part
bin = dec2bin(DHTdcc(1:4),8);
totbin{end+1} = bin;   
bin = dec2bin(DHTdcc(5:6),4);
totbin{end+1} = bin; 

% append bits and huffman values of DCC to the binary cell
bin = dec2bin(K4,8);
totbin{end+1} = bin;
  
% Define Huffman table marker (ACC):FFC4=255,196
%        Length:two bytes that indicate the number of bytes that this header contains=000,181
%        Class=001
%        Identifier=001
%        Bits=The next 16 bytes from an array of unsigned 1-byte integers whose elements give the number of Huffman codes for each possible code length (1-16).
%        Huffman values=001,002,003,...,249,250

% save DHT marker, length, class and identifier
DHTacc = [255 196 000 181 001 001];  

% append them to the binary cell, using the right number of bits for every part
bin = dec2bin(DHTacc(1:4),8);
totbin{end+1} = bin;   
bin = dec2bin(DHTacc(5:6),4);
totbin{end+1} = bin; 

% append bits and huffman values of ACC to the binary cell
bin = dec2bin(K6,8);
totbin{end+1} = bin; 


% Start of Scan marker: FFDA=255,218
%        Length:two bytes that indicate the number of bytes that this header contains=000,012
%        Number of components=003
%        1st Component parameters:
%           ID=001
%           DC and AC table numbers=000 and 000 (luma)
%        2nd Component parameters:
%           ID=002
%           DC and AC table numbers=001 and 001 (chroma)
%        3rd Component parameters:
%           ID=003
%           DC and AC table numbers=001 and 001 (chroma)
%        Ss=000 (Baseline DCT)
%        Se=063 (Baseline DCT)
%        Ah and Al=000 (Baseline DCT)

SOS = [255 218 000 012 003 001 000 000 002 001 001 003 001 001 000 063 000 000];

% append SOS to the binary cell, using the right number of bits for every part

bin = dec2bin(SOS(1:6),8);
totbin{end+1} = bin; 
bin = dec2bin(SOS(7:8),4);
totbin{end+1} = bin; 
bin = dec2bin(SOS(9),8);
totbin{end+1} = bin; 
bin = dec2bin(SOS(10:11),4);
totbin{end+1} = bin; 
bin = dec2bin(SOS(12),8);
totbin{end+1} = bin; 
bin = dec2bin(SOS(13:14),4);
totbin{end+1} = bin; 
bin = dec2bin(SOS(15:16),8);
totbin{end+1} = bin; 
bin = dec2bin(SOS(17:18),4);
totbin{end+1} = bin; 

% convert the cell to a vector of 1s and 0s
cbin=char(totbin);
STR=convertCharsToStrings(cbin');
STR = strrep(STR,' ','');
stream=str2double(regexp(num2str(STR),'\d','match'));



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


% convert the image to a binary vector
dataStream = compressData(img,subimg,qTableL,qTableC,DCL,ACL,DCC,ACC);


% append it to the rest of the stream vector
stream=[stream str2double(regexp(num2str(dataStream),'\d','match'))];


% End of Image (EOI) marker: FFD9=255,217
EOI = [255 217];

% append EOI to the binary stream vector
binstr=dec2bin(EOI,8);
STR=convertCharsToStrings(binstr');
stream=[stream str2double(regexp(num2str(STR),'\d','match'))];

% this is the final JPEGencStream
JPEGencStream = stream;


% save the JPEG encoded stream in a jpeg file
fileID = fopen('encImage.jpeg','w');
l = length(stream);
for i=1:8:l-7
    str = stream(i+7:-1:i);
    fwrite(fileID,str,'ubit1');
end
fclose(fileID);


% show the JPEG encoded image
% im = imread('encImage.jpeg','jpeg');
% figure
% imshow(im)

end


function dataStream = compressData(img,subimg,qTableL,qTableC,DCL,ACL,DCC,ACC)

% DC and AC tables will be global variables and will contain either luma or chroma values as needed
global DCtable ACtable

% convert the image from RGB to YCbCr format
[imgY,imgCb,imgCr] = convert2ycbcr(img,subimg);


% typecast to double in order to contain the shifted values
imgY = double(imgY);
imgCb = double(imgCb);
imgCr = double(imgCr);

% apply level shifting by substracting the 128 offset
imgY = imgY - 128;
imgCb = imgCb - 128;
imgCr = imgCr - 128;


imageSizeY = size(imgY);   % size of the Y image component
nrowY = imageSizeY(1);     % number of Y rows
ncolY = imageSizeY(2);     % number of Y columns

% same for the Cb component
imageSizeCb = size(imgCb);
nrowCb = imageSizeCb(1);
ncolCb = imageSizeCb(2);

% same for the Cr component
imageSizeCr = size(imgCr);
nrowCr = imageSizeCr(1);
ncolCr = imageSizeCr(2);


% initialize 3 cells that will contain all Y,Cb,Cr blocks in right order for encoding
Yblocks = cell(nrowY/8*ncolY/8,1);
CBblocks = cell(nrowCb/8*ncolCb/8,1);
CRblocks = cell(nrowCr/8*ncolCr/8,1);

% initialize the 3 indexes for the above cells
indY = 1;
indCb = 1;
indCr = 1;

if(isequal(subimg,[4 2 0]))   % for [4 2 0] subsampling,
    % the data units of the Y component will be ordered starting from the top left corner 
    % zig-zaging across a 4 unit block, then moving the same way across the first row of 4 unit blocks
    % and likewise across the next row of 4 unit blocks etc
    for k = 0:2:nrowY/8-2
        for t = 0:2:ncolY/8-2
            for i = k+1:k+2
                for j = t+1:t+2
                    Yblocks{indY} = imgY((i-1)*8+1:i*8,(j-1)*8+1:j*8);
                    indY = indY + 1;
                end
            end
        end
    end
else                          % for [4 4 4] or [4 2 2] subsampling,
    % the data units of the Y component will be ordered starting from the top left corner,
    % across the first row and then likewise across the second row etc
    for i=1:8:nrowY-7
        for j=1:8:ncolY-7                    
            Yblocks{indY} = imgY(i:i+7,j:j+7,:); 
            indY = indY + 1;
        end
    end
end

% the data units of the Cb component will be ordered starting from the top left corner,
% across the first row and then likewise across the second row etc
for i=1:8:nrowCb-7
    for j=1:8:ncolCb-7   
        CBblocks{indCb} = imgCb(i:i+7,j:j+7,:);
        indCb = indCb + 1;
    end
end

% the data units of the Cr component will be ordered in the same way
for i=1:8:nrowCr-7
    for j=1:8:ncolCr-7 
        CRblocks{indCr} = imgCr(i:i+7,j:j+7,:);
        indCr = indCr + 1;
    end
end


% initialize 3 cells that will contain the ordered Y,Cb,Cr blocks in binary huffman encoded format
y = {};
cb = {};
cr = {};


% update global tables to luma DC and AC values
DCtable = DCL;
ACtable = ACL;

for ind=1:nrowY/8*ncolY/8                       % iterate across the Yblocks cell
    block = Yblocks{ind};                       % extract the next 8x8 block             
    dctBlock = blockDCT(block);                 % apply DCT to the block
    qBlock = quantizeJPEG(dctBlock,qTableL,1);  % quantize the dc transformed block (using the luma quantization table)

    % compute the run lengths of the quantized dc transformed block
    if ind == 1                      
        runSymbols = runLength(qBlock,0);       % for the first block use 0 as the predictor for the DC differential encoding
        DCpred = qBlock(1,1);                   % save the block's DC element as the predictor for the next block
    else
        runSymbols = runLength(qBlock,DCpred);  % for the rest of the blocks use the saved predictor
        DCpred = qBlock(1,1);                   % save the block's DC element as the predictor for the next block
    end        
    huffStream = huffEnc(runSymbols);           % apply huffman encoding to the run length symbols

    y{end+1} = huffStream;                      % save the huffman encoded stream in the y cell 
end


% update global tables to chroma DC and AC values
DCtable = DCC;
ACtable = ACC;

% exactly the same procedure for the Cb component
% the only difference is that the chroma quantization table is used 
for ind=1:nrowCb/8*ncolCb/8                              
        block = CBblocks{ind};       
        dctBlock = blockDCT(block);
        qBlock = quantizeJPEG(dctBlock,qTableC,1);
        
        if ind == 1
            runSymbols = runLength(qBlock,0);
            DCpred = qBlock(1,1);          
        else
            runSymbols = runLength(qBlock,DCpred);
            DCpred = qBlock(1,1);
        end
        huffStream = huffEnc(runSymbols);
        
        cb{end+1} = huffStream;
end

% exactly the same procedure for the Cr component
% using the chroma quantization table like above
for ind=1:nrowCr/8*ncolCr/8                               
        block = CRblocks{ind};
        dctBlock = blockDCT(block);
        qBlock = quantizeJPEG(dctBlock,qTableC,1);
        
        if ind == 1
            runSymbols = runLength(qBlock,0);
            DCpred = qBlock(1,1);
        else
            runSymbols = runLength(qBlock,DCpred);
            DCpred = qBlock(1,1);
        end
        huffStream = huffEnc(runSymbols);
        
        cr{end+1} = huffStream;
end


if(isequal(subimg,[4 4 4]))       % for [4 4 4] subsampling,
    
    % the data units will be interleaved in the sequence: Y - Cb - Cr
    dataStream = [];
    i=1;
    j=1;
    k=1;
    while k<=nrowCr/8*ncolCr/8
        dataStream=[dataStream y{i} cb{j} cr{k}];
        i=i+1;
        j=j+1;
        k=k+1;
    end
    
elseif(isequal(subimg,[4 2 2]))   % for [4 2 2] subsampling,
    
    % the data units will be interleaved in the sequence: YY - Cb - Cr
    dataStream = [];
    i=1;
    j=1;
    k=1;
    while k<=nrowCr/8*ncolCr/8
        dataStream=[dataStream y{i:i+1} cb{j} cr{k}];
        i=i+2;
        j=j+1;
        k=k+1;
    end
    
elseif(isequal(subimg,[4 2 0]))   % for [4 2 0] subsampling,
    
    % the data units will be interleaved in the sequence: YYYY - Cb - Cr
    dataStream = [];
    i=1;
    j=1;
    k=1;
    while k<=nrowCr/8*ncolCr/8
        dataStream=[dataStream y{i:i+3} cb{j} cr{k}];
        i=i+4;
        j=j+1;
        k=k+1;
    end
    
end

% perform byte alignment by padding with 1-bits
l = length(dataStream);
while rem(l,8) ~= 0 
    dataStream = append(dataStream,'1');
    l = length(dataStream);
end

% convert the binary encoded stream to a decimal array
l = length(dataStream);
ar = [];
for i=1:8:l-7
    str = dataStream(i:i+7);
    str = num2str(str);
    ar = [ar bin2dec(str)];
end

k = find(ar == 255);  % find 0xFF occurences
nar = ar(1:k(1));
for i=1:length(k)-1
    nar = [nar 0 ar(k(i)+1:k(i+1))];   % stuff 0x00 after them
end
nar = [nar 0 ar(k(i+1)+1:end)];

% convert the decimal array back to a binary stream
nar=dec2bin(nar);
dataStream=num2str(convertCharsToStrings(nar'));


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


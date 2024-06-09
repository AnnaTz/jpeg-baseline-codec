function imgCmp = JPEGdecodeStream(JPEGencStream)


% convert the input vector of 1s and 0s to a string
stream = num2str(JPEGencStream);
stream =  strrep(stream,' ','');


SOI = [255 216];                    % SOI marker
cbin = dec2bin(SOI,8);              % convert it to binary
str=convertCharsToStrings(cbin');   % convert it to string
i = strfind(stream,str);            % find it in the data string
if isempty(i)
    error('SOI not found')
end
i = i(1);   

DQT = [255	219];                   % DQT marker
cbin = dec2bin(DQT,8);              % convert it to binary
str=convertCharsToStrings(cbin');   % convert it to string
j = strfind(stream(i+16:end),str);  % find it in the data string (after the SOI)
j = j + i + 15;
if isempty(j)
    error('DQT not found')
end
qTable0 = [];   % initialize two empty quantization tables
qTable1 = [];
DQTid = [];     % initialize a vector that will contain the id of the quantization tables

for z = 1:2                                         % for the two quantization tables
    i = j(z);
    i = i + 32;
    dqtprec  = stream(i:i+3);                       % load the quantization table precision
    dqtprec = bin2dec(dqtprec);                     % convert it to decimal
    if dqtprec == 0                                 % translate it to number of bits
        dqtprec = 8;
    elseif dqtprec == 1
        dqtprec = 16;
    else
        error('Invalid DQT precision')
    end
    i = i + 4;
    dqtid = stream(i:i+3);                          % load the quantization table id
    i = i + 4;
    dqtid = bin2dec(dqtid);                         % convert it to decimal
    DQTid = [DQTid dqtid];                          % save it in the id vector
    qTable = stream(i:i+64*dqtprec-1);              % load the quantization table values
    qTable = cellstr(reshape(qTable,dqtprec,[])');  % separate and save them in a binary cell 
    qTable = char(qTable);                          % convert the cell to a char array
    qTable = bin2dec(qTable);                       % convert the array to a decimal vector
    qTable = izigzag(qTable,8,8);                   % create the quantization table by inverting the zig zag order
    if isempty(qTable0)                             % save it in the first empty table
        qTable0 = qTable;
    elseif isempty(qTable1)
        qTable1 = qTable;
    end 
end

i = i + 64*dqtprec;


SOF = [255 192];                    % SOF marker
cbin = dec2bin(SOF,8);              % convert it to binary
str=convertCharsToStrings(cbin');   % convert it to string
j = strfind(stream(i:end),str);     % find it in the data string (after the DQT section)
j = j + i - 1;
if isempty(j)
    error('SOF not found')
end
i = j(1);

i = i + 8*5;
y = stream(i:i+15);     % load the number of rows in the image
nrows = bin2dec(y);     % convert it to decimal
i = i + 16;
x = stream(i:i+15);     % load the number of columns in the image
ncols = bin2dec(x);     % convert it to decimal
i = i + 8*2 - 1;
nf = stream(i:i+8);     % load the number of components
nf = bin2dec(nf);       % convert it to decimal
i = i + 8 + 1;

Cid = [];     % initialize a vector that will contain the components' identifiers
hv = [];      % initialize a vector that will contain the sampling factors
Tq = [];      % initialize a vector that will contain the quantization table selectors

for j=1:nf                  % for every component
    ci = stream(i:i+7);     % load the identifier
    ci = bin2dec(ci);       % convert it to decimal
    Cid = [Cid ci];         % save it in the identifiers' vector
    i = i + 8;

    hi = stream(i:i+3);     % load the horizontal sampling factor
    hi = bin2dec(hi);       % convert it to decimal
    hv = [hv hi];           % save it in the sampling factors' vector
    i = i + 4;

    vi = stream(i:i+3);     % load the vertical sampling factor
    vi = bin2dec(vi);       % convert it to decimal
    hv = [hv vi];           % save it in the sampling factors' vector
    i = i + 4;
  
    tqi = stream(i:i+7);    % load the quantization table selector
    tqi = bin2dec(tqi);     % convert it to decimal
    Tq = [Tq tqi];          % save it in the quantization selectors' vector
    i = i + 8;
end

DHT = [255 196];                     % DHT marker
cbin = dec2bin(DHT,8);               % convert it to binary
str=convertCharsToStrings(cbin');    % convert it to string
j = strfind(stream(i:end),str);      % find it in the data string (after the SOF section)
j = j + i - 1;
if isempty(j)
    error('DCT not found')
end

dhTable1 = {};    % initialize 4 empty huffman entropy coding tables
dhTable2 = {};
dhTable3 = {};
dhTable4 = {};
Tc = [];      % initialize a vector that will contain the tables' class (DC or AC)
Th = [];      % initialize a vector that will contain the tables' identifiers

for z = 1:4                         % for every table
    i = j(z);
    i = i + 2*16;
    tci = stream(i:i+3);            % load the class
    tci = bin2dec(tci);             % convert it to decimal
    Tc = [Tc tci];                  % save it in the classes' vector
    i = i + 4;

    thi  = stream(i:i+3);           % load the identifier
    thi = bin2dec(thi);             % convert it to decimal
    Th = [Th thi];                  % save it in the identifiers' vector
    i = i + 4;
    
    nvals = 0;                      % initialize the total number of huffman codes
    bits = zeros(16,1);             % initialize a vector that will contain the number of codes of each length
    for k=1:16                      % for each possible code length
        ibits = stream(i:i+7);      % load the number of huffman codes with this bit length
        bits(k) = bin2dec(ibits);   % turn it to decimal
        i = i + 8;
        nvals = nvals + bits(k);    % increment the sum number of huffman codes
    end
     
    dhTable = stream(i:i+nvals*8-1);             % load the huffman coding values
    dhTable = cellstr(reshape(dhTable,8,[])');   % separate and save them in a binary cell
    dhTable = char(dhTable);                     % convert the cell to a char array
    dhTable = bin2dec(dhTable);                  % convert the array to a decimal vector

    % construct the huffman coding table using the bits and values
    % more details are explained inside the local functions GenCodeTab, findHsize and ConstrHcode
    Bits = bits;
    HuffVal = dhTable;
    [Huffsize,HuffVal,HuffCode] = ConstrHcode(Bits, HuffVal);
    dhTable = {dec2hex(HuffVal,2) Huffsize HuffCode'};
    
    if isempty(dhTable1)      % save it in the first empty table
        dhTable1 = dhTable;
    elseif isempty(dhTable2)
        dhTable2 = dhTable;
    elseif isempty(dhTable3)
        dhTable3 = dhTable;
    elseif isempty(dhTable4)
        dhTable4 = dhTable;
    end  
end

i = i + nvals*8;


SOS = [255 218];                     % SOS marker
cbin = dec2bin(SOS,8);               % convert it to binary
str=convertCharsToStrings(cbin');    % convert it to string
j = strfind(stream(i:end),str);      % find it in the data string (after the DHT section)
j = j + i - 1;
if isempty(j)
    error('SOS not found')
end
i = j(1);

i = i + 2*16;
ns = stream(i:i+7);    % load the number of image components in the scan
ns = bin2dec(ns);      % convert it to decimal
i = i + 8;

Cs = [];     % initialize a vector that will contain the order of the components in the scan
Td = [];     % initialize a vector that will contain their DC table selectors
Ta = [];     % initialize a vector that will contain their AC table selectors

for j=1:ns
    ci = stream(i:i+7);    % load the identifier of the component that will be next in the scan sequence
    ci = bin2dec(ci);      % convert it to decimal
    Cs = [Cs ci];          % save it in the component ordering table
    i = i + 8;

    tdi = stream(i:i+3);   % load the DC selector of the component
    tdi = bin2dec(tdi);    % convert it to decimal
    Td = [Td tdi];         % save it in the DC selectors' vector
    i = i + 4;

    tai = stream(i:i+3);   % load the AC selector of the component
    tai = bin2dec(tai);    % convert it to decimal
    Ta = [Ta tai];         % save it in the AC selectors' vector
    i = i + 4;
end

if ns ~= 3
    error('This codec only supports 3 components (YCbCr)')
end


oCid = [];
for v=1:ns
    oCid = [oCid find(Cid == Cs(v))];
end

% reorder the quantization table selectors according to the component order in the scan
oTq = Tq(oCid);

% convert the sampling factors' vector to an array containing the horizontal factors in the first column and the vertical in the second
ohv = reshape(hv,2,length(hv)/2);

% reorder the sampling factors according to the component order in the scan
ohv = ohv(:,oCid)';

% translate the sampling factors to [4 4 4], [4 2 2] or [4 2 0] subsampling
if isequal(ohv,[1 1; 1 1; 1 1])
    subimg = [4 4 4];
elseif isequal(ohv,[2 1; 1 1; 1 1])
    subimg = [4 2 2];
elseif isequal(ohv,[2 2; 1 1; 1 1])
    subimg = [4 2 0];
else
    error('Unsupported subsampling factors')
end


% find which of the saved quantization tables is for luminance
z = find(DQTid == oTq(1));
if z == 1
    qTableL = qTable0;
elseif z == 2
    qTableL = qTable1;
end

% find which of the saved quantization tables is for chrominance
z = find(DQTid == oTq(2));
if z ~= find(DQTid == oTq(3))
    error('Both Cb and Cr should be chroma')
end
if z == 1
    qTableC = qTable0;
elseif z == 2
    qTableC = qTable1;
end

% find which of the saved huffman encoding tables is the DCL
z = find(Th == Td(1));
zi = find(Tc(z) == 0);
z = z(zi);
if z == 1
    DCL = dhTable1;
elseif z == 2
    DCL = dhTable2;
elseif z == 3
    DCL = dhTable3;
elseif z == 4
    DCL = dhTable4;
end

% find which of the saved huffman encoding tables is the ACL
z = find(Th == Ta(1));
zi = find(Tc(z) == 1);
z = z(zi);
if z == 1
    ACL = dhTable1;
elseif z == 2
    ACL = dhTable2;
elseif z == 3
    ACL = dhTable3;
elseif z == 4
    ACL = dhTable4;
end

% find which of the saved huffman encoding tables is the DCC
z = find(Th == Td(2));
zi = find(Tc(z) == 0);
z = z(zi);
if z == 1
    DCC = dhTable1;
elseif z == 2
    DCC = dhTable2;
elseif z == 3
    DCC = dhTable3;
elseif z == 4
    DCC = dhTable4;
end

% find which of the saved huffman encoding tables is the ACC
z = find(Th == Ta(2));
zi = find(Tc(z) == 1);
z = z(zi);
if z == 1
    ACC = dhTable1;
elseif z == 2
    ACC = dhTable2;
elseif z == 3
    ACC = dhTable3;
elseif z == 4
    ACC = dhTable4;
end


EOI = [255 217];                     % EOI marker
cbin = dec2bin(EOI,8);               % convert it to binary
str=convertCharsToStrings(cbin');    % convert it to string
j = strfind(stream(i:end),str);      % find it in the data string (after the SOS section)
j = j + i - 1;
if isempty(j)
    error('EOI not found')
end
j = j(end) - 1;     
i = i + 3*8;

% the rest of the binary stream is the image data
streamData = stream(i:j);


% convert the binary encoded stream to a decimal array
l = length(streamData);
ar = [];
for i=1:8:l-7
    str = streamData(i:i+7);
    str = num2str(str);
    ar = [ar bin2dec(str)];
end

k = find(ar == 255);   % find 0xFF occurences
nar = ar(1:k(1));
for i=2:length(k)
    nar = [nar ar(k(i-1)+2:k(i))];   % delete 0x00 after them
end
nar = [nar ar(k(i)+2:end)];

% convert the decimal array back to a binary stream
nar=dec2bin(nar);
streamData=num2str(convertCharsToStrings(nar'));



% decode the compressed image data
imgCmp = decompressData(streamData,subimg,nrows,ncols,qTableL,qTableC,DCL,ACL,DCC,ACC);


end


function img = decompressData(data,subimg,nrows,ncols,qTableL,qTableC,DCL,ACL,DCC,ACC)


ncolsY = ncols;    % the Y component has the same number of columns with the RGB image
nrowsY = nrows;    % and the same number of rows
if(isequal(subimg,[4 4 4]))         % for [4 4 4] subsampling
    ncolsCb = ncols;                % the Cb component has the same number of columns with the RGB image
    nrowsCb = nrows;                % and the same number of rows
    ncolsCr = ncols;                % the Cr component has the same number of columns with the RGB image
    nrowsCr = nrows;                % and the same number of rows
elseif(isequal(subimg,[4 2 2]))     % for [4 2 2] subsampling
    ncolsCb = ncols / 2;            % the Cb component has half the number of columns of the RGB image
    nrowsCb = nrows;                % and the same number of rows
    ncolsCr = ncols / 2;            % the Cr component has half the number of columns with the RGB image
    nrowsCr = nrows;                % and the same number of rows
else                                % else (for [4 2 0] subsampling)
    ncolsCb = ncols / 2;            % the Cb component has half the number of columns of the RGB image
    nrowsCb = nrows / 2;            % and half the number of rows
    ncolsCr = ncols / 2;            % the Cr component has half the number of columns of the RGB image
    nrowsCr = nrows / 2;            % and half the number of rows
end   


% initialize the Y,Cb and Cr image components with size as described above
% and type double in order to contain the level shifted values
imgY = zeros(nrowsY,ncolsY,'double');
imgCb = zeros(nrowsCb,ncolsCb,'double');
imgCr = zeros(nrowsCr,ncolsCr,'double');

% DC and AC tables will be global variables and will contain either luma or chroma values as needed
global DCtable ACtable

% initialize 3 cells that will contain the runSymbols of the Y,Cb and Cr components
Yblocks = {};
CBblocks = {};
CRblocks = {};


if(isequal(subimg,[4 4 4]))         % for [4 4 4] subsampling
    m = 1;                          % 1 Y data unit is being followed by 1 Cb and 1 Cr data unit in the interleaving sequence
elseif(isequal(subimg,[4 2 2]))     % for [4 2 2] subsampling
    m = 2;                          % 2 Y data units are being followed by 1 Cb and 1 Cr data unit in the interleaving sequence
else                                % else (for [4 2 0] subsampling)
    m = 4;                          % 4 Y data units are being followed by 1 Cb and 1 Cr data unit in the interleaving sequence
end


while(1)     % iterating through the data
    
    % update global tables to luma DC and AC values
    DCtable = DCL;
    ACtable = ACL;
    
    for k = 1:m                                 % as many times as defined by the interleaving pattern described above
        [runSymbols,index] = huffDec(data);     % decode the next block's runSymbols
        Yblocks{end+1} = runSymbols;            % save them
        data = data(index:end);                 % keep the rest of the data
    end

    % update global tables to chroma DC and AC values
    DCtable = DCC;
    ACtable = ACC; 
    
    [runSymbols,index] = huffDec(data);     % decode the next block's runSymbols
    CBblocks{end+1} = runSymbols;           % save them
    data = data(index:end);                 % keep the rest of the data

    [runSymbols,index] = huffDec(data);     % decode the next block's runSymbols
    CRblocks{end+1} = runSymbols;           % save them
    if index-1 == strlength(data)           % if the data index has reached the end,
        break                               % stop
    else                                    % else
        data = data(index:end);             % keep the rest of the data
    end
    
    for i = 1:length(data)
        if ~strcmp(data(i),'1')   % check if the remaining data stream consists of 1-bits
            break
        end
    end
    if i==length(data)   % if it only contains 1-bits,
        break            % stop, because those remaining bits are byte alignment padding bits
    end
        
end


if (isequal(subimg,[4 4 4])) || (isequal(subimg,[4 2 2]))     % for [4 4 4] or [4 2 2] subsampling
    k = 1;     % initialize the block index
    for i=1:8:nrowsY-7                       % iterate rows and columns of the Y component by 8
        for j=1:8:ncolsY-7

            runSymbols = Yblocks{k};   % load the runSymbols of the current Y block

            % compute the quantized block from the run length symbols
            if k == 1   
                qBlock = irunLength(runSymbols, 0);           % for the first block use 0 as the predictor for the DC differential encoding
                DCpred = qBlock(1,1);                         % save the quantized block's DC element as the predictor for the next block
            else
                qBlock = irunLength(runSymbols, DCpred);      % for the rest of the blocks use the saved predictor
                DCpred = qBlock(1,1);                         % save the quantized block's DC element as the predictor for the next block
            end
            dctBlock = dequantizeJPEG(qBlock, qTableL, 1);    % dequantize the quantized block, using the luma quantization table
            block = iBlockDCT(dctBlock);                      % apply inverse DCT to the dequantized block

            imgY(i:i+7,j:j+7,:) = block;                      % save the reconstructed block in the Y image component 

            k = k + 1;                                        % increment the block index
        end
    end
else                                      % else (for [4 2 0] subsampling)
    k = 1;     % initialize the block index
    for w = 0:2:nrowsY/8-2
        for t = 0:2:ncolsY/8-2                      % iterate rows and columns of the Y component 
            for i = w+1:w+2                         % in the right order in which the blocks were encoded
                for j = t+1:t+2

                    runSymbols = Yblocks{k};   % load the runSymbols of the current Y block

                    % compute the quantized block from the run length symbols
                    if k == 1   
                        qBlock = irunLength(runSymbols, 0);           % for the first block use 0 as the predictor for the DC differential encoding
                        DCpred = qBlock(1,1);                         % save the quantized block's DC element as the predictor for the next block
                    else
                        qBlock = irunLength(runSymbols, DCpred);      % for the rest of the blocks use the saved predictor
                        DCpred = qBlock(1,1);                         % save the quantized block's DC element as the predictor for the next block
                    end
                    dctBlock = dequantizeJPEG(qBlock, qTableL, 1);    % dequantize the quantized block, using the luma quantization table
                    block = iBlockDCT(dctBlock);                      % apply inverse DCT to the dequantized block

                    imgY((i-1)*8+1:i*8,(j-1)*8+1:j*8) = block;        % save the reconstructed block in the Y image component 

                    k = k + 1;                                        % increment the block index
                end
            end
        end
    end   
end


% exactly the same procedure for the Cb component, as described above in the [4 4 4] or [4 2 2] subsampling case
% the only difference is that the chroma quantization table is used 
k = 1;
for i=1:8:nrowsCb-7
    for j=1:8:ncolsCb-7
       
        runSymbols = CBblocks{k};
        
        if k == 1
            qBlock = irunLength(runSymbols, 0);
            DCpred = qBlock(1,1);
        else
            qBlock = irunLength(runSymbols, DCpred);
            DCpred = qBlock(1,1);
        end
        dctBlock = dequantizeJPEG(qBlock, qTableC, 1);
        block = iBlockDCT(dctBlock);
        
        imgCb(i:i+7,j:j+7,:) = block;
                
        k = k + 1;       
    end
end


% exactly the same procedure for the Cr component
% using the chroma quantization table like above
k = 1;
for i=1:8:nrowsCr-7
    for j=1:8:ncolsCr-7
       
        runSymbols = CRblocks{k};
        
        if k == 1
            qBlock = irunLength(runSymbols, 0);
            DCpred = qBlock(1,1);
        else
            qBlock = irunLength(runSymbols, DCpred);
            DCpred = qBlock(1,1);
        end
        dctBlock = dequantizeJPEG(qBlock, qTableC, 1);
        block = iBlockDCT(dctBlock);
        
        imgCr(i:i+7,j:j+7,:) = block;
         
        k = k + 1;       
    end
end


% invert level shifting by adding the 128 offset
imgY = imgY + 128;
imgCb = imgCb + 128;
imgCr = imgCr + 128;


% convert the image from YCbCr to RGB format
img = convert2rgb(imgY,imgCr,imgCb,subimg);


end

% this is the same procedure as described in the irunLength function
function output = izigzag(in, vmax, hmax)

h = 1;      % initialize the column index of qBlock
v = 1;      % initialize the row index of qBlock

vmin = 1;   % minimum number of row index
hmin = 1;   % minimum number of column index

output = zeros(vmax, hmax);    % initialize the output matrix

i = 1;    % initialize the zzStream element index

while ((v <= vmax) && (h <= hmax))          % while the row and column indexes remain inside the size of the block
    if (mod(h + v, 2) == 0)                 % going up
        if (v == vmin)                      % if it's the first row
            output(v, h) = in(i);           % save the current zzStream element to the current position of qBlock
            if (h == hmax)                  % if it's the last column
                v = v + 1;                  % increment to next row
            else                            % else
                h = h + 1;                  % increment to next column
            end
            i = i + 1;                      % increment to next zzStream element
        elseif ((h == hmax) && (v < vmax))  % if it's the last column and not the last row yet
            output(v, h) = in(i);           % save the current zzStream element to the current position of qBlock
            v = v + 1;                      % increment to next row
            i = i + 1;                      % increment to next zzStream element
        elseif ((v > vmin) && (h < hmax))   % if it's not the first row neither the last column
            output(v, h) = in(i);           % save the current zzStream element to the current position of qBlock
            v = v - 1;                      % decrement to previous row
            h = h + 1;                      % increment to next column
            i = i + 1;                      % increment to next zzStream element
        end
        
    else                                    % going down
       if ((v == vmax) && (h < hmax))       % if it's the last row and not the last column yet
            output(v, h) = in(i);           % save the current zzStream element to the current position of qBlock
            h = h + 1;                      % increment to next column
            i = i + 1;                      % increment to next zzStream element
        
       elseif (h == hmin)                   % if it's the first column
            output(v, h) = in(i);           % save the current zzStream element to the current position of qBlock
            if (v == vmax)                  % if it's the last row
                h = h + 1;                  % increment to next column
            else                            % else
                v = v + 1;                  % increment to next row
            end
            i = i + 1;                      % increment to next zzStream element
       elseif ((v < vmax) && (h > hmin))    % if it's not the last row yet neither the first column
            output(v, h) = in(i);           % save the current zzStream element to the current position of qBlock
            v = v + 1;                      % increment to next row
            h = h - 1;                      % decrement to previous column
            i = i + 1;                      % increment to next zzStream element
       end
    end
    if ((v == vmax) && (h == hmax))         % if it's the last column and last row
        output(v, h) = in(i);               % save the current zzStream element to the current position of qBlock
        break                               % and stop
    end
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


% this is almost the same procedure as in huffDec.m
% the difference is that, along with the runSymbols matrix, this huffDec returns an index to the start of the next block's stream in huffStream
% this is needed to separately decode each block's binary stream, since all have been appended continuously in the huffStream
function [runSymbols,index] = huffDec(huffStream)


% DC and AC tables are global and currently have luma or chroma values as needed
global DCtable ACtable

index = -1;    % initialize index that will point after this block's stream in huffStream

runSymbols = repmat(-1,64,2);   % initialize runSymbols matrix with -1 values
                                % this is done to save time by avoid dynamic expanding inside loops
                                % unused rows will be erased at the end of the function

q = 1;    % initialize row index for the runSymbols 

streamlen=length(huffStream);    % length of huffStream


symb = DCtable{:,1};    % symbols of the DC table
clen = DCtable{:,2};    % bit lengths of the DC symbols
code = DCtable{:,3};    % codes of the DC symbols (in decimal format)

n=length(symb);     % number of DC symbols
maxlen=max(clen);   % the maximum bit length of DC codes

bincodes = cell(n,1);   % cell that will contain all the binary codes
for m = 1:n
    bincodes{m} = dec2bin(code(m),clen(m));   % save all codes in binary format with the corresponding bit length
end

j=0;    % initialize the length of the part that will be selected from HuffStream
while j<maxlen                                              % gradually expand the selected part of HuffStream 
                                                            % until it reaches the maximum length of possible binary codes
    
    c=huffStream(1:1+j);                                    % get the stream part
    ind=1;                                                  % initialize the index of the table of binary codes
    while (ind<=n && ~isequal(bincodes{ind},c))             % search a binary code equal to the current stream part 
        ind=ind+1;                                          % increment the table's index
    end
    if ind<=n                                               % if equality has been found
        rs = symb(ind,:);                                   % get the run/size symbol of the recognised code
        if strcmp(rs,'00')                                  % if the rs is '00'
            runSymbols(q,:) = [0 0];                        % then just add EOB to the runSymbols (no vli has been encoded)
            q = q + 1;                                      % increment the row index of runSymbols
            cat = 0;                                        % the category is zero
        else                                                % in all other cases
            cat = hex2dec(rs(2));                           % the second hex letter of the rs is the category
            vli = huffStream(2+j:2+j+cat-1);                % get the x next bits of huffstream (where x is the category)
                                                            % this is the variable-length integer
            if cat ~= (floor(log2(abs(bin2dec(vli)))) + 1)  % if the category formula doesn't give the right result,
                                                            % it means that the vli came from a negative number
                stuff = '';                                 % initialize the stuffing string that will make vli 16-bits
                for k = 1:(16-cat)
                    stuff = append(stuff,'1');              % create the binary string of 1s 
                end
                vli = append(stuff,vli);                    % and append it to the start of vli
                
                                                            % if the vli came from a negative number, 
               dc = bin2dec(vli) - 65535;                   % the dc is the signed integer version of vli + 1
            else                                            % else
               dc = bin2dec(vli);                           % the ac is simply the integer version of vli
            end
            runSymbols(q,:) = [0 dc];                       % save the dc as a length in the runSymbols matrix, along with a zero run
            q = q + 1;                                      % increment the row index of runSymbols
        end
        break                                               % stop increasing the size of the selected stream part
    else                                                    % else
        j=j+1;                                              % increment the length of the selected HuffStream part
    end
end


symb = ACtable{:,1};    % symbols of the AC table
clen = ACtable{:,2};    % bit lengths of the AC symbols
code = ACtable{:,3};    % codes of the AC symbols (in decimal format)

n=length(symb);     % number of AC symbols
maxlen=max(clen);   % the maximum bit length of AC codes

bincodes = cell(n,1);   % cell that will contain all the binary codes
for m = 1:n
    bincodes{m} = dec2bin(code(m),clen(m));   % save all codes in binary format with the corresponding bit length
end

counter = 1;    % counter of the encoded coefficients currently included in runSymbols 
                % initialized to 1, since runSymbols already includes the dc coefficient

i = 2+j+cat;    % save the huffStream's index after the dc huffman code
while (i<=streamlen)                                            % gradually shift the selected part of HuffStream
                                                                % until it reaches the end of this block's stream
    j=0;                                                        % initialize the length of the part 
    while j<maxlen                                              % gradually expand the selected part of HuffStream 
                                                                % until it reaches the maximum length of possible binary codes
        c=huffStream(i:i+j);                                    % get the stream part
        ind=1;                                                  % initialize the index of the table of binary codes
        while (ind<=n && ~isequal(bincodes{ind},c))             % search a binary code equal to the current stream part 
            ind=ind+1;                                          % increment the table's index
        end
        if ind<=n                                               % if equality has been found
            rs = symb(ind,:);                                   % get the run/size symbol of the recognised code
            if strcmp(rs,'00')                                  % if the rs is '00'
                runSymbols(q,:) = [0 0];                        % then just add EOB to the runSymbols (no vli has been encoded)
                q = q + 1;                                      % increment the row index of runSymbols
                cat = 0;                                        % the category is zero
                index = i;                                      % save index of where EOB was found
            elseif strcmp(rs,'F0')                              % if the rs is 'F0'
                runSymbols(q,:) = [15 0];                       % then just add ZRL to the runSymbols (no vli has been encoded)
                q = q + 1;                                      % increment the row index of runSymbols
                cat = 0;                                        % the category is zero
            else                                                % in all other cases
                cat = hex2dec(rs(2));                           % the second hex letter of the rs is the category
                vli = huffStream(i+j+1:i+j+1+cat-1);            % get the x next bits of huffstream (where x is the category)
                                                                % this is the variable-length integer
                if cat ~= (floor(log2(abs(bin2dec(vli)))) + 1)  % if the category formula doesn't give the right result,
                                                                % it means that the vli came from a negative number
                    stuff = '';                                 % initialize the stuffing string that will make vli 16-bits
                    for k = 1:(16-cat)
                        stuff = append(stuff,'1');              % create the binary string of 1s 
                    end
                    vli = append(stuff,vli);                    % and append it to the start of vli
                    
                                                                % if the vli came from a negative number, 
                   ac = bin2dec(vli) - 65535;                   % the ac is the signed integer version of vli + 1
                else                                            % else
                   ac = bin2dec(vli);                           % the ac is simply the integer version of vli
                end
                runSymbols(q,:) = [hex2dec(rs(1)) ac];          % save the pair of run-length in the runSymbols matrix
                                                                % the first hex letter of the rs is the zero run
                                                                % and the ac number is the length
                q = q + 1;                                      % increment the row index of runSymbols
            end
            
            counter = counter + runSymbols(q-1,1) + 1;          % increment the counter of the encoded coefficients in runSymbols
            
            break                                               % stop increasing the size of the selected stream part
        else                                                    % else
            j=j+1;                                              % increment the length of the selected HuffStream part
        end
    end
    if (counter == 64) || (index ~= -1)                         % if EOB has been found or runSymbols includes 64 encoded coefficients,
        index = i+j+1+cat;                                      % save where the next block's stream will start in huffStream
        break                                                   % stop
    end
    i=i+j+1+cat;                                                % increment the stream index to the position after the explored part
end

% erase all [-1 -1] rows that have remained in the final runSymbols matrix
runSymbols((runSymbols(:,1) == -1) & (runSymbols(:,2) == -1),:) = [];
runSymbols = reshape(runSymbols, [], 2);

end

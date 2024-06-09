function huffStream = huffEnc(runSymbols)

% DC and AC tables are global and currently have luma or chroma values as needed
global DCtable ACtable

huffStream = '';  % initialize huffStream as an empty string

sizeRunSymbols = size(runSymbols);  % size of the run-length symbols matrix


symbols = DCtable{:,1};   % symbols of the DC table
clens = DCtable{:,2};     % bit lengths of the DC symbols
codes = DCtable{:,3};     % codes of the DC symbols (in decimal format)

dc = runSymbols(1,2);   % load the dc from the first pair of run-length symbols

% calculate the dc's category
if dc == 0 
    cat = 0;  
else
    cat = floor(log2(abs(dc))) + 1;
end
if cat > 11 
     error('Huffman DC category out of bounds')
end


n = length(codes);  % number of DC codes

% find the DC symbol that corresponds to the category
for j=1:n
    if symbols(j,:) == sprintf('%X%X',0,cat)
        break
    end
end

% convert the corresponding code to binary format with the corresponding bit length
vlc = dec2bin(codes(j),clens(j));   % (variable-length code)

% append vlc to the huffStream
huffStream = append(huffStream,vlc);

% if dc is not zero, then we need to append vli (variable-length integer) to the vlc
if dc ~= 0
    if dc < 0                                               % if dc < 0 (and is of category x), then take the x least significant
        vli = sprintf('%d',bitget(int16(dc-1),cat:-1:1));   % bits of the two's complement representation of dc - 1
    else                                                    % else, take the x least significant bits of the
        vli = sprintf('%d',bitget(int16(dc),cat:-1:1));     % binary representation of dc
    end

    % append vli to the huffStream
    huffStream = append(huffStream,vli);
end


symbols = ACtable{:,1};   % symbols of the AC table
clens = ACtable{:,2};     % bit lengths of the AC symbols
codes = ACtable{:,3};     % codes of the AC symbols (in decimal format)

n = length(codes);  % number of AC codes

for i=2:sizeRunSymbols   % for each run-length symbol
    
    ac = runSymbols(i,2);   % load the ac from the current pair of run-length symbols
    
    % calculate the ac's category
    if ac == 0
        cat = 0;
    else
        cat = floor(log2(abs(ac))) + 1;
    end
    if cat > 10
         error('Huffman AC category out of bounds')
    end

    % find the AC symbol that corresponds to the category
    for j=1:n
        if symbols(j,:) == sprintf('%X%X',runSymbols(i,1),cat)
            break
        end
    end

    % convert the corresponding code to binary format with the corresponding bit length
    vlc = dec2bin(codes(j),clens(j));   % (variable-length code)

    % append vlc to the huffStream
    huffStream = append(huffStream,vlc);
    
    % if ac is not zero, then we need to append vli (variable-length integer) to the vlc
    if ac ~= 0 
        if ac < 0                                               % if ac < 0 (and is of category x), then take the x least significant 
            vli = sprintf('%d',bitget(int16(ac-1),cat:-1:1));   % bits of the two's complement representation of ac - 1
        else                                                    % else, take the x least significant bits of the
            vli = sprintf('%d',bitget(int16(ac),cat:-1:1));     % binary representation of ac
        end

        % append vli to the huffStream
        huffStream = append(huffStream,vli);
    end
    
end
end


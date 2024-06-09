function runSymbols = huffDec(huffStream)

% DC and AC tables are global and currently have luma or chroma values as needed
global DCtable ACtable

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


i = 2+j+cat;    % save the huffStream's index after the dc huffman code
while i<=streamlen                                              % gradually shift the selected part of HuffStream
                                                                % until it reaches the end
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
            break                                               % stop increasing the size of the selected stream part
        else                                                    % else
            j=j+1;                                              % increment the length of the selected HuffStream part
        end
    end
    i=i+j+1+cat;                                                % increment the stream index to the position after the explored part
end

% erase all [-1 -1] rows that have remained in the final runSymbols matrix
runSymbols((runSymbols(:,1) == -1) & (runSymbols(:,2) == -1),:) = [];
runSymbols = reshape(runSymbols, [], 2);

end

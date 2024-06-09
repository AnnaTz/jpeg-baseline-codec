function qBlock = irunLength(runSymbols,DCpred)


sizeSymbols = size(runSymbols);     % size of the runSymbols matrix
R = sizeSymbols(1);     % number of rows of runSymbols, which is the number of run-length pairs

j=1;    % initialize the zzStream (zig zag stream) index
for i=1:R                         % for each [run length] pair of runSymbols
    for k=0:runSymbols(i,1)-1
        zzStream(j)=0;            % add a 'run' number of zeros to the zzStream
        j=j+1; 
    end
    zzStream(j)=runSymbols(i,2);  % add the 'length' element to the zzStream
    j=j+1;
end

% additional zeros might be needed after the EOB, so
while length(zzStream) < 64     % if zzStream isn't full yet, 
    zzStream(end+1) = 0;        % add zeros until full
end
                                      % the first element is differentially encoded, so 
zzStream(1) = zzStream(1) + DCpred;   % add the given DC predictor



h = 1;      % initialize the column index of qBlock
v = 1;      % initialize the row index of qBlock

vmin = 1;   % minimum number of row index
hmin = 1;   % minimum number of column index
vmax=8;     % maximum number of row index
hmax=8;     % maximum number of column index

qBlock = zeros(vmax, hmax);    % initialize the 8x8 qBlock

i = 1;    % initialize the zzStream element index

while ((v <= vmax) && (h <= hmax))          % while the row and column indexes remain inside the size of the block
    if (mod(h + v, 2) == 0)                 % going up
        if (v == vmin)                      % if it's the first row
            qBlock(v, h) = zzStream(i);     % save the current zzStream element to the current position of qBlock
            if (h == hmax)                  % if it's the last column
                v = v + 1;                  % increment to next row
            else                            % else
                h = h + 1;                  % increment to next column
            end
            i = i + 1;                      % increment to next zzStream element
        elseif ((h == hmax) && (v < vmax))  % if it's the last column and not the last row yet
            qBlock(v, h) = zzStream(i);     % save the current zzStream element to the current position of qBlock
            v = v + 1;                      % increment to next row
            i = i + 1;                      % increment to next zzStream element
        elseif ((v > vmin) && (h < hmax))   % if it's not the first row neither the last column
            qBlock(v, h) = zzStream(i);     % save the current zzStream element to the current position of qBlock
            v = v - 1;                      % decrement to previous row
            h = h + 1;                      % increment to next column
            i = i + 1;                      % increment to next zzStream element
        end
        
    else                                    % going down
       if ((v == vmax) && (h < hmax))       % if it's the last row and not the last column yet
            qBlock(v, h) = zzStream(i);     % save the current zzStream element to the current position of qBlock
            h = h + 1;                      % increment to next column
            i = i + 1;                      % increment to next zzStream element
        
       elseif (h == hmin)                   % if it's the first column
            qBlock(v, h) = zzStream(i);     % save the current zzStream element to the current position of qBlock
            if (v == vmax)                  % if it's the last row
                h = h + 1;                  % increment to next column
            else                            % else
                v = v + 1;                  % increment to next row
            end
            i = i + 1;                      % increment to next zzStream element
       elseif ((v < vmax) && (h > hmin))    % if it's not the last row yet neither the first column
            qBlock(v, h) = zzStream(i);     % save the current zzStream element to the current position of qBlock
            v = v + 1;                      % increment to next row
            h = h - 1;                      % decrement to previous column
            i = i + 1;                      % increment to next zzStream element
       end
    end
    if ((v == vmax) && (h == hmax))         % if it's the last column and last row
        qBlock(v, h) = zzStream(i);         % save the current zzStream element to the current position of qBlock
        break                               % and stop
    end
end

end



function runSymbols = runLength(qBlock,DCpred)

% construct the vector of indexes in zig zag order
ind = reshape(1:numel(qBlock), size(qBlock));    % indices of elements
ind = fliplr( spdiags( fliplr(ind) ) );          % get the anti-diagonals
ind(:,1:2:end) = flipud( ind(:,1:2:end) );       % reverse order of odd columns
ind(ind==0) = [];                                % keep non-zero indices


zzStream = round(qBlock(ind));   % save the zig zag stream after rounding the qBlock elements

                                                  % the first element must be differentially encoded, so 
zzStream(1) = round(qBlock(1)) - round(DCpred);   % substract the given DC predictor


% the run length symbols are initially saved in a vector

% save the first zig zag element in the second position of the vector, with a zero run before it
c(1)=0;
c(2) = zzStream(1);

if zzStream(2)~=0     % if the first ac zig zag element is not zero,
    c(3)=0;           % save a zero run in the run-length vector before starting the iteration
    j=4;              % and initialize the run-length vector's index 
else                  % else
    j=3;              % just initialize the index accordingly
end
a=length(zzStream);     % length of the zig zag stream
count=0;    % counter of continuous zeros
flag=0;     % initialize the flag that will indicate if there where zeros preceding an element
for n=2:a                                       % for every zig zag element (after the dc)
    b=zzStream(n);                              % load the current zzStream element
    if b==0                                     % if it's zero
        if n==a                                 % if it's the last element
            count=count+1;                      % increase the zeros counter
            c(j)=count;                         % and save its value in the run-length vector
        elseif zzStream(n)==zzStream(n+1)       % else if the next zzStream element is zero as well
            count=count+1;                      % just increase the zeros counter
        else                                    % else  (the next zzStream element is a non zero number)
            count=count+1;                      % increase the zeros counter
            c(j)=count;                         % save its value in the run-length vector
            j=j+1;                              % increment the run-length vector's index
            count=0;                            % re-initialize the zeros counter
        end
        flag=0;                                 % set the flag to preceding zeros 
    else                                        % else if the current zzStream element is not zero
        if flag==1                              % if there were no preceding zeros
            c(j)=0;                             % save a zero run in the run-length vector
            j=j+1;                              % increment the run-length vector's index
        end
        c(j)=b;                                 % save the non zero element in the run-length vector
        j=j+1;                                  % increment the run-length vector's index
        flag=1;                                 % set the flag to no preceding zeros
    end
end
if flag==0                   % if the run-length vector ended in zeros
    c(end) = c(end) - 1;     % substract a zero from the counter
    c(end+1) = 0;            % and save it as a length in the run-length vector
end

% the elements in odd positions of the runSymbols vector will be the first column of the runSymbols matrix
% and those in even positions will be the second column
c1=c(1:2:length(c));
c2=c(2:2:length(c));
runSymbols=[c1' c2'];

if runSymbols(end,2) == 0           % no matter how many zeros there were at the end, 
    runSymbols(end,:) = [0 0];      % there will be only one EOB
end

% adjust the runSymbols matrix so that the maximum size of a run is 15

runSymb=[;];    % initialize a new empty runSymbols matrix
lasti = 0;      % initialize the index that will point to the last found long run
for i = 1:size(runSymbols,1)                          % for every run-length pair
    if runSymbols(i,1) > 15                           % if there is a run longer than 15
        if isempty(runSymb)                           % if it's the first long run found in the matrix
            runSymb(1:i-1,:) = runSymbols(1:i-1,:);   % save all the previous run-lengths to the new matrix
            k = i;                                    % and update the index of the found long run in the new matrix
        else                                          % else
            k = k + i - lasti - 1;                    % just update the index accordingly
        end
        int = fix(runSymbols(i,1)/16);                % how many times is the found long run greater than 16 ([15 0] -> 15 + 1)
        rem = mod(runSymbols(i,1),16);                % how many zeros remain from the division
        for j=0:int-1   
            runSymb(k,:) = [15 0];                    % add the longest possible run in the new matrix as many times as needed
            k = k + 1;                                % increment the new matrix's index
        end
        runSymb(k,:) = [rem runSymbols(i,2)];         % save the pair of the few remaining zeros and the length
        k = k + 1;                                    % increment the new matrix's index
        runSymb(k:k+size(runSymbols,1)-i-1,:) = runSymbols(i+1:size(runSymbols,1),:);   % save all the next run-lengths to the new matrix
        lasti = i;                                    % save the index to where the long run was found in the initial matrix
    end
end

 
if ~isempty(runSymb)        % if changes were made,
    runSymbols = runSymb;   % update the runSymbols matrix
end

end


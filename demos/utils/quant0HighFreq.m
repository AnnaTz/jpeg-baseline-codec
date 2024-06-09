% load image1
data1 = load('demos/images/img1_down.mat');
rgbImage1 = data1.img1_down;

% load image2
data2 = load('demos/images/img2_down.mat');
rgbImage2 = data2.img2_down;


% In order to zero the blocks' 20,40,50,60 or 63 highest frequency elements
% each of those elements will be divided by a large number contained at the same position in the quantization tables.
% To alter the quantization tables without changing the rest of the project's functions
% those large numbers will firstly be saved in qMatrix
% and then passed to JPEGencode as qScale (since qTables are getting multiplied by qScale anyway)
% the true qScale value for the lower frequency elements (those that will not be zeroed) will be 1


% create the qScale that will zero the 20 highest frequency elements
v20 = ones(1,64);
v20(45:64) = 1000000000;
qMatrix(:,:,1) = izigzag(v20, 8, 8);

% create the qScale that will zero the 40 highest frequency elements
v40 = ones(1,64);
v40(25:64) = 1000000000;
qMatrix(:,:,2) = izigzag(v40, 8, 8);

% create the qScale that will zero the 50 highest frequency elements
v50 = ones(1,64);
v50(15:64) = 1000000000;
qMatrix(:,:,3) = izigzag(v50, 8, 8);

% create the qScale that will zero the 60 highest frequency elements
v60 = ones(1,64);
v60(5:64) = 1000000000;
qMatrix(:,:,4) = izigzag(v60, 8, 8);

% create the qScale that will zero the 63 highest frequency elements
v63 = ones(1,64);
v63(2:64) = 1000000000;
qMatrix(:,:,5) = izigzag(v63, 8, 8);


for i = 1:5     % for every case of experimentation

% load qScale
qScale = qMatrix(:,:,i);


% encode image1 
JPEGenc1 = JPEGencode(rgbImage1,[4 4 4],qScale);

% reconstruct image1 by decoding
imgRec1 = JPEGdecode(JPEGenc1);

% show reconstructed image1
figure()
imshow(imgRec1)


% encode image2
JPEGenc2 = JPEGencode(rgbImage2,[4 4 4],qScale);

% reconstruct image2 by decoding
imgRec2 = JPEGdecode(JPEGenc2);

% show reconstructed image2
figure()
imshow(imgRec2)

end


% this is the same procedure as described in the irunLength function
% and also used in JPEGdecodeStream.m
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

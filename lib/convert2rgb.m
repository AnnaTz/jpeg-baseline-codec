function imageRGB = convert2rgb(imageY,imageCr,imageCb,subimg)


image_size = size(imageY);  % size of the Y image component
nrow = image_size(1);       % number of Y rows
ncol = image_size(2);       % number of Y columns


if isequal(subimg,[4 4 4])          % if the subsampling factors are [4 4 4]
    imageYCbCr(:,:,1) = imageY;     % all rows and columns of the Y component are given
    imageYCbCr(:,:,2) = imageCb;    % all rows and columns of the Cb component are given
    imageYCbCr(:,:,3) = imageCr;    % all rows and columns of the Cr component are given
    
elseif isequal(subimg,[4 2 2])                      % if the subsampling factors are [4 2 2]
    imageYCbCr(:,:,1) = imageY;                     % all rows and columns of the Y component are given
    imageYCbCr(:,1:2:ncol-1,2) = imageCb(:,:);      % only the odd columns of the Cb component are given
    imageYCbCr(:,1:2:ncol-1,3) = imageCr(:,:);      % only the odd columns of the Cr component are given
    
    % compute each even column of Cb and Cr as the interpolation of the previous and next column
    imageYCbCr(:,2:2:ncol-2,2:3) = (double(imageYCbCr(:,1:2:ncol-3,2:3))+ double(imageYCbCr(:,3:2:ncol-1,2:3)))/2;
    % the last even column doesn't have a next one so it's just equal to its previous column
    imageYCbCr(:,ncol,2:3) = imageYCbCr(:,ncol-1,2:3);
      
elseif isequal(subimg,[4 2 0])                              % if the subsampling factors are [4 2 0]
    imageYCbCr(:,:,1) = imageY;                             % all rows and columns of the Y component are given
    imageYCbCr(1:2:nrow-1,1:2:ncol-1,2) = imageCb(:,:);     % only the odd rows and columns of the Cb component are given
    imageYCbCr(1:2:nrow-1,1:2:ncol-1,3) = imageCr(:,:);     % only the odd rows and columns of the Cr component are given
    
    % for the odd rows of Cb and Cr
    % compute each even column as the interpolation of the previous and next column
    imageYCbCr(1:2:nrow-1,2:2:ncol-2,2:3) = (double(imageYCbCr(1:2:nrow-1,1:2:ncol-3,2:3))+ double(imageYCbCr(1:2:nrow-1,3:2:ncol-1,2:3)))/2;
    % the last even column doesn't have a next one so it's just equal to its previous column
    imageYCbCr(1:2:nrow-1,ncol,2:3) = imageYCbCr(1:2:nrow-1,ncol-1,2:3);

    % compute each even row of Cb and Cr as the interpolation of the previous and next row
    imageYCbCr(2:2:nrow-2,:,2:3) = (double(imageYCbCr(1:2:nrow-3,:,2:3)) + double(imageYCbCr(3:2:nrow-1,:,2:3)))/2;
    % the last even row doesn't have a next one so it's just equal to its previous row
    imageYCbCr(nrow,:,2:3) = imageYCbCr(nrow-1,:,2:3);  
else
    error('Unsupported subsampling factors')
end

% the transformation matrix for the YCbCr to RGB conversion
T = [1 0 1.402;...
     1 -0.344136 -0.714136;...
     1 1.772 0];

% the transformation offset for the YCbCr to RGB conversion
offset = T * [0;128;128];

% initialize the uint8 output RGB image with the same size as the input YCbCr image
imageRGB = zeros(size(imageYCbCr), 'uint8');


% every component p of the RGB image equals to T(p,1)*Y + T(p,2)*Cb + T(p,3)*Cr - offset(p)
for p = 1:3
    imageRGB(:,:,p) = imlincomb(T(p,1),imageYCbCr(:,:,1),T(p,2),imageYCbCr(:,:,2), ...
        T(p,3),imageYCbCr(:,:,3),-offset(p));
end


end


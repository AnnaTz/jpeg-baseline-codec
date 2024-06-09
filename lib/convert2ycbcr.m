function [imageY,imageCb,imageCr] = convert2ycbcr(imageRGB,subimg)

image_size = size(imageRGB);    % size of the RGB image
nrow = image_size(1);           % number of RGB rows
ncol = image_size(2);           % number of RGB columns


% the transformation matrix for the RGB to YCbCr conversion
T = [0.299 0.587 0.114;...
     -0.168736 -0.331264 0.5;...
     0.5 -0.418688 -0.081312];
          
% the transformation offset for the RGB to YCbCr conversion
offset = [0; 128; 128];


% initialize the output YCbCr image with the same size and type as the input RGB image 
imageYCbCr = zeros(size(imageRGB), 'like', imageRGB);


% every component p of the YCbCr image equals to T(p,1)*R + T(p,2)*G + T(p,3)*B + offset(p)
for p = 1:3
    imageYCbCr(:,:,p) = imlincomb(T(p,1),imageRGB(:,:,1),T(p,2),imageRGB(:,:,2), ...
        T(p,3),imageRGB(:,:,3),offset(p));
end


if isequal(subimg,[4 4 4])                % if the subsampling factors are [4 4 4]
    imageY(:,:) = imageYCbCr(:,:,1);      % keep all the rows and columns of the Y component
    imageCb(:,:) = imageYCbCr(:,:,2);     % keep all the rows and columns of the Cb component
    imageCr(:,:) = imageYCbCr(:,:,3);     % keep all the rows and columns of the Cr component
    
elseif isequal(subimg,[4 2 2])                  % if the subsampling factors are [4 2 2]
    imageY(:,:) = imageYCbCr(:,:,1);            % keep all the rows and columns of the Y component
    imageCb(:,:) = imageYCbCr(:,1:2:ncol,2);    % keep only the odd columns of the Cb component
    imageCr(:,:) = imageYCbCr(:,1:2:ncol,3);    % keep only the odd columns of the Cr component
      
elseif isequal(subimg,[4 2 0])                          % if the subsampling factors are [4 2 0]
    imageY(:,:) = imageYCbCr(:,:,1);                    % keep all the rows and columns of the Y component
    imageCb(:,:) = imageYCbCr(1:2:nrow,1:2:ncol,2);     % keep only the odd rows and columns of the Cb component
    imageCr(:,:) = imageYCbCr(1:2:nrow,1:2:ncol,3);     % keep only the odd rows and columns of the Cr component
else
    error('Unsupported subsampling factors')
end

end


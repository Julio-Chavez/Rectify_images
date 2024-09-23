function [C, RC] = merge_images(A,RA,B,RB,varargin)
% Adapted from imfuse function

AisRGB = size(A,3) > 1;
BisRGB = size(B,3) > 1;

[A,B,A_mask,B_mask,RC] = calculateOverlayImages(A,RA,B,RB);

if isempty(varargin)
    scaling = 'none';
    method = 'mergeMin';
elseif length(varargin) == 1
    scaling = varargin{1};
    method = 'mergeMin';
elseif length(varargin) == 2
    scaling = varargin{1};
    method = varargin{2};
else
    error('Incorrect number of imputs.')
end

switch lower(method)
    case 'mergeavg'
        C = ImgMergeAvg(A, AisRGB, B, BisRGB, scaling);
    case 'mergemax'
        C = ImgMergeMax(A, AisRGB, B, BisRGB, scaling);
    case 'mergemin'
        C = ImgMergeMin(A, AisRGB, B, BisRGB, scaling);
    case 'blend'
        C = ImgBlend(A, A_mask, AisRGB, B, B_mask, BisRGB, scaling);
    otherwise
        error('Unidentified method.')
end

% C = mergeImgBlend(A, A_mask, AisRGB, B, B_mask, BisRGB, scaling);
% C = mergeImg(A, AisRGB, B, BisRGB, scaling);

end



function [A_padded,B_padded,A_mask,B_mask,R_output] = calculateOverlayImages(A,RA,B,RB)

% First calculate output referencing object. World limits are minimum
% bounding box that contains world limits of both images. Resolution is the
% minimum resolution in each dimension. We don't want to down sample either
% image.
outputWorldLimitsX = [min(RA.XWorldLimits(1),RB.XWorldLimits(1)),...
                      max(RA.XWorldLimits(2),RB.XWorldLimits(2))];
                  
outputWorldLimitsY = [min(RA.YWorldLimits(1),RB.YWorldLimits(1)),...
                      max(RA.YWorldLimits(2),RB.YWorldLimits(2))];                 
                  
goalResolutionX = min(RA.PixelExtentInWorldX,RB.PixelExtentInWorldX);
goalResolutionY = min(RA.PixelExtentInWorldY,RB.PixelExtentInWorldY);

widthOutputRaster  = ceil(diff(outputWorldLimitsX) / goalResolutionX);
heightOutputRaster = ceil(diff(outputWorldLimitsY) / goalResolutionY);

R_output = imref2d([heightOutputRaster, widthOutputRaster]);
R_output.XWorldLimits = outputWorldLimitsX;
R_output.YWorldLimits = outputWorldLimitsY;

fillVal = 0;
A_padded = resampleImageToNewSpatialRef(A,RA,R_output,'bilinear',fillVal);
B_padded = resampleImageToNewSpatialRef(B,RB,R_output,'bilinear',fillVal);

[outputIntrinsicX,outputIntrinsicY] = meshgrid(1:R_output.ImageSize(2),1:R_output.ImageSize(1));
[xWorldOverlayLoc,yWorldOverlayLoc] = intrinsicToWorld(R_output,outputIntrinsicX,outputIntrinsicY);
A_mask = contains(RA,xWorldOverlayLoc,yWorldOverlayLoc);
B_mask = contains(RB,xWorldOverlayLoc,yWorldOverlayLoc);

end







function [B,RB] = resampleImageToNewSpatialRef(A,RA,RB,interpMethod,fillValue)
%resampleImageToNewSpatialRef Resample spatially referenced image to new grid.
%
%   [B,RB,MASK] = resampleImageToNewSpatialRef(A,RA,RB,FILLVALUE,METHOD)
%   resamples the spatially referenced image A,RA to a new grid whose world
%   extent and resolution is defined by RB. During resampling, the output
%   image B is assigned the value fillValue. METHOD defines the
%   interpolation method used during resampling, valid values are
%   'nearest','bilinear', and 'bicubic'. The output B,RB is an output
%   spatially referenced image. MASK is a logical image which is false
%   where fillValues were assigned and true where values of A were used to
%   define B.

% Copyright 2012 The MathWorks, Inc.

[bIntrinsicX,bIntrinsicY] = meshgrid(1:RB.ImageSize(2),1:RB.ImageSize(1));

[xWorldOverlayLoc,yWorldOverlayLoc] = RB.intrinsicToWorld(bIntrinsicX,bIntrinsicY);

% Convert these locations to intrinsic coordinates to be used in
% interpolation in image A.
[xIntrinsicLoc,yIntrinsicLoc] = RA.worldToIntrinsic(xWorldOverlayLoc,yWorldOverlayLoc);

% Resample to form output image
if isa(A,'double')
    B = interp2(A,xIntrinsicLoc,yIntrinsicLoc,interpMethod,fillValue);
else
    B = interp2(single(A),single(xIntrinsicLoc),single(yIntrinsicLoc),...
            interpMethod,fillValue);
end

% Move B back to original type
B = cast(B,class(A));
end


function result = ImgBlend(A, A_mask, AisRGB, B, B_mask, BisRGB, scaling)
        % Creates a transparent overlay image
        
        % make the images similar in size and type
        [A,B] = makeSimilar(A,B,AisRGB,BisRGB,scaling);
        
        % compute regions of overlap
        onlyA = A_mask & ~B_mask;
        onlyB = ~A_mask & B_mask;
        bothAandB = A_mask & B_mask;
        % AorB = A_mask | B_mask;
        
        % weight each image equally
        weight1 = 1;
        weight2 = 1;
        
        % allocate result image
        result = zeros([size(A,1) size(A,2) size(A,3)], class(A));
        
        % for each color band, compute blended output band
        for i = 1:size(A,3)
            a = A(:,:,i);
            b = B(:,:,i);
            r = result(:,:,i);
            r(onlyA) = a(onlyA);
            r(onlyB) = b(onlyB);
            r(bothAandB) = uint8( weight1 .* single(a(bothAandB)) + weight2 .* single(b(bothAandB)));
            result(:,:,i) = r;
        end
        
end % mergeImgBlend

function result = ImgMergeAvg(A, AisRGB, B, BisRGB, scaling)
        % make the images similar in size and type
        [A,B] = makeSimilar(A,B,AisRGB,BisRGB,scaling);
        
        % allocate result image
        r = zeros([size(A,1) size(A,2) 2], class(A));
        r(:,:,1) = A;
        r(:,:,2) = B;
        
        % average both images
        result = sum(r, 3) ./ sum(r~=0, 3);
        result = uint8(result);
end % mergeImgAvg

function result = ImgMergeMax(A, AisRGB, B, BisRGB, scaling)
        % make the images similar in size and type
        [A,B] = makeSimilar(A,B,AisRGB,BisRGB,scaling);
        
        % allocate result image
        r = zeros([size(A,1) size(A,2) 2], class(A));
        r(:,:,1) = A;
        r(:,:,2) = B;
        
        % get max intensity of both images
        result = max(r, [], 3);
        result = uint8(result);
end % mergeImgMax

function result = ImgMergeMin(A, AisRGB, B, BisRGB, scaling)
        % make the images similar in size and type
        [A,B] = makeSimilar(A,B,AisRGB,BisRGB,scaling);
        
        % allocate result image
        r = zeros([size(A,1) size(A,2) 2], class(A));
        r(:,:,1) = A;
        r(:,:,2) = B;
        
        % get min intensity of both images
        r(r == 0) = 255; % convert zeros to 255 to use min func properly
        result = min(r, [], 3);
        result(result == 255) = 0;
end % mergeImgMin



function image_data = scaleGrayscaleImage(image_data)

if (islogical(image_data))
    return
end
% convert to floating point
image_data = single(image_data);
minData = min(image_data(:));
maxData = max(image_data(:));

if (minData == maxData)
    return
end

% Scale to range [0 1]
image_data = (image_data - minData)/(maxData - minData);

end % scaleGrayscaleImage

function [A, B] = scaleTwoGrayscaleImages(A,B)
% takes 2 same size images and scales them as a single dataset

% convert to floating point composite image
image_data = [single(A) single(B)];

% scale as a single image
image_data = scaleGrayscaleImage(image_data);

% split images
A = image_data(:,1:size(A,2),:);
B = image_data(:,size(A,2)+1:end,:);

end % scaleTwoGrayscaleImages

function [A,B] = makeSimilar(A,B,AisRGB,BisRGB,scaling)

% scale the images
if ~AisRGB && ~BisRGB
    % both are grayscale
    switch lower(scaling)
        case 'none'
        case 'joint'
            [A,B] = scaleTwoGrayscaleImages(A,B);
        case 'independent'
            A = scaleGrayscaleImage(A);
            B = scaleGrayscaleImage(B);
    end
    
elseif AisRGB && BisRGB
    % both are RGB
    
elseif AisRGB && ~BisRGB
    % convert B to pseudo RGB "gray" image
    if strcmpi(scaling,'None')
        B = repmat(im2uint8(B),[1 1 3]);
    else
        % Scale the sole grayscale image alone
        B = scaleGrayscaleImage(B);
        B = repmat(im2uint8(B),[1 1 3]);
    end
    
elseif ~AisRGB && BisRGB
    % convert A to pseudo RGB "gray" image
    if strcmpi(scaling,'None')
        A = repmat(im2uint8(A),[1 1 3]);
    else
        % Scale the sole grayscale image alone
        A = scaleGrayscaleImage(A);
        A = repmat(im2uint8(A),[1 1 3]);
    end
end

% convert to uint8
A = im2uint8(A);
B = im2uint8(B);

end % makeSimilar
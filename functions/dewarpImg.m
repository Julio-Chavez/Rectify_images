function [dewarpedImg, xl, yl] = dewarpImg(image, varargin)
% A better algorithm to dewearp images
% Assuming image is the raw image in image coordinates 
% image = calImg{1}; 
% rectify_quad{1} = calibrate_camera(I{1}, W{1}, 2);
if length(varargin) == 1
    rectify_quad = varargin{1};
elseif length(varargin) == 2
    % Obtain rectification function based on provided calibration points
    I = varargin{1};
    W = varargin{2};
    rectify_quad = calibrate_camera(I, W, 2);
else
    error('Usage: [dewarpedImg, xl, yl] = dewarpImg(image, rectify_quad) or [dewarpedImg, xl, yl] = dewarpImg(image, I, W).')
end


% Get the size of the image
[imageHeight, imageWidth, numChannels] = size(image);

% Get grid of image
[X, Y] = meshgrid(1:imageWidth, 1:imageHeight);

% Apply rectification to all pixel grids 
imagePoints2 = rectify_quad([X(:), Y(:)]); 

% Reshape rectified image back to original dimensions 
rectified_image_x = reshape(imagePoints2(:,1), imageHeight, imageWidth, []);
rectified_image_y = reshape(imagePoints2(:,2), imageHeight, imageWidth, []);

XWorldLimits = [min(rectified_image_x(:)), max(rectified_image_x(:))];
YWorldLimits = [min(rectified_image_y(:)), max(rectified_image_y(:))];

% Flatten the image and grid coordinates for scattered interpolation
imagePixels = reshape(image, [], numChannels);

% Create a grid for the dewarped image
xl = linspace(XWorldLimits(1), XWorldLimits(2), imageWidth);
yl = linspace(YWorldLimits(1), YWorldLimits(2), imageHeight);
[Xl, Yl] = meshgrid(xl, yl);

% Interpolate each channel separately
dewarpedImg = zeros([imageHeight, imageWidth, numChannels], class(image));

for c = 1:numChannels
    F = scatteredInterpolant(imagePoints2(:,1), imagePoints2(:,2), double(imagePixels(:,c)), 'linear', 'none');
    dewarpedImg(:,:,c) = F(Xl, Yl);
end

% Convert the dewarped image back to the original class
% dewarpedImg = cast(dewarpedImg, class(image));
end





function rectify = calibrate_camera(I,W,order)
% calculate the transformation function to convert image coordinates to
% world (physical) coordinates, including calibration (converting pixels to
% meters), undistortion, and rectification
%
% rectify: function handle to map image coordinates to world coordinates 
% I: set of known calibration points in image coordinates (n x 2 vector) [px]
% W: set of known calibration points in world coordinates (n x 2 vector) [m]
% order: 1 for linear transformation (corrects for camera viewing angle but not lens distortion), 
%        2 for quadratic transformation (corrects for camera viewing angle and lens distortion)
%
% references: Fujita et al 1998 (Water Res.), Creutin et al 2003 (J. Hydrol.)
%
% to use function handle: points_m = rectify(points_px)
% points_px: set of points in image coordinates (m x 2 vector) [px]
% points_m: set of points converted to world coordinates (m x 2 vector) [px]


% find transformation coefficients

if order == 2
    A = [I.^2, I, ones(size(I,1),1), -I(:,1).^2.*W(:,1), -I(:,2).^2.*W(:,1), -I(:,1).*W(:,1), -I(:,2).*W(:,1), zeros(size(I,1),5);
        zeros(size(I,1),5), -I(:,1).^2.*W(:,2), -I(:,2).^2.*W(:,2), -I(:,1).*W(:,2), -I(:,2).*W(:,2), I.^2, I, ones(size(I,1),1)];
else
    A = [I, ones(size(I,1),1), -I(:,1).*W(:,1), -I(:,2).*W(:,1), zeros(size(I,1),3);
        zeros(size(I,1),3), -I(:,1).*W(:,2), -I(:,2).*W(:,2), I, ones(size(I,1),1)];
end
Z = [W(:,1); W(:,2)];
B = (A'*A)^-1*A'*Z;
    

% function to map image coords to world coords

if order == 2
    rectify = @(I) [(B(1)*I(:,1).^2 + B(2)*I(:,2).^2 + B(3)*I(:,1) + B(4)*I(:,2) + B(5))./ ...
        (B(6)*I(:,1).^2 + B(7)*I(:,2).^2 + B(8)*I(:,1) + B(9)*I(:,2) + 1), ...
        (B(10)*I(:,1).^2 + B(11)*I(:,2).^2 + B(12)*I(:,1) + B(13)*I(:,2) + B(14))./ ...
        (B(6)*I(:,1).^2 + B(7)*I(:,2).^2 + B(8)*I(:,1) + B(9)*I(:,2) + 1)];

else
    rectify = @(I) [(B(1)*I(:,1) + B(2)*I(:,2) + B(3))./ ...
        (B(4)*I(:,1) + B(5)*I(:,2) + 1), ...
        (B(6)*I(:,1) + B(7)*I(:,2) + B(8))./ ...
        (B(4)*I(:,1) + B(5)*I(:,2) + 1)];
end

end
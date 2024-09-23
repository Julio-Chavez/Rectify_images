function [croppedImage, modifiedRef] = cropImageWithRef(image, ref2d, pixelDecrease, orientation)
    % Get the image size
    imageSize = size(image);
    
    % Initialize the cropping limits
    xStart = 1;
    xEnd = imageSize(2);
    yStart = 1;
    yEnd = imageSize(1);
    
    % Modify the cropping limits based on the orientation
    switch lower(orientation)
        case 'left'
            xStart = xStart + pixelDecrease;
        case 'right'
            xEnd = xEnd - pixelDecrease;
        case 'top'
            yStart = yStart + pixelDecrease;
        case 'bottom'
            yEnd = yEnd - pixelDecrease;
        otherwise
            error('Invalid orientation. Use "left", "right", "top", or "bottom".');
    end
    
    % Crop the image
    croppedImage = image(yStart:yEnd, xStart:xEnd, :);
    
    % Modify the imref2d structure
    modifiedRef = ref2d;
    
    % Update XWorldLimits and YWorldLimits based on the cropping
    if any(strcmpi(orientation, {'left', 'right'}))
        modifiedRef.XWorldLimits = modifiedRef.XWorldLimits + ...
            [xStart - 1, xEnd - imageSize(2)] * modifiedRef.PixelExtentInWorldX;
        modifiedRef.ImageSize(2) = size(croppedImage, 2);
    end
    
    if any(strcmpi(orientation, {'top', 'bottom'}))
        modifiedRef.YWorldLimits = modifiedRef.YWorldLimits + ...
            [yStart - 1, yEnd - imageSize(1)] * modifiedRef.PixelExtentInWorldY;
        modifiedRef.ImageSize(1) = size(croppedImage, 1);
    end
end

function [calImg, I, W, rectify_quad, cams] = get_calibration_params(locationRunParams, locationCalParams, locationCalImg)

plot_calibration_on = false;

warning off
run_params = readtable(locationRunParams);
cal_params = readtable(locationCalParams);
warning on

% display the current experimental run
% fprintf('\nwindspeed = %2.f m/s, particle type = %s\n', run_params.WindSpeed_m_s(n), run_params.ParticleType{n});

% image parameters 
cams = cell2mat(cal_params.Cam)';

% preallocate image rectification function handles for each camera view
rectify_quad = cell(length(cams),1);

xyimg = cell(4,1); calIdx_rect = xyimg; calImg = xyimg;
 imagePoints2 = xyimg; W = xyimg; I = xyimg;

for cam = 1:length(cams)   
    
    cam_left = (cam <= 2); % 'true' for the left two cameras, 'false' for right two cameras

    % load background and calibration images
    bkgd = cam_imread(fullfile(locationCalImg, ['Cam' cams(cam) '-bkgd.tif']), cam_left);
    cal = cam_imread(fullfile(locationCalImg, ['Cam' cams(cam) '-cal.tif']), cam_left);
    mask = cam_imread(fullfile(locationCalImg, ['Cam' cams(cam) '-mask.tif']), cam_left);

    % background image must be "double" type
    bkgd = double(bkgd);

    % image dimensions
    % img_ix = size(bkgd,2); img_iy = size(bkgd,1);

    % CALIBRATION: GET MAPPING FUNCTION FROM IMAGE TO WORLD COORDINATES
    
    % subtract background from calibration image
    cal = double(cal) - bkgd; 

    % shift image intensities so that all are positive
    cal = uint8(cal - min(cal(:)));

    % invert image intensities
    cal = 255 - cal;

    calImg{cam} = cal;
    
    % binarize calibration image (adaptative binarization method)
    B = imbinarize(cal,'adaptive','Sensitivity', cal_params.calAdBinSens(cam));   

    % remove false positive pixels by eroding and dilating twice
    B = imerode(B,[0 1 0; 1 1 1; 0 1 0]); 
    B = imerode(B,[0 1 0; 1 1 1; 0 1 0]); 
    B = imdilate(B,[0 1 0; 1 1 1; 0 1 0]); 
    B = imdilate(B,[0 1 0; 1 1 1; 0 1 0]); 

    % number of rows and columns of dots in the calibration grid 
    n_rows_cal = cal_params.nRowsCal(cam); 
    n_cols_cal = cal_params.nColsCal(cam);

    % mask top and edges of calibration image (to remove dots that are not 
    % part of the n_rows_cal x n_cols_cal grid)
    B = B.*logical(mask);

    % view calibration images
    if plot_calibration_on
        figure; subplot(161); pcolor_img(cal); title('bkgd sub and inverted')
        subplot(162); pcolor_img(B); title('masked and binarized')
    end
    
    % detect calibration dots
    CC = bwconncomp(B);
    cal_dots = regionprops('table', CC, cal, 'Centroid', 'Area', 'EquivDiameter');

    % remove false detections based on dot size
    idx = cal_dots.Area > cal_params.calAreaThres1_px(cam) & cal_dots.Area < cal_params.calAreaThres2_px(cam);
    cal_dots = cal_dots(idx,:);    

    % view detected dots
    if plot_calibration_on
        subplot(163); pcolor_img(cal); hold on
        viscircles(cal_dots.Centroid,cal_dots.EquivDiameter/2); title('dots detected')
    end
    
    % known coords of the dots in physical space (world coordinates):
    % generate grid
    W{cam} = generateCircleGridPoints([n_rows_cal, n_cols_cal], cal_params.spacing_m(cam), "PatternType","symmetric") + ...
        [cal_params.worldOffset_x(cam), cal_params.worldOffset_y(cam)]*cal_params.spacing_m(cam);

    % corresponding coordinates of detected dots in the image (image coordinates): 
    % sort dots into ascending order from lower left corner of image by
    % separating the dots by x-coordinate into vertical bins
    [~,top_row] = sort(cal_dots.Centroid(:,2),'descend'); 
    top_row = top_row(1:n_cols_cal);  % top row of dots
    bin_lim = cal_dots.Centroid(top_row,1); 
    bin_lim = sort(bin_lim,'ascend');
    bin_lim = [bin_lim - diff(bin_lim(1:2))/2; inf];  % bin limits
    
    I{cam} = zeros(n_rows_cal*n_cols_cal,2); 
    for j = 1:n_cols_cal    
        if plot_calibration_on; line([bin_lim(j) bin_lim(j)],[0 img_iy]); end  % plot bin limits
        cal_dots0 = cal_dots(cal_dots.Centroid(:,1) > bin_lim(j) & cal_dots.Centroid(:,1) < bin_lim(j+1),:);  % dots in a vertical bin
        [~,sort_idx] = sort(cal_dots0.Centroid(:,2),'ascend');  % sort by y-coord
        I{cam}(n_rows_cal*(j-1)+1 : n_rows_cal*j,:) = cal_dots0.Centroid(sort_idx,:);  % image points
    end
    if plot_calibration_on  % plot detected image points
        point_colors = jet(size(I{cam},1));
        subplot(164); scatter(I{cam}(:,1),I{cam}(:,2),15,point_colors,"filled"); grid on; axis tight equal; title('detected points')
        subplot(165); scatter(W{cam}(:,1),W{cam}(:,2),15,point_colors,"filled"); grid on; axis tight equal; title('known coords')
    end
    
    % CALIBRATION: compute the image rectification function for this camera
    % from the world coords and image coords of the calibration dots (using
    % a quadratic transformation: order = 2)
    trans_order = 2;
    rectify_quad{cam} = calibrate_camera(I{cam},W{cam},trans_order);    
    imagePoints2{cam} = rectify_quad{cam}(I{cam});
    % view calibration to confirm that calibration dots are mapped to the correct world coords
    if plot_calibration_on
        subplot(166); scatter(imagePoints2{cam}(:,1),imagePoints2{cam}(:,2),15,point_colors,"filled"); 
        grid on; axis tight equal; title('corrected points')
    end

    [ximg, yimg] = meshgrid((1:size(cal,2)), (1:size(cal,1)));
    xyimg{cam} = [ximg(:), yimg(:)];
    calIdx_rect{cam} = rectify_quad{cam}(xyimg{cam});

end

end
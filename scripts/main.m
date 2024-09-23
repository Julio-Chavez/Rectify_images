clc
clear


locationCalParams = '../data/raw/tif_files/cal_parameters_demo.xlsx';
locationRunParams = '../data/raw/tif_files/run_parameters_demo.xlsx';
locationCalImg = '../data/raw/tif_files/calibration_and_bkgd';
locationRunImg = '../data/raw/tif_files/run1';

locationCal = '/Users/juliochavez/Desktop/ShadowTracking/data/data_demo1/calibration_and_bkgd/';
% locationCalImg = '/Users/juliochavez/Desktop/ShadowTracking/data/data_demo/calibration_and_bkgd';
locationImg = '/Users/juliochavez/Desktop/shadowtrackingLuci/data_demo/run1';
locationSave = './';
folderSave = 'run1';

fs = 30; % sampling frequency

% Define the input directory and the output directory
outputDir = fullfile(locationSave, folderSave); % Directory to save the output images

% Check if the output directory exists, if not, create it
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% get calibration images and rectification function
[calImg, I, W, rectify_quad, cams] = get_calibration_params(locationRunParams, locationCalParams, locationCalImg);

%% load and dewarp calibration images
locationSaveCalDewarped = '../data/output/tif_files/cal_dewarped';

dewarpedCal = cell(4,2); imgRef = cell(4,2);
xlCal = dewarpedCal; ylCal = dewarpedCal;

for cam = 1:length(cams)
    % Put camera name in first column
    dewarpedCal{cam,1} = cams(cam);
    xlCal{cam,1} = cams(cam);
    ylCal{cam,1} = cams(cam);
    imgRef{cam,1} = cams(cam);

    % dewarp calibration images
    [dewarpedCal{cam,2}, xlCal{cam,2}, ylCal{cam,2}] = dewarpImg(calImg{cam}, rectify_quad{cam});

    % reference the calibration images
    imgRef{cam,2} = imref2d(size(dewarpedCal{cam,2}), [min(xlCal{cam,2}) max(xlCal{cam,2})], [min(ylCal{cam,2}) max(ylCal{cam,2})]);
    imwrite( dewarpedCal{cam,2}, fullfile(locationSaveCalDewarped, ['Cam' cams(cam) '-cal.tif']) )
end

save( fullfile(locationSaveCalDewarped, 'calRef.mat'), 'xlCal', 'ylCal', 'imgRef' )

%% load and dewarp bkgd
dewarpedBkgd = cell(4,2);

fprintf('Dewarping background images...')
tt = tic;
parfor cam = 1:length(cams) 
    cam_left = (cam <= 2); % 'true' for the left two cameras, 'false' for right two cameras

    dewarpedCal{cam,1} = cams(cam);
    img_name = ['Cam' cams(cam) '-bkgd.tif'];
    img = cam_imread(fullfile(locationCalImg, img_name), cam_left);

    % dewarp images
    [dewarpedBkgd{cam,2}, ~, ~] = dewarpImg(img, rectify_quad{cam});

    % save images
    imwrite( dewarpedBkgd{cam,2}, fullfile(locationSaveCalDewarped, ['Cam' cams(cam) '-bkgd.tif']) )
end
toc(tt)

%% plot dewarped Calibration and Background images
figure(1)
montage({dewarpedBkgd{1,2}, dewarpedBkgd{2,2},dewarpedBkgd{3,2},dewarpedBkgd{4,2}}, "Size",[1 4])
title('Background')

figure(2)
montage({dewarpedCal{1,2}, dewarpedCal{2,2},dewarpedCal{3,2},dewarpedCal{4,2}}, "Size",[1 4])
title('Calibration plate')

%% load and dewarp run images
locationSaveDewarped = '../data/output/tif_files/run1_dewarped';
frames=(1000:1002);

% n= 19; % for seq files
% [~,headerInfo] = loadseqframe(fullfile(locationRunImg,'CamA', ['CamA_', run_params.ID{n}, '.seq']),1);
% run_name = run_params.ID{n};

numFrames = numel(frames);
numCams = length(cams);

tt = tic;

% Flatten the parallelization so that each iteration works on a unique combination of frame and camera
parfor idx = 1:(numFrames * numCams)
    % Compute frame and camera indices
    frame = ceil(idx / numCams);
    cam = mod(idx - 1, numCams) + 1;
    cam_left = (cam <= 2); % 'true' for the left two cameras, 'false' for right two cameras

    % Perform dewarping
    fprintf('Dewarping image...')
    t = tic;
    % load tif image
    img = cam_imread_seq(fullfile(locationRunImg, ...
            ['Cam', cams(cam),'-', sprintf('%04d', frames(frame)), '.tif']), cam_left, []);
    % load seq image
    % img = cam_imread_seq(fullfile(locationRunImg,['Cam', cams(cam)], ...
    %         ['Cam', cams(cam),'_', run_name, '.seq']), cam_left, frames(frame));
    temp3 = dewarpImg(img, rectify_quad{cam});
    imwrite(temp3, fullfile(locationSaveDewarped, ['Cam', cams(cam),'-', sprintf('%04d', frames(frame) ), '.tif']) )
    toc(t)
end

toc(tt)

clear temp_cams temp_frames temp3 frame cam cam_left img tt

% save([locationSaveDewarped,'metadata_dewarped.mat'],'headerInfo', 'imgRefMerged')

%% ======================MERGE IMAGES================================

%% merge background images (to do cropping and translation)
locationSaveCalMerged = '../data/output/tif_files/cal_merged';

scaling = 'none'; method = 'mergeMin';
fprintf('Merging background images...')

cropPix = [0 0 0 0 0 0]; %[90, 478, 352, 114, 525, 112];
tranmm = [0 0 0 0 0 0 0 0]; %[0 0 -.009 -.00125 0 0 0 .0011];


tt = tic;

dewarpedBkgd_cropped = cell(4,1); imRefCropped = imgRef;
% crop images to get rid of grey bands at edge of images (make sure you
% have enough overlap)
% CamA
[dewarpedBkgd_cropped{1}, imRefCropped{1,2}] = cropImageWithRef(dewarpedBkgd{1,2}, imgRef{1,2},cropPix(1), 'right');
% CamC
[dewarpedBkgd_cropped{2}, imRefCropped{2,2}] = cropImageWithRef(dewarpedBkgd{2,2}, imgRef{2,2}, cropPix(2), 'left');
[dewarpedBkgd_cropped{2}, imRefCropped{2,2}] = cropImageWithRef(dewarpedBkgd_cropped{2}, imRefCropped{2,2}, cropPix(3), 'right');
% CamB
[dewarpedBkgd_cropped{3}, imRefCropped{3,2}] = cropImageWithRef(dewarpedBkgd{3,2}, imgRef{3,2}, cropPix(4), 'left');
[dewarpedBkgd_cropped{3}, imRefCropped{3,2}] = cropImageWithRef(dewarpedBkgd_cropped{3}, imRefCropped{3,2}, cropPix(5), 'right');
% CamD
[dewarpedBkgd_cropped{4}, imRefCropped{4,2}] = cropImageWithRef(dewarpedBkgd{4,2}, imgRef{4,2}, cropPix(6), 'left');



% correct the position of the images by translation
imgRefTransd = imRefCropped; % to avoid modifying this variable in each run
imgRefTransd{1,2}.XWorldLimits = imRefCropped{1,2}.XWorldLimits + tranmm(1);
imgRefTransd{1,2}.YWorldLimits = imRefCropped{1,2}.YWorldLimits + tranmm(2);
imgRefTransd{2,2}.XWorldLimits = imRefCropped{2,2}.XWorldLimits + tranmm(3);
imgRefTransd{2,2}.YWorldLimits = imRefCropped{2,2}.YWorldLimits + tranmm(4);
imgRefTransd{3,2}.XWorldLimits = imRefCropped{3,2}.XWorldLimits + tranmm(5);
imgRefTransd{3,2}.YWorldLimits = imRefCropped{3,2}.YWorldLimits + tranmm(6);
imgRefTransd{4,2}.XWorldLimits = imRefCropped{4,2}.XWorldLimits + tranmm(7);
imgRefTransd{4,2}.YWorldLimits = imRefCropped{4,2}.YWorldLimits + tranmm(8);


[bAC, RbAC] = merge_images(dewarpedBkgd_cropped{1},imgRefTransd{1,2},dewarpedBkgd_cropped{2},imgRefTransd{2,2}, scaling, method );
% RbAC.XWorldLimits = RbAC.XWorldLimits - .009;
% RbAC.YWorldLimits = RbAC.YWorldLimits - .00125;
[bACB, RbACB] = merge_images(bAC, RbAC,dewarpedBkgd_cropped{3},imgRefTransd{3,2}, scaling, method );
[Hb, mergedImgRef] = merge_images(bACB, RbACB,dewarpedBkgd_cropped{4},imgRefTransd{4,2}, scaling, method );
toc(tt)
Hb = flipud(Hb); % WASIRF specific
imwrite(Hb, fullfile(locationSaveCalMerged,'bkgd.tif'));
disp('All background images processed and saved.');

x = (1:mergedImgRef.ImageSize(2))*mergedImgRef.PixelExtentInWorldX;
y = -(1:mergedImgRef.ImageSize(1))*mergedImgRef.PixelExtentInWorldY;
imagesc(x,y,Hb)%(300:700, 1500:2300))
colormap(gray)
axis equal tight xy

%% Merge calibration plate images

scaling = 'none'; method = 'mergeMin';
fprintf('Merging calibration images...')
tt = tic;
[cAC, RcAC] = merge_images(dewarpedCal{1,2},imgRefTransd{1,2},dewarpedCal{2,2},imgRefTransd{2,2}, scaling, method );
[cACB, RcACB] = merge_images(cAC, RcAC,dewarpedCal{3,2},imgRefTransd{3,2}, scaling, method );
[Hc, RHc] = merge_images(cACB, RcACB,dewarpedCal{4,2},imgRefTransd{4,2}, scaling, method );
toc(tt)
Hc = flipud(Hc); % WASIRF specific
imwrite(Hc, fullfile(locationSaveCalMerged,'calPlate.tif'));
disp('All background images processed and saved.');

imshow(Hc)


%% merge run images
locationSaveMerged = '../data/output/tif_files/run1_merged';
frames=(1000:1002);
numFrames = numel(frames); numCams = length(cams);
fprintf(['Merging ' num2str(numFrames) ' frames from run ' num2str(n) '...\n'])
temp_frames = nan(numFrames,1);

% % Check if the output directory exists, if not, create it
% if ~exist(fullfile(saveloc, outputDir), 'dir')
%     mkdir(fullfile(saveloc, outputDir));
% end
scaling = 'none'; method = 'mergeMin';

tt = tic;
for ii =1:numFrames
    temp_frames(ii) = frames(ii);

    dewA = imread( fullfile(locationSaveDewarped,['CamA-' sprintf('%04d', frames(ii)) '.tif']) );
    dewC = imread( fullfile(locationSaveDewarped,['CamC-' sprintf('%04d', frames(ii)) '.tif']) );
    dewB = imread( fullfile(locationSaveDewarped,['CamB-' sprintf('%04d', frames(ii)) '.tif']) );
    dewD = imread( fullfile(locationSaveDewarped,['CamD-' sprintf('%04d', frames(ii)) '.tif']) );

    dewarped_cropped = cell(4,1); imRefCropped = imgRef;
    % crop images to get rid of grey bands at edges
    % CamA
    [dewarped_cropped{1}, ~] = cropImageWithRef(dewA, imgRef{1,2}, cropPix(1), 'right');
    % CamB
    [dewarped_cropped{2}, ~] = cropImageWithRef(dewC, imgRef{2,2}, cropPix(2), 'left');
    [dewarped_cropped{2}, ~] = cropImageWithRef(dewarped_cropped{2}, imRefCropped{2,2}, cropPix(3), 'right');
    % CamC
    [dewarped_cropped{3}, ~] = cropImageWithRef(dewB, imgRef{3,2}, cropPix(4), 'left');
    [dewarped_cropped{3}, ~] = cropImageWithRef(dewarped_cropped{3}, imRefCropped{3,2}, cropPix(5), 'right');
    % CamD
    [dewarped_cropped{4}, ~] = cropImageWithRef(dewD, imgRef{4,2}, cropPix(6), 'left');

    %merge
    [AC, RAC] = merge_images( dewarped_cropped{1}, imgRefTransd{1,2},dewarped_cropped{2}, imgRefTransd{2,2}, scaling, method );
    [ACB, RACB] = merge_images( AC, RAC, dewarped_cropped{3}, imgRefTransd{3,2}, scaling, method );
    [ACBD, imgRefMerged] = merge_images( ACB, RACB, dewarped_cropped{4}, imgRefTransd{4,2}, scaling, method );

    % flip final image (WASIRF specific)
    ACBD = flipud(ACBD);

    % save image as tif 
    imgName = ['img_' sprintf('%04d', frames(ii)) '.tif']; % Create the file names
    imwrite(ACBD, fullfile(locationSaveMerged, imgName)); % Save the images as TIFF files
    
end
toc(tt)

% see merged images
% ii = 1;
figure(1)
imagesc(x,y,ACBD)%(1500:2000,3000:4000))
colormap(gray)
axis equal tight xy
drawnow

















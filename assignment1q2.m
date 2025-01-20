%MCEN90018 - Advanced Fluid Dynamics - Assignment 1 - Q1
% This script answers Question 1 in the Assignment 1 sheet.

% Author: Francesco Mascadri
% Contact: fmascadri@student.unimelb.edu.au
% April 2022

%% Clear workspace
clear all;
clc;
close all;

%% Plotting flag
PLOTTING = true; % change to true to generate plots

%% Tell script where to find assignment images
imagedirectory = './data/calibration_images/';

%% Number of calibration dots
numDotsRows = 17;
numDotsCols = 21;

%% Setup template
dsize = 15;
sigi = dsize./2.5;
sigj = dsize./2.5;
[it,jt] = meshgrid([-dsize:1:dsize],[-dsize:1:dsize]);
template = 255*exp(-(((it).^2)./(2*sigi.^2) + ((jt).^2)./(2*sigj.^2)));

%plotting for template image
if PLOTTING
    figure
    imagesc(template)
    set(gca, 'ydir', 'normal');
    colormap("gray");
    daspect([1 1 1]);
end

fprintf('Template instantiated.\n')
%% Question 2a
%Read images
leftImages = dir([imagedirectory, 'cal_image_left_*.tiff']);
rightImages = dir([imagedirectory, 'cal_image_right_*.tiff']);

%Initalise array for storing dot data
leftImagesXLocations = zeros(1,numDotsCols*numDotsRows,length(leftImages));
leftImagesYLocations = zeros(1,numDotsCols*numDotsRows,length(leftImages));
rightImagesXLocations = zeros(1,numDotsCols*numDotsRows,length(rightImages));
rightImagesYLocations = zeros(1,numDotsCols*numDotsRows,length(rightImages));

%Left images
for i = 1:length(leftImages)
    fprintf(['Processing image: ' leftImages(i).name '\n'])
    calL = imread(leftImages(i).name);
    calL = flipud(mean(calL,3));

    [xmax, ymax] = get_dot_locations(template,calL,numDotsCols,numDotsRows);

    [leftImagesXLocations(:,:,i), leftImagesYLocations(:,:,i)] = sort_matrices(xmax,ymax,numDotsRows,numDotsCols);    
end

fprintf('Q2a Left images complete.\n')

%Right images
for i = 1:length(rightImages)
    fprintf(['Processing image: ' rightImages(i).name '\n'])
    calR = imread(rightImages(i).name);
    calR = flipud(mean(calR,3));

    [xmax, ymax] = get_dot_locations(template,calR,numDotsCols,numDotsRows);

    [rightImagesXLocations(:,:,i), rightImagesYLocations(:,:,i)] = sort_matrices(xmax,ymax,numDotsRows,numDotsCols);    
end

fprintf('Q2a Right images complete.\n')

%% Plot results from all calibration images for Q2b
if PLOTTING
    figure
    ax1 = subplot(1,2,1);
    hold on
    for n=1:length(leftImages)
        scatter(ax1,leftImagesXLocations(:,:,n),leftImagesYLocations(:,:,n),7,'filled');
    end
    title(ax1,'Left image')
    xlabel(ax1,'x-distance (pixels)');
    ylabel(ax1,'y-distance (pixels)');
    daspect(ax1,[1 1 1]);
    axis(ax1,[0 4667 0 3500]);
    grid on
    
    ax2 = subplot(1,2,2);
    hold on
    for m=1:length(rightImages)
        scatter(ax2,rightImagesXLocations(:,:,m),rightImagesYLocations(:,:,m),7,'filled');
    end
    title(ax2,'Right image')
    xlabel(ax2,'x-distance (pixels)');
    ylabel(ax2,'y-distance (pixels)');
    daspect(ax2,[1 1 1]);
    axis(ax2,[0 4667 0 3500]);
    grid on
    
    legend('1900 mm', ...
           '1920 mm', ...
           '1940 mm', ...
           '1960 mm', ...
           '1980 mm', ...
           '2000 mm', ...
        'Location', 'best', ...
        'Orientation','horizontal')
    
    sgtitle('Q2a Calibration results')
    
    hold off
end

%% Part 2b

% Reshape data into appropriate format
calibration = [[leftImagesXLocations(:,:,1), leftImagesXLocations(:,:,2),...
                leftImagesXLocations(:,:,3), leftImagesXLocations(:,:,4),...
                leftImagesXLocations(:,:,5), leftImagesXLocations(:,:,6),];
               [leftImagesYLocations(:,:,1), leftImagesYLocations(:,:,2),...
                leftImagesYLocations(:,:,3), leftImagesYLocations(:,:,4),...
                leftImagesYLocations(:,:,5), leftImagesYLocations(:,:,6),];
               [rightImagesXLocations(:,:,1), rightImagesXLocations(:,:,2),...
                rightImagesXLocations(:,:,3), rightImagesXLocations(:,:,4),...
                rightImagesXLocations(:,:,5), rightImagesXLocations(:,:,6),];
               [rightImagesYLocations(:,:,1), rightImagesYLocations(:,:,2),...
                rightImagesYLocations(:,:,3), rightImagesYLocations(:,:,4),...
                rightImagesYLocations(:,:,5), rightImagesYLocations(:,:,6),];
              ];

% Column format required for nlinfit
calibration = calibration';

% Model function
modelfun = @(b,X)(b(1)+...
                 b(2).*X(:,1)+...
                 b(3).*X(:,2)+...
                 b(4).*X(:,3)+...
                 b(5).*X(:,4)+...
                 b(6).*X(:,1).*X(:,2)+...
                 b(7).*X(:,1).*X(:,3)+...
                 b(8).*X(:,1).*X(:,4)+...
                 b(9).*X(:,2).*X(:,3)+...
                 b(10).*X(:,2).*X(:,4)+...
                 b(11).*X(:,3).*X(:,4)+...
                 b(12).*X(:,1).^2+...
                 b(13).*X(:,2).^2+...
                 b(14).*X(:,3).^2+...
                 b(15).*X(:,4).^2);

% Initial guesses for coefficients 
beta0_x = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
beta0_y = [1000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
beta0_z = [2000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];

% True location data setup in the way that matches the sort direction of
% the image x y z data.
true_x = repmat(repelem(-500:50:500,1,numDotsRows),1,6)';
true_y = repmat(800:-50:0,1,numDotsCols*6)';
true_z = repelem(1900:20:2000, 1, numDotsRows*numDotsCols)';

% Fitted coefficients
xcoeff = nlinfit(calibration, true_x, modelfun, beta0_x);
ycoeff = nlinfit(calibration, true_y, modelfun, beta0_y);
zcoeff = nlinfit(calibration, true_z, modelfun, beta0_z);

% Compare with known calibration image
sample_image = [leftImagesXLocations(:,:,1);
                leftImagesYLocations(:,:,1);
                rightImagesXLocations(:,:,1);
                rightImagesYLocations(:,:,1)]';

x_sample = modelfun(xcoeff, sample_image);
y_sample = modelfun(ycoeff, sample_image);

%% Plot estimated vs true for report
if PLOTTING
    figure
    hold on
    scatter(x_sample, y_sample, 7, 'filled');
    scatter(true_x, true_y, 7, 'filled');
    legend("Estimated locations", "True locations")
    daspect([1 1 1])
    axis([-700 700 -200 1000])
    title('Estimated dot locations vs true calibration locations (z=1900mm)')
    xlabel("x-distance (mm)")
    ylabel("y-distance (mm)")
    grid on
    hold off
end

%% Part 2b) Test images
% Read images
image_left = imread('./data/test_images/test_left.tiff');
image_right = imread('./data/test_images/test_right.tiff');

image_left = mean(image_left, 3);
image_right = mean(image_right, 3);

% image_left = flip(image_left);
% image_right = flip(image_right);

% Calibration dot parameters
numDotsCols = 25;
numDotsRows = 17;

% Left image processing
[xmax, ymax] = get_dot_locations(template, image_left, numDotsCols, numDotsRows);
[leftImageXLocations, leftImageYLocations] = sort_matrices(xmax,ymax,numDotsRows,numDotsCols);

% Right image processing
[xmax, ymax] = get_dot_locations(template, image_right, numDotsCols, numDotsRows);
[rightImageXLocations, rightImageYLocations] = sort_matrices(xmax,ymax,numDotsRows,numDotsCols);

% Convert to physical space
image = [leftImageXLocations;...
         leftImageYLocations;...
         rightImageXLocations;...
         rightImageYLocations]';

x_output = modelfun(xcoeff,image);
y_output = modelfun(ycoeff,image);
z_output = modelfun(zcoeff,image);

% Reshape data for surface plot
xgrid = reshape(x_output, numDotsRows, numDotsCols);
ygrid = reshape(y_output, numDotsRows, numDotsCols);
zdist = reshape(z_output, numDotsRows, numDotsCols);

%% Plotting for Q2b
if true
    figure;
    s = surf(xgrid, zdist, ygrid);
    daspect([1 1 1]);
    xlabel("X (mm)")
    ylabel("Z (mm)")
    zlabel("Y (mm)")
    set(s,"CData",zdist);
    c = colorbar();
    c.Label.String = "z distance from camera plane (mm)";

    figure
    s = surf(xgrid, zdist, ygrid);
    xlabel("X (mm)")
    ylabel("Z (mm)")
    zlabel("Y (mm)")
    set(s,"CData",zdist);
    c = colorbar();
    c.Label.String = "z distance from camera plane (mm)";
    daspect([1 1 1]);
    axis([-500 500 0 2200 0 1000]);

end

%% Part 2c

% Read images
left2c = imread("./data/test_images/mystery_left_2022.tiff");
right2c = imread("./data/test_images/mystery_right_2022.tiff");

left2c = mean(left2c, 3);
right2c = mean(right2c, 3);

% % Apply gaussian filter to image
% left2c = imgaussfilt(left2c,2);
% right2c = imgaussfilt(right2c,2);

%left2c = flipud(left2c);
%right2c = flipud(right2c);

% Cross-correlation parameters
winsize = 64;
wsize = [winsize, winsize];

[xgrid_l, ygrid_l] = meshgrid((winsize+1):winsize:(4667-(winsize-1)),(winsize+1):winsize:(3500-(winsize-1)));

deltax_l = zeros(size(ygrid_l,1),size(xgrid_l,2));
deltay_l = zeros(size(ygrid_l,1),size(xgrid_l,2));

%Left image
hLeft = figure;
hLeft.OuterPosition = [200 150 1152 788];
hLeft.InnerPosition = [200 150 1136 695];
ax1 = gca;
hold on
imagesc(ax1, left2c);
colormap gray;
daspect([1 1 1]);
axis([0 4667 0 3500]);
title('Left image');
xlabel('x (pixels)');
ylabel('y (pixels)');

x_points = reshape(xgrid_l,1,[]);
y_points = reshape(ygrid_l,1,[]);
scatter(x_points, y_points, 'g');

xline(xgrid_l(1,:), 'LineStyle','-','Color','w','LineWidth',1.2);
yline(ygrid_l(:,1), 'LineStyle','-','Color','w','LineWidth',1.2);

[dpx_est_decimal, dpy_est_decimal] = calculateDisplacementSubpixel(left2c,right2c,wsize,xgrid_l,ygrid_l);
dpx_est_decimal(isnan(dpx_est_decimal)) = 0;
dpy_est_decimal(isnan(dpy_est_decimal)) = 0;
dpx_est = round(dpx_est_decimal);
dpy_est = round(dpy_est_decimal);

% Add the calculated displacements to the left grid to get the right grid
xgrid_r = xgrid_l + dpx_est;
ygrid_r = ygrid_l + dpy_est;

% Then sort the data into column form
[x_left, y_left] = sort_matrices(xgrid_l, ygrid_l, size(xgrid_l,1), size(xgrid_l,2));
[x_right, y_right] = sort_matrices(xgrid_r, ygrid_r, size(xgrid_l,1), size(xgrid_l,2));

% Then run it through the model
image_mystery = [x_left; y_left; x_right; y_right]';
x_mystery_data = modelfun(xcoeff, image_mystery);
y_mystery_data = modelfun(ycoeff, image_mystery);
z_mystery_data = modelfun(zcoeff, image_mystery);

%% Plot surface
% Reshape data
xgrid_final = reshape(x_mystery_data, size(xgrid_l,1), size(xgrid_l,2));
ygrid_final = reshape(y_mystery_data, size(xgrid_l,1), size(xgrid_l,2));
zdist_final = reshape(z_mystery_data, size(xgrid_l,1), size(xgrid_l,2));

figure
s = surf(xgrid_final, zdist_final, ygrid_final);
xlabel("X (mm)")
ylabel("Z (mm)")
zlabel("Y (mm)")
set(s,"CData",zdist_final);
c = colorbar();
c.Label.String = "z distance from camera plane (mm)";
daspect([1 1 1]);
axis([-800 800 1700 2200 -100 1000]);

%% Functions
function [xmax, ymax] = get_dot_locations(template, search, numDotsCols, numDotsRows)

    R = normxcorr2(template, search);

    ymax = zeros(1,numDotsCols*numDotsRows);
    xmax = zeros(1,numDotsCols*numDotsRows);
    i = 1;

    % find peaks in correlation matrix, extract peak, and exclude peak
    % until all peaks exhausted
    while find(R>0.8)
        [ymax(i), xmax(i)] = find(R==max(R(:)));
        dpx = xmax(i) - floor(size(template,2)/2);
        dpy = ymax(i) - floor(size(template,1)/2);
        R(dpy:dpy+size(template,1),dpx:dpx+size(template,2)) = -1;
        i = i + 1;
    end
end

function [final_x, final_y] = sort_matrices(xmax, ymax, num_rows, num_cols)
    % sort matrices
    [semi_sorted_x, idx] = sort(xmax, 'ascend');
    semi_sorted_y = ymax(idx);
    
    final_x = [];
    final_y = [];
 
    for j=0:num_cols-1
        [temp_y, idx] = sort(semi_sorted_y(...
            1+(j*num_rows):num_rows+(j*num_rows)), 'ascend');
        temp_x = semi_sorted_x(1+(j*num_rows):num_rows+(j*num_rows));
        temp_x = temp_x(idx);
        
        final_x = [final_x temp_x];
        final_y = [final_y temp_y];
    end
    
end

function [dpx, dpy] = calculateDisplacementSubpixel(imagea, imageb, wsize, xgrid, ygrid)

dpx = zeros(size(xgrid,1),size(xgrid,2));
dpy = zeros(size(ygrid,1),size(ygrid,2));

for i = 1:size(ygrid,1)
    m = ygrid(i,1);
    for j = 1:size(xgrid,2)
        n = xgrid(1,j);
        x1 = n - (wsize(1))/2;
        x2 = n + (wsize(1))/2 - 1;
        y1 = m - (wsize(2))/2;
        y2 = m + (wsize(2))/2 - 1;
        
        % normalization
        image1 = imagea(y1:y2,x1:x2);
        image2 = imageb(y1:y2,x1:x2);
        image1n = image1 - mean(mean(image1));
        sigma_1 = sqrt(mean(mean(image1n.^2)));
        image2n = image2 - mean(mean(image2));
        sigma_2 = sqrt(mean(mean(image2n.^2)));
        if sigma_1 == 0
            sigma_1 = 0.01;
        end
        if sigma_2 == 0
            sigma_2 = 0.01;
        end
        image1n = image1n./sigma_1;
        image2n = image2n./sigma_2;
        
        % cross correlation of image pairs
        C = xcorr2(image2n,image1n);
        [ymax, xmax] = find(C==(max(max(C))));
        
        ymax = ymax(1);
        xmax = xmax(1);
        
        if ymax == 1
            ymax = 2;
        end
        if xmax == 1
            xmax = 2;
        end
        if ymax == 2*wsize(1)-1
            ymax = 2*wsize(1)-2;
        end
        if xmax == 2*wsize(1)-1
            xmax = 2*wsize(1)-2;
        end

        % Gaussian sub pixel fitting
        y2 = C(ymax,xmax);
        y1 = C((ymax - 1), xmax);
        y3 = C((ymax + 1), xmax);
        
        x2 = C(ymax,xmax);
        x1 = C(ymax,(xmax-1));
        x3 = C(ymax,(xmax+1));
        
        deltay = (log(y1) - log(y3))/(log(y1)+log(y3)-2*log(y2))/2;
        deltax = (log(x1) - log(x3))/(log(x1)+log(x3)-2*log(x2))/2;
        
        dpx(i,j) = xmax + deltax - wsize(1);
        dpy(i,j) = ymax + deltay - wsize(2);
    end
end

end

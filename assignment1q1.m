%MCEN90018 - Advanced Fluid Dynamics - Assignment 1 - Q1
% This script answers Question 1 in the Assignment 1 sheet.

% Author: Francesco Mascadri
% Contact: fmascadri@student.unimelb.edu.au
% April 2022

%% Clear workspace
clear all;
clc;
close all;

%% Flags for generating plots
PLOTTING = true;

%% Tell script where to find assignment images
imagedirectory = './data/';

%% Question 1a
%Read images
left = imread([imagedirectory, 'mars/mars_L.jpg']);
right = imread([imagedirectory, 'mars/mars_R.jpg']);

left = double(left);
left = mean(left,3);
right = double(right);
right = mean(right,3);

left = flipud(left);
right = flipud(right);

%Interrogation window inputs
winsize = 32;
wsize = [winsize,winsize];
[xgrid, ygrid] = meshgrid((winsize+1):winsize:(1024-(winsize-1)),(winsize+1):winsize:(1024-(winsize-1)));

deltax = zeros(size(ygrid,2),size(xgrid,1));
deltay = zeros(size(ygrid,2),size(xgrid,1));

%Left image
hLeft = figure;
hLeft.OuterPosition = [200 150 1152 788];
hLeft.InnerPosition = [200 150 1136 695];
ax1 = subplot(2,2,1);
hold on
imagesc(ax1, left);
colormap gray;
daspect([1 1 1]);
axis([0 1024 0 1024]);
title('Left image');
xlabel('x (pixels)');
ylabel('y (pixels)');

x_points = reshape(xgrid,1,[]);
y_points = reshape(ygrid,1,[]);
scatter(x_points, y_points, 'g');

xline(xgrid(1,:), 'LineStyle','-','Color','w','LineWidth',1.2);
yline(ygrid(:,1), 'LineStyle','-','Color','w','LineWidth',1.2);

%live contour plot updating as loop runs
testfig = subplot(2,2,4);
testfig.Title.String = "Live updating contour plot";

template = subplot(2,2,2);
daspect(template, [1 1 1]);
template.Title.String = "Template image";

sgtitle({'Q1 Cross-correlation processing',['Window size = ' num2str(winsize) 'px']}) 

% peak R value
peak = zeros(size(xgrid));

for i=1:size(ygrid,2)
    intervalY = (winsize*i+1):(winsize*(i+1));
    for j=1:size(xgrid,2)

        figure(hLeft)
        intervalX = (winsize*j+1):(winsize*(j+1));

        line(ax1, [intervalX(1) intervalX(1);...
              intervalX(end) intervalX(1);...
              intervalX(1) intervalX(end);...
              intervalX(end) intervalX(end)],...
              [intervalY(1) intervalY(end);...
              intervalY(1) intervalY(1);...
              intervalY(end) intervalY(end);...
              intervalY(end) intervalY(1)],...
              'LineStyle','-','color','r','LineWidth',1.5);

        tempL = left(intervalY,intervalX);
        tempL = tempL - mean(mean(tempL));

        %Plot template image
        imagesc(template, tempL);
        set(template,'ydir','normal');
        daspect(template, [1 1 1])
        template.Title.String = "Template image";

        R = normxcorr2(tempL,right);
        peak(i,j) = max(max(R));

        [ymax,xmax] = find(R==peak(i,j));

        % Plot contour map
        corrplot = subplot(2,2,3);
        [tempx, tempy] = meshgrid(1:size(R,1),1:size(R,1));
        surf(corrplot, tempx, tempy, R, 'EdgeColor', 'none');
        view([-1 -1 1])
        xlim([0 1024]);
        ylim([0 1024]);
        zlim([-1 1]);
        colormap(corrplot, jet(64))
        title(corrplot, ['Correlation map. Peak R=', num2str(peak(i,j))]);
        colorbar;
        caxis([-1 1])
        axis([0 1024 0 1024]);

        %normxcorr2 adds padding equal to window size. Look at size(R) vs
        %size(image) for reference.
        yOffset = ymax - winsize;
        xOffset = xmax - winsize;

        deltax(i,j) = xOffset - winsize*j+1; 
        deltay(i,j) = yOffset - winsize*i+1;

        %live contour plot updating as loop runs
        surf(testfig,xgrid,ygrid,deltax);
        view(testfig,2)
        daspect(testfig,[1 1 1])
        axis(testfig,[0 1024 0 1024])
        caxis(testfig,[-150 0])
        colormap(testfig,jet(64))
        testfig.Title.String = "Live updating contour plot";
    end
end

%% clean up of potentially spurious vectors
CORRELATION_CUTOFF = 0.3;
grid_idx = find(peak<CORRELATION_CUTOFF);
deltax_cleaned = deltax;
deltay_cleaned = deltay;
deltax_cleaned(grid_idx) = NaN;
deltay_cleaned(grid_idx) = NaN;
dxpadded = padarray(deltax_cleaned,[1 1],'replicate','both'); % pad array to make the next mean of neighbours step easier. No index out of bounds problem
dypadded = padarray(deltay_cleaned,[1 1],'replicate','both');

for a = 1:size(dxpadded,1)-2 %2 accounts for unit padding on both sides
    m = a + 1; % +1 accounts for initial padding on left/top
    for b = 1:size(dxpadded,2)-2
        n = b + 1;
        if isnan(deltax_cleaned(a,b))
            deltax_cleaned(a,b) = mean(dxpadded(m-1:m+1,n-1:n+1),'all', 'omitnan');
            deltay_cleaned(a,b) = mean(dypadded(m-1:m+1,n-1:n+1),'all', 'omitnan');
        end
    end
end

%% Image plotting
if PLOTTING
    hLeft = figure;
    ax1 = axes;
    imagesc(ax1, left);
    set(ax1, 'ydir', 'normal')
    colormap gray;
    daspect([1 1 1]);
    axis([0 1024 0 1024]);
    title('Left image');
    xlabel('x (pixels)');
    ylabel('y (pixels)');
    
    xline(xgrid(1,:), 'LineStyle','-','Color','w','LineWidth',1.2);
    yline(ygrid(:,1), 'LineStyle','-','Color','w','LineWidth',1.2);
    
    hRight = figure;
    ax2 = axes;
    imagesc(ax2, right);
    set(ax2, 'ydir', 'normal')
    colormap gray;
    daspect([1 1 1]);
    axis([0 1024 0 1024]);
    title('Right image');
    xlabel('x (pixels)');
    ylabel('y (pixels)');
    
    xline(xgrid(1,:), 'LineStyle','-','Color','w','LineWidth',1.2);
    yline(ygrid(:,1), 'LineStyle','-','Color','w','LineWidth',1.2);
end

%% Plotting of x-contour plot
if PLOTTING
    figure
    ax1 = axes;
    im = imagesc(ax1, left);
    colormap(ax1, "gray")
    set(ax1, 'ydir', 'normal')
    hold all;
    axis square;
    xline(xgrid(1,:), 'LineStyle','-','Color','w','LineWidth',1.2);
    yline(ygrid(:,1), 'LineStyle','-','Color','w','LineWidth',1.2);

    title('Contour plot, x-displacement');
    xlabel('x (pixels)');
    ylabel('y (pixels)');

    ax2 = axes;
    s = surf(ax2, xgrid, ygrid, deltay_cleaned);
    view(2)
    axis square;
    hold off;
    s.FaceAlpha = 0.5;
    s.EdgeColor = 'none';
    ax2.GridLineStyle = 'none';
    ax2.Color = 'none';
    colormap(ax2, jet(64));
    caxis([-150 0]);
    c = colorbar('position', [0.85, 0.1, 0.025, 0.8]);
    c.Label.String = 'x-displacement (pixels)';
    linkaxes([ax1,ax2])
end

%% Plotting showing feature on each image
if PLOTTING
    figure
    ax1 = axes;
    im = imagesc(ax1, left);
    hold on
    colormap(ax1, "gray")
    set(ax1, 'ydir', 'normal')
    axis square;
    rectangle('Position',[321, 449, 32, 32],'EdgeColor','r','LineWidth',1.2);

    title('Left image with hidden feature');
    xlabel('x (pixels)');
    ylabel('y (pixels)');


    figure
    ax2 = axes;
    im = imagesc(ax2, right);
    hold on
    colormap(ax2, "gray")
    set(ax2, 'ydir', 'normal')
    axis square;
    rectangle('Position',[191,449, 32, 32],'EdgeColor','r','LineWidth',1.2);

    xline(xgrid(1,:), 'LineStyle','-','Color','w','LineWidth',1.2);
    yline(ygrid(:,1), 'LineStyle','-','Color','w','LineWidth',1.2);

    title('Right image with hidden feature');
    xlabel('x (pixels)');
    ylabel('y (pixels)');

end



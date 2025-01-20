% calulateDisplacements
% PIV package
% Author: Francesco Mascadri
% Contact: fmascadri@student.unimelb.edu.au
% April 2022

function displacementGrid = calculateDisplacementsv2(imagea, imageb, wsize, xgrid, ygrid)
    %calculateDisplacements

    arguments
        %imagea The first image at time t
        imagea (:,:) double {mustBeNumeric, mustBeNonempty}
        %imageb The second image at time t+dt
        imageb (:,:) double {mustBeNumeric, mustBeNonempty}
        %wsize x,y interrogation window size [pixels]
        wsize (1,2)  {mustBeInteger, mustBeNonempty, mustBeNonzero}
        %xgrid Array of x positions dictating the centrepoints of all
        %interrogation windows
        xgrid (:,:) {mustBeNumeric, mustBeNonempty}
        %ygrid Array of y positions dictating the centerpoints of all
        %interrogation windows
        ygrid (:,:) {mustBeNumeric, mustBeNonempty}
    end

    %Intialise displacement grid
    displacementGrid = crosscorr.DisplacementGrid();
    displacementGrid.dpx = zeros(size(xgrid,1),size(xgrid,2));
    displacementGrid.dpy = zeros(size(ygrid,1),size(ygrid,2));

    %Cross-correlation
    for i = 1:size(ygrid,1)
        m = ygrid(i,1);
        for j = 1:size(xgrid,2)
            n = xgrid(1,j);

            x1 = n - wsize(1)/2;
            x2 = n + wsize(1)/2 - 1;
            y1 = m - wsize(2)/2;
            y2 = m + wsize(2)/2 - 1;

            kernel = imagea(x1:x2,y1:y2);
            kernel = kernel - mean(mean(kernel));

            % cross-correlation of image pairs
            R = normxcorr2(kernel,imageb);
            [ymax, xmax] = find(R==(max(max(R))));

            ymax = ymax(1);
            xmax = xmax(1);

            displacementGrid.dpx(i,j) = xmax - wsize(1);
            displacementGrid.dpy(i,j) = ymax - wsize(2);

            logger = fx.log4m.getLogger;
            logger.debug('calculateDisplacements', ['i = ' num2str(i) ', j = ' num2str(j)]);
        end
    end

    %Log ouput
    logger = fx.log4m.getLogger;
    logger.trace('PIV:calculateDisplacements','calculateDisplacements called');
end
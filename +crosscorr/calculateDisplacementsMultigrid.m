% calulateDisplacementsMultigrid
% PIV package
% Author: Francesco Mascadri
% Contact: fmascadri@student.unimelb.edu.au
% April 2022

function displacementGrid = calculateDisplacementsMultigrid(imagea, imageb, wsize, xgrid, ygrid, dpx_est, dpy_est)
    %calculateDisplacementsMultigrid Calculates displacements with a
    %multigrid approach

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
        %dpx_est Array of estimates (to the nearest integer) of the
        %x-displacement of particles in each interrogation window
        dpx_est (:,:) {mustBeInteger, mustBeNonempty}
        %dpy_est Array of estimates (to the nearest integer) of the
        %y-displacement of particles in each interrogation window
        dpy_est (:,:) {mustBeInteger, mustBeNonempty}
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

            % modification for multi-grid approach
            x1_b = x1 + dpx_est(i,j);
            x2_b = x2 + dpx_est(i,j);
            y1_b = y1 + dpy_est(i,j);
            y2_b = y2 + dpy_est(i,j);

            % catch error if displacements put array indices out of bounds
            try 
                imageb(y1_b:y2_b,x1_b:x2_b);
            catch
                x1_b = x1;
                x2_b = x2;
                y1_b = y1;
                y2_b = y2;
            end

            %normalisation
            image1 = imagea(y1:y2,x1:x2);
            image2 = imageb(y1_b:y2_b,x1_b:x2_b); % note modified interrogation window
            image1n = image1 - mean(mean(image1));
            sigma_1 = sqrt(mean(mean(image1n.^2)));
            image2n = image2 - mean(mean(image2));
            sigma_2 = sqrt(mean(mean(image2n.^2)));
            image1n = image1n./sigma_1;
            image2n = image2n./sigma_2;

            % cross-correlation of image pairs
            C = xcorr2(image2n, image1n)./(wsize(1)*wsize(2));
            [ymax, xmax] = find(C==(max(max(C))));

            ymax = ymax(1);
            xmax = xmax(1);

            % guards against index out of bounds error
            if xmax == 1 
                y2 = C(ymax, xmax);
                y1 = C((ymax - 1), xmax);
                y3 = C((ymax + 1), xmax);

                x2 = C(ymax, xmax);
                x1 = C(ymax, xmax); % replace value outside array bounds with max x
                x3 = C(ymax, (xmax + 1));
            elseif xmax == size(C,2)
                y2 = C(ymax, xmax);
                y1 = C((ymax - 1), xmax);
                y3 = C((ymax + 1), xmax);

                x2 = C(ymax, xmax);
                x1 = C(ymax, (xmax - 1));
                x3 = C(ymax, xmax); % replace value outside array bounds with max x
            elseif ymax == 1
                y2 = C(ymax, xmax);
                y1 = C(ymax, xmax); % replace value outside array boundsd with max y
                y3 = C((ymax + 1), xmax);

                x2 = C(ymax, xmax);
                x1 = C(ymax, (xmax - 1));
                x3 = C(ymax, (xmax + 1));
            elseif ymax == size(C,2)
                y2 = C(ymax, xmax);
                y1 = C((ymax - 1), xmax); 
                y3 = C(ymax, xmax); % replace value outside array boundsd with max y

                x2 = C(ymax, xmax);
                x1 = C(ymax, (xmax - 1));
                x3 = C(ymax, (xmax + 1));
            else 
                %Default gaussian sub-pixel fitting
                y2 = C(ymax, xmax);
                y1 = C((ymax - 1), xmax);
                y3 = C((ymax + 1), xmax);
    
                x2 = C(ymax, xmax);
                x1 = C(ymax, (xmax - 1));
                x3 = C(ymax, (xmax + 1));
            end

            % guards against spurious cross-correlation value next to the
            % peak correlation point causing the deltax and/or deltay to be
            % complex (due to taking logarithm of negative number)
            if x1 <= 0
                x1 = 1e-6;
            end
            if x2 <= 0
                x2 = 1e-6;
            end
            if x3 <= 0
                x3 = 1e-6;
            end
            if y1 <= 0
                y1 = 1e-6;
            end
            if y2 <= 0
                y2 = 1e-6;
            end
            if y3 <= 0
                y3 = 1e-6;
            end

            deltay = (log(y1) - log(y3))/(log(y1)+log(y3)-2*log(y2))/2;
            deltax = (log(x1) - log(x3))/(log(x1)+log(x3)-2*log(x2))/2;

            % error checking - better to throw error and kill processing
            % asap instead of finding out after it's done processing
            if imag(deltax)
                imag(deltax)
                error('complex number in x')
            end

            if imag(deltay)
                imag(deltay)
                error('complex number in y')
            end

            displacementGrid.dpx(i,j) = xmax + deltax - wsize(1);
            displacementGrid.dpy(i,j) = ymax + deltay - wsize(2);
        end
    end

    %Log ouput
    logger = fx.log4m.getLogger;
    logger.trace('PIV:calculateDisplacementsMultigrid',...
        'calculateDisplacementsMultigrid called');
end
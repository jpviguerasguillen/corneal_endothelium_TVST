function [imEdges, imEdgesAlt] = applyCellWatershedFreq(imgProb, fc, flagBorder, flagCorrt, k_sigm)
% It computes the segmentation of a (corneal endothelium) image by applying 
% watershed to an already probabilistic input image. It follows a simplified 
% version of the approach suggested by Selig et al (2015).
% It needs to receive (as input) the characteritic frequency.
%
% STEPS:
% 1. Apply a Gaussian filter (input variable required) to smooth the image.
%    This is a very delicate part if a good segmentation is desired.
% 2. The classical watershed is applied to the image. This provides the 
%    final segmentation.
% 3. (OPTIONAL). A final step to improve the segmentation is done. 
%
% This function requires the following functions: 
%   - DIP Toolbox: watershed (https://diplib.org/)
%
% INPUT:
%   imgProb:    Probability image of a cell image
%   fc:         Characteristic frequency
%   flagCorrt:  It applies 'edge correction' at the end.
%   k_sigm:     Parameter related to the sigma of the Gaussian filter.
%   flagBorder: To add a border to the Probability image, which can help
%               to close cells (in the border of the image) during the 
%               watershed
%
% OUTPUT:
%   imEdges:    (Binary image) Image with the detected edges (1 = edge).
%   imEdgesAlt: (Binary image) If flagCorrt was used, the corrected image 
%               is returned in 'imEdges' and the non-corrected image in 
%               'imEdgesAlt' (just for visual analysis).


% This code needs DIP toolbox 
% run('C:\Program Files\DIPimage 2.9\dipstart.m')


%% Parameters
if length(size(imgProb)) ~= 2
  error('Error. \nInput image must be 2D.')
end

% Parameters for watershed. This is the default input arguments
if (nargin < 3 || isempty(flagBorder)); flagBorder = true;  end
if (nargin < 4 || isempty(flagCorrt));  flagCorrt = true;   end  % Apply cell correction.
if (nargin < 5 || isempty(k_sigm));     k_sigm    = 0.20;   end  % Value for the Gaussian smoothing filter

[N, M] = size(imgProb);


%% Add a border in the image if indicated
%  Add 1 pixel around the image and set a value half of a normal edge.
%  Benefits:
%  1. We can close "open cells" that are close to the border, which help
%     with the watershed segmentation.
%  2. If there is a strong, true edge close to the border, setting a border
%     with half intensity would not distort/move the good edge when the
%     smoothing filter is applied.
if flagBorder
  imgProbX = 128*ones(N+2, M+2);
  imgProbX(2:end-1, 2:end-1) = imgProb;
  imgProb = imgProbX;
end


%% Apply a Gaussian filter
sigPDF  = k_sigm/fc;                            % This makes sigmaPDF independent of the cell size
ax      = -3*sigPDF : 1 : 3*sigPDF;          
Gauss1  = normpdf(ax,0,sigPDF);                 % This is a 1D Gaussian
imGauss = conv2(Gauss1,Gauss1,imgProb,'same');  % Apply Gaussian (2D) in the image


% Apply the classic watershed to the smoothed likelihood map
connectWS = 1; % 1 means 4-connected neighbours
max_depth = 2; 
max_size  = 2;
image_out = watershed(imGauss, connectWS, max_depth, max_size); % This watershed function is from DIPimage
imEdges = uint8(image_out);


%% If a border was added, reduce the image to its original size.
if flagBorder
    imEdges = imEdges(2:end-1,2:end-1);
    imEdges(1,:)   = 1;
    imEdges(:,1)   = 1;
    imEdges(end,:) = 1;
    imEdges(:,end) = 1;
end


%% Correcting segmentation borders
% The purpose is to shrink the cells a 20% and repeat the classic
% seeded wathershed from those shrinked cells as seeds.

if flagCorrt
    imEdgesAlt = imEdges; % Place the default image as the 2nd returned argument
    
    % We need to create a temporal image (imEdgesX) where the perimeter is  
    % set to zero. This is used only to compute the distance to the edges.
    imEdgesX = imEdges;
    imEdgesX(1,:)   = 0;
    imEdgesX(:,1)   = 0;
    imEdgesX(end,:) = 0;
    imEdgesX(:,end) = 0;
        
    % First, compute a mask so the inner 80% of the cells are found.
    imDT2Ed = bwdist(imEdgesX);                 % The distance to the edges
    imMarkr = zeros(N,M);                       % The marker image
    imCells = abs(1-imEdges);                   % The image where the cells are 1s.
    [imLabel, numLabel] = bwlabel(imCells,4);   % Label each cell
    for ii = 1 : numLabel
        iCell = find(imLabel == ii);            % All the coordinate of the pixels of one cell
        if numel(iCell) > 1                     % If a cell is only 1 pixel, this is a noisy cell. Do not include it.
            [R,C] = ind2sub([N,M],iCell);       % The coordinates in Rows and Cols
            mR = mean(R);                       % Look for the center of the cell
            mC = mean(C);
            dC = sqrt(sum(([R, C] - repmat([mR,mC], [numel(R),1])).^2, 2)); % The Euclidean distance of each pixel to the center 
            dE = imDT2Ed(iCell);                                            % The Euclidean distance of each pixel to the nearest edge
            ratSh = dC./(dC+dE);                                            % The ratio to determine if a pixel is too close to the edge
            imMarkr(iCell(ratSh < 0.8)) = 1;                                % The inner pixels are marked with 1s in the marker image
        end
    end
    
    % Before applying waterseed, see that we expect that the cells have a 4-conn property. 
    % Thus, the markers need to be also connected in a 4-conn way.
    % Remove any pixel-marker that has no connections in a 4-conn
    filtA = [0 1 0; 1 0 1; 0 1 0];
    imMarkA = logical(conv2(imMarkr,filtA, 'same'));
    imMarkB = imMarkr .* imMarkA;

    % Second, smooth a little the original image (imProb) before employing the watershed
    sigS = 2;  % Sigma to smooth the original image (value from Selig's paper)
    ax     = -3*sigS : 1 : 3*sigS;          
    Gauss2 = normpdf(ax,0,sigS);
    imgIs  = conv2(Gauss2,Gauss2,imgProb,'same');  % Apply Gaussian
    
    % If we used flagBorder, 'imgProb' had a lager size and thus imgIs is 
    % also larger. Just remove the border.
    if flagBorder
      imgIs = imgIs(2:end-1,2:end-1);
    end
    
    % Apply the seeded watershed
    image_out = waterseed(logical(imMarkB), imgIs, connectWS, max_depth, max_size); % Apply seeded watershed
    imEdges = uint8(image_out);
    
    % Finally, note that the border of the image is zero. Make it 1 to keep
    % consistency with the method.
    imEdges(1,:)   = 1;
    imEdges(:,1)   = 1;
    imEdges(end,:) = 1;
    imEdges(:,end) = 1;
    
else
    imEdgesAlt = 0;
end

end


function [fc1, sg1] = obtain_freqAnalysis(imgProb, flagFit, flagPrint)
% It obtains the characteristic frequency (fc) from a probabilistic image. 
%
% INPUT:
%   imgProb:    Probability image of a cell image
%   flagFit:    Use a fitting equation to find the characteristic
%               frequency. By default, it takes the frequency with the 
%               highest modeFM value. 
% OUTPUT:
%   fc:         The characteristic frequency.
%   sg:         The sigma in the Gaussian used for fitting in the frequency
%               analysis.
%
% This function requires the following functions: 
%   - Estimation of cell density: freqAnalysis_radialMean

%% Parameters
if (nargin < 3 || isempty(flagPrint)); flagPrint = true;  end

if length(size(imgProb)) ~= 2
  error('Error. \nInput image must be 2D.')
end

% Parameters for watershed. This is the default input arguments
% By default, use fitting in the Fourier Analysis 
if (nargin < 2 || isempty(flagFit));    flagFit  = true;   end  
  
% Some input variables for the Frequency Analysis (see such function!)
recDilation = true;   %-
plotGraphs  = false;  %-
HPfilter    = true;   % Depending on the image, use a High pass filter


%% Apply the frequency analysis 
[fc1, sg1] = freqAnalysis_radialMean(imgProb, recDilation, plotGraphs, HPfilter, [], flagFit);
if flagPrint
  disp(['   The analysis gives freq = ' num2str(fc1, 5)]); 
  disp(['   The analysis gives sigm = ' num2str(sg1, 5)]); 
end

% If the f is too small (or large), most probably the fitting failed.
% Re-do the analysis without fitting.
if (fc1 < 0.02 || fc1 > 0.08) 
  disp('... redoing the frequency analysis without fitting... ');
  fc1 = freqAnalysis_radialMean(imgProb, recDilation, plotGraphs, HPfilter);
  sg1 = 0;
  disp(['   The NEW analysis gives freq = ' num2str(fc1)]);
end



end


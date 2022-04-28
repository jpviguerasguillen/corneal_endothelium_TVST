function [fc, sg, modeRM, freq] = freqAnalysis_radialMean(img, recDilation, plotGraphs, HPfilter, zeroPad, flagFit)
% freqAnalysis_radialMean(img) computes and plots the radial mean of the frequency 
% spectrum of image 'img', and provides the characteristic frequency (f).
% The radial mean is the integral over theta (the polar coordinate, from 0 to 
% 2*pi) of the frequency spectrum, so a 2D signal is converted into a 1D signal
%
% For more details, see paper: Selig 2015, "Fully automatic evaluation
% of the corneal endothelium from in vivo confocal microscopy".
%            
% The steps:
%   1. A 2D Fourier transform of the image is done.
%   2. A vector of frequencies is defined (freq = 0, 1/N,...,1/2), being N
%      the width of the image.
%   3. For each point in the Fourier image, we round it to its neirest freq
%      in the previous vector.
%   4. We compute the mean for each value in freq. This is signal modeRM.
%   5. The characteristic frequency (f) is the highest peak in this signal.
%
% Some options are included:
%   - Apply reconstruction by dilation between steps 3 and 4.
%   - Apply zero-padding to the image before doing the 2D Fourier transform
%   - Apply a high-pass filter to modeRM in order to remove the lowest
%     frequencies, which can be higher than the desired characteristic freq.
%   - Plot the graphs in order to visualize the signals (Fourier and modeRM)
%
% IMPORTANT: It is adviced NOT to zero-pad a conventional image. 
%     Reconstruction by dilation will not work in a zero-padded Fourier
%     image. However, if the input image is the PDF of the stochastic watershed, 
%     zero-padding might be an option if also the high-pass filter is
%     included, but this needs to be supervised. 
%
% NOTE2: High pass filter is adviced depending on the image. Selig images
% have low cell density, which means very low characteristic frequency.
% Thus, HPfilter might remove the actual frequency point. However, for
% other images with higher cell density, HPfilter is necessary.
%   
%
% INPUT
%   img:            Input image
%   recDilation:    (Boolean) Activates the reconstruction by dilation
%   plotGraphs:     (Boolean) Plots the radial mean of the freq. analysis   
%   HPfilter:       (Boolean) Apply a high pass filter so the low
%                   frequencies (less than f_min = 580) are removed. This 
%                   is helpful for blurred images.
%   zeroPad:        (Boolean) Zeropad the input image
%   flagFit:        (Boolean) A flag to launch the fitting function and
%                   find a better frequency peak. 
%   
% OUTPUT
%   fc:             Frequency peak in the radial mean (modeRM). If fitting
%                   was used, fc is the mu of the Gaussian.
%   sg:             The sigma of the Gaussian in the fitting. If fitting
%                   was not applied, it returns 0.
%   modeRM:         The whole modeRM
%   freq:           The array of frequencies related to the modeRM
% 
%
% How is it computed?
% - First, we compute the 2D Fourier Transform of the image.
%   That gives a Fourier image in Cartesian coordinates.
% - In order to approximate to the image in Polar coordinates, compute for
%   each pixel the distance to the center of the Fourier image (this is rho).
% - Round rho and sum all pixels with the same rho.
%
% (C) J.P. Vigueras-Guillen, 2015


%% Process the input arguments

% Error if number of inputs is less than 1 or more than 2
narginchk(1,6); 

% This is the default input arguments
if (nargin < 2 || isempty(recDilation));    recDilation = true;  end
if (nargin < 3 || isempty(plotGraphs));     plotGraphs  = false; end 
if (nargin < 4 || isempty(HPfilter));       HPfilter    = false; end 
if (nargin < 5 || isempty(zeroPad));        zeroPad     = false; end
if (nargin < 6 || isempty(flagFit));        flagFit     = false; end


%% Fourier Transform
% Process image size information
[N, M] = size(img);

% Zero-padding if indicated in the input arguments
if zeroPad
    if isa(img,'uint8')
        Q = 2^nextpow2(max(N,M));    
        img2 = uint8(zeros(Q,Q));    
        img2(1:N,1:M) = img;
        img = img2;
        M = Q;
        N = Q;
    else
        Q = 2^nextpow2(max(N,M));    
        img2 = zeros(Q,Q);    
        img2(1:N,1:M) = im2double(img);
        img = img2;
        M = Q;
        N = Q;
    end    
end

% Compute the Fourier Transform 
imgF  = fftshift(fft2(img));
imgF = (abs(imgF)/(N*M)) ;              % Normalization of the magnitude   


%% Make a Cartesian grid: Adjust the grid to locate correctly f=0
%  If the original image has an EVEN number of elements in one axis, the
%  Fourier transform has also an even number of points in that axis. The
%  f=0 lies then in the previous element of the middle point.
%  If the original image has an ODD number of elements in one axis, the
%  Fourier transform has also an odd number of points, and f=0 lies in the
%  middle point.

if ~mod(N,2)                            % Even rows
    gridN = -N/2 : N/2-1;
else                                    % Odd rows
    gridN = -N/2+0.5 : N/2-0.5;
end
if ~mod(M,2)                            % Even columns
    gridM = -M/2 : M/2-1;
else                                    % Odd columns
    gridM = -M/2+0.5 : M/2-0.5;
end


% Correct the ellipse-effect:
% If the image was NOT square, then the expected ring in the Fourier image
% (where the characteristic frequency f lies) will be elliptical. We just 
% need to adapt the grid.
% However, if zero-padding was done, for sure the image is square.
if N ~= M
    if N > M
        gridM = gridM * N/M;
    else
        gridN = gridN * M/N;
    end
end

% Make the final 2D grid from the 1D grids already computed
[X, Y] = meshgrid(gridM, gridN);        % Remember: meshgrid uses X & Y, opposite of rows & cols
[~, rho] = cart2pol(X, Y);              % Convert to polar coordinates
rho = round(rho);                       % Round rho (a way of sampling the frequency domain)


%% Apply reconstruction by dilation
if recDilation
    H = zeros(N, M);                        % Define the dilated image
    cnt = find(rho == 0);                   % Find the index for f = 0
    H(cnt) = imgF(cnt);                     % Set the initial value (f = 0)
    IMrec = imreconstruct(H, imgF, 4);
    F = imgF - IMrec;
else
    F = imgF;
end


%% Compute radially average frequency spectrum
halfRho = floor(max(N,M)/2) + 1;        % It is enough to set Rho as half of the (largest) side of the image
ind = cell(halfRho, 1);                                                     
for r = 0 : halfRho - 1                 % Take the indices of those with the same rho
    ind{r+1} = find(rho == r);
end
modeRM = zeros(halfRho,1);
for r = 0 : halfRho - 1                 % Compute the mean for each rho value                                       
    modeRM(r+1) = mean(F(ind{r+1}));
end


%% (For study purposes) Plot the spectrum (zoomed) pre and post reconstruction

% Gn1 = round(N*3/10);
% Gn2 = round(N*7/10);
% Gm1 = round(M*3/10);
% Gm2 = round(M*7/10);
% figure; imshow(imgF(Gn1:Gn2, Gm1:Gm2))
% figure; imshow(F(Gn1:Gn2, Gm1:Gm2))


%% Result: Find the peak  by choosing the highest peak.
% This is a simple method to obtain the peak, but so far it's enough.
% For images with very low density or blurred images, the real f can be 
% overshadowed by the low frequencies. 
% Simple approach: Set a minimum frequency (in this case the 4% of the 
% frequencies), so frequencies below that value are set to zero.
% Depending on the type of images, you might tune that limit.
freqLim = 0.04; 

freq = ((1: halfRho)/halfRho/2)';
modeRM(1:2) = 0;         % Remove the two lowest frequencies by default 
modeRMoriginal = modeRM; % Save the original if you want to plot it
freqOriginal = freq;

if HPfilter
    fmin = round(freqLim*length(modeRM));  % This is the frequency at low limit.
    modeRM = modeRM(fmin+1:end);           % Remove those first points  
    freq = freq(fmin+1:end);
end

% If fitting is not used, just take the highest peak in the data
[~,ind] = max(modeRM);
fc = freq(ind);          % This is the most common (mode) freq
hp_val = modeRM(ind);    % This is the value of the highest peak


%% Result: Find the peak by fitting a model to the signal.
try
  if flagFit 
      % Try to fit function that contains an exponential and a Gaussian 
      % (the most close distribution that fits the data). Use the already
      % estimated f1 as a starter for the fitting.
      %
      % NOTE 1: It seems that setting the same sigma parameter in the two  
      % places of the Gaussian makes the convergence fail, but setting
      % different parameters (scaling and spread is not related) make
      % a better convergence. In the end, we want a good estimation of the
      % Gaussian mean, so only think about that!
      %
      % Note 2: Using the default trust-region-reflective algorithm.
      % Levenberg-Marquardt algorithm does not handle bound constraints.

      funY = @(x,xdata) x(1)*exp(-x(2)*xdata) + (1/(x(3)*sqrt(2*pi)))*exp(-(1/2)*((xdata - x(4))/x(5)).^2) + x(6); % The function to fit
      x0 = [0.5; 1; 0.1; fc; 0.005; 0];           % Initial values of the parameters to find (think better initialization!)
      lb = [0;0;0;fc/2;0;0];                      % Lower bounds
      ub = [max(modeRM); Inf; Inf; 2*fc; Inf; 1]; % Upper bounds. (It is important to set good bounds to have good convergence).
      options = optimoptions('lsqcurvefit','Display','off', 'MaxIterations', 1000);
      [xcoeff, ~,~,~,~] = lsqcurvefit(funY,x0,freq,modeRM,lb,ub,options); % Fit the function!
      fc = xcoeff(4);                             % The mean of the Gaussian in the function funY
      %fprintf('The ''trust-region-reflective'' algorithm took %d function evaluations,\n',output.funcCount)
      %fprintf('The ''trust-region-reflective'' algorithm has residual norm %f,\n',resnorm)


      % ####################################################################
      % Update: Redo a refinement in the fitting. 
      % 0) First, checked whether the fitting failed (it failed if the highest 
      %    peak is far from the estimated mu). 
      %    - If succeed: Point to fit new exponential is 2*sigma + mu (after
      %      the Gaussian).
      %    - If failed:  Point to fit is 1.25*highest_peak (sigma is usually 
      %      4 times less than the mu).
      uppLim = 1.2*freq(ind); % The upper limit
      lowLim = 0.8*freq(ind); % The lower limit
      if fc < lowLim || fc > uppLim
          f_sup = 1.25*freq(ind);
      else
          f_sup = xcoeff(4) + 2* xcoeff(5);
      end

      % 1) Fit the exponential from the previous computed point.
      idx = find(freq >= f_sup);
      freqN = freq(idx);
      modeN = modeRM(idx);    
      funZ = @(x,xdata)x(1)*exp(-x(2)*xdata); % Only the exponential
      x0 = [xcoeff(1); xcoeff(2)]; % Use previous values as initial
      lb = [0;0];
      ub = [max(modeRM); Inf];
      [xcoeffB, ~,~,~,~] = lsqcurvefit(funZ,x0,freqN,modeN,lb,ub,options); 

      % 2) Substract the exponential to the signal, but set the values 
      %    previous to the point (2*sigma + mu) to the value at that point.
      valEx = funZ(xcoeffB,freq);
      valEx(1:idx(1)) = valEx(idx(1));
      modeO = modeRM - valEx;
      modeO(modeO < 0) = 0;

      % 3) Fit only the Gaussian to the remaining signal.
      funX = @(x,xdata) (1/(x(1)*sqrt(2*pi)))*exp(-(1/2)*((xdata - x(2))/x(3)).^2); 
      x0 = [0.1; freq(ind); 0.005]; % Initialization is important 
      lb = [0; 0; 0]; 
      ub = [Inf; Inf; Inf];
      [xcoeffC, ~,~,~,~] = lsqcurvefit(funX,x0,freq,modeO,lb,ub,options); 

      % Return the two parameters, mu and sigma
      fc = xcoeffC(2);   
      sg = xcoeffC(3);

      % For study purposes, display all the estimated f and sigmas:
  %     disp(['   The highest peak gives f = ' num2str(freq(ind), 5)]); 
  %     disp(['   The 1st fitting gives  f = ' num2str(xcoeff(4), 5)]); 
  %     disp(['   The 2nd fitting gives  f = ' num2str(fc, 5)]); 
  %     disp(['   The 1st fitting gives sigma = ' num2str(xcoeff(5), 5)]); 
  %     disp(['   The 2nd fitting gives sigma = ' num2str(sg, 5)]); 
  else
      % Without fitting, there is no sigma:
      sg = 0;
  end
catch
  disp('   ----- ERROR in fitting. Using highest peak instead.')
  sg = 0;
end


%% Plot the results
if plotGraphs
%     figure; imshow(log(imgF),[])    % The Fourier image!
%     figure; plot(freq(1:round(2/5*length(modeRM))), modeRM(1:round(2/5*length(modeRM)))) % Plot the results (only until f=0.2)
        
    if flagFit
      figure; hold on;
      plot(freq, modeRM,'b-'); 
      %plot(freqOriginal, modeRMoriginal,'b-'); 
      plot(freq,funY(xcoeff,freq),'r-')
      legend('Data', 'Fitting') % , 'Characteristic frequency'
      fontSize = 14;
      xlabel('\it f\rm (\mum^{-1})','FontSize',fontSize)
      ylabel('\it F_R_M','FontSize',fontSize)
      title(['Radial mean of the frequency spectrum. f=' num2str(fc,'%.4f')],'FontSize',fontSize)
      
      % The exponential
      figure; hold on;
      plot(freq, modeRM,'b-'); 
      plot(freq,valEx,'r-')
      legend('Data', 'Fitting')
      xlabel('\it f\rm (\mum^{-1})','FontSize',fontSize)
      ylabel('\it F_R_M','FontSize',fontSize)
      
      % The Gaussian alone
      figure; hold on;
      plot(freq, modeO,'b-'); 
      plot(freq,funX(xcoeffC,freq),'r-')
      legend('Data - Gaussian', 'Fitting')
      xlabel('\it f\rm (\mum^{-1})','FontSize',fontSize)
      ylabel('\it F_R_M','FontSize',fontSize)

    else
      figure; hold on;
      plot(freqOriginal, modeRMoriginal,'b-')
      scatter(fc,modeRM(ind))
      legend('Data') % , 'Characteristic frequency'
      fontSize = 14;
      xlabel('\it f\rm (\mum^{-1})','FontSize',fontSize)
      ylabel('\it F_R_M','FontSize',fontSize)
      title(['Radial mean of the frequency spectrum. f=' num2str(fc,'%.4f')],'FontSize',fontSize)
    end
end


%% Density is related to f by a factor (alpha) which depends on the shape
% and regularity. For the sake of simplicity, alpha can be considered 1. 
% The relation between pixel and size in our images is: 1 pixel = 1.038 um.

% delta_f = (f/0.001038)^2;             % This is cells/mm^2 (for our images where )
% l_width = 1/f;                        % This is the most common cell width
% alpha = 1;                            % In the paper, they set alpha to 1
% A_averg = alpha.*l_width^2;           % ... so they estimate the average area such this

% Check manually that width
% figure(); imshow(img)
% [x,y] = ginput(2);
% witdh = sqrt((x(1)-x(2))^2 + (y(1)-y(2))^2);


end
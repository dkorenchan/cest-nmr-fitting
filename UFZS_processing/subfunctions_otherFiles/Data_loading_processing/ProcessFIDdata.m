% ProcessFIDdata: Takes 2D NMR data and performs apodization/filtering,
% zero-filling, and Fourier transform based upon specified parameters.
%
%   INPUTS:
%       fid         -   2D array containing FIDs to be processed, with the 
%                       time axis along the 2nd dimension
%       ppars       -   Struct containing numerical/string processing 
%                       parameters 
%       pflgs       -   Struct containing values of logical flags for data
%                       processing options and specifications
%       datapars    -   Struct containing values extracted from data
%                       parameter files
%
%   OUTPUTS:
%       spec        -   2D array containing processed spectra, with the
%                       frequency axis along the 2nd dimension
%
function spec=ProcessFIDdata(fid,ppars,pflgs,datapars)
% Calculate the filter to apply to the FID
if strcmp(ppars.filter,'exponential')
    % First, we'll need the dwell time (in s) to generate a time array
    if pflgs.jeol
        np=datapars(1).x_X_CURR_POINTS;
        aqStart=datapars(1).x_X_START;
        aqEnd=datapars(1).x_X_STOP;
        dw=(aqEnd-aqStart)./(np-1);
    else
        dw = 1/datapars(1).sw/2; %dwell time 
    end
    
    t=(0:ceil(size(fid,2)/2-1))*dw;
    expvec = exp((t-max(t))*ppars.ap*pi);
    expvec(floor(size(fid,2)/2+1):size(fid,2)) = exp(-t*ppars.ap*pi);
    apdmat = repmat(expvec,size(fid,1),1); %exponential
elseif strcmp(ppars.filter,'gaussian')
    % Calculate parameters for Gaussian filter, generate
    gmean = size(fid,2)/2 + 0.5;
    sigma = sqrt(5/log(10)*(gmean - 1)^2 / ppars.edge); %calculate sigma values for Gaussian
    apdmat = repmat(exp(-1 * ((1:size(fid,2)) - gmean).^2 ./ 2 ./ (sigma^2)),...
        [size(fid,1) 1]);
end

% Apodize, zerofill, Fourier transform in spectral dimension only
fidap = fid .* apdmat; %apodize
%     fidap(:,[1:ppars.zeropts (size(fid,2)-ppars.zeropts+1):size(fid,2)]) = 0; 
%         %zero out selected points
fidzf = padarray(fidap,[0,floor((size(fid,2) * (ppars.zf - 1))/2)],'both'); 
    %zerofill in spectral dimension
spec = fftshift(fft(fftshift(fidzf,2),[],2),2); 
    %cut out 1st points for group delay compensation
end
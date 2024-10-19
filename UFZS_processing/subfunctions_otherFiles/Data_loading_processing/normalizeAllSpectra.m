% normalizeAllSpectra: Performs normalization along the vertical scale for
% all z-spectra by one of two ways, indicated in the input params: 
% (1)   Use the two non-saturated z-spectra to determine scaling factors to 
%       get them to overlap, then cluster all spectra with either one or 
%       the other non-saturated z-spectrum and apply the appropriate 
%       scaling; or
% (2)   Normalize all z-spectra such that they have the same value at a
%       specified ppm value
%
%   INPUTS:
%       results     -   Struct containing processed results, including 
%                       spectra     
%       nosatidx    -   Vector containing the indices of the spectra
%                       corresponding with zero saturation amplitude. These
%                       correspond with the first dimension of
%                       results.spec.
%       params      -   Struct containing the relevant parameters for 
%                       normalization, based upon params.normtype:
%                           'nosatfix'  -   performs (1) above (if it can)
%                           'ppmval'    -   performs (2) above, using ppm
%                                           value params.ppmnormind
%
%   OUTPUTS:
%       results     -   Struct containing processed results, now updated  
%                       with normalized spectra 
%       normpars    -   Struct containing parameters relevant to the
%                       'nosatfit' normalization, if performed
%
function [results,normpars]=normalizeAllSpectra(results,nosatidx,params)
normpars=struct;
if strcmp(params.normtype,'nosatfit') && length(nosatidx)>=2 
    %normalize using scaling between non-saturated spectra, plus 
    %k-means clustering
    disp('Fitting non-saturated spectra for normalization...')
    nosat=abs(results.spec(nosatidx([1 end]),:));        
    nosatopts=optimoptions("lsqnonlin","DiffMinChange",1/size(results.spec,2)/1E3); 
        %THIS IS IMPORTANT! If min change in x is too low, it won't be able to
        %detect the gradient for x(2) and x(3) due to rounding
    % First, ID how to scale the 1st non-saturated spectrum to nearly
    % match the 2nd
    normpars.adjFactors=lsqnonlin(@(x) fitNoSatProfiles(x,results.specppm,...
        params.ppmwdw,nosat(1,:),nosat(2,:)),[1 1 0],[.9 .8 -.003],...
        [1.1 1.2 .003],nosatopts);       
    
    % Then cluster the z-spectra using spectral values outside 
    % params.ppmwdw, and adjust all those which cluster with the 2nd 
    % non-saturated spectrum (DON'T INCLUDE OFFSET! Peaks should stay 
    % aligned)
%         evalppm=15;
    normpars.segVec=kmeans(abs(results.spec(:,results.specppm<-params.ppmwdw|...
        results.specppm>params.ppmwdw)),2);
    % Enforce the non-sat spectra as being in different groups
    if normpars.segVec(nosatidx(1)) == normpars.segVec(nosatidx(2))
        if normpars.segVec(nosatidx(1)+1) ~= normpars.segVec(nosatidx(2)-1)
            normpars.segVec(nosatidx(1)) = normpars.segVec(nosatidx(1)+1);
            normpars.segVec(nosatidx(2)) = normpars.segVec(nosatidx(2)-1);
        else
            warning('Non-sat spectra were clustered together!');
        end
    end
    for i=1:length(normpars.segVec)
        if normpars.segVec(i)==normpars.segVec(nosatidx(1))
            results.spec(i,:)=adjZProfile(results.spec(i,:),...
                [normpars.adjFactors(1:2) 0]);
        end
    end
elseif strcmp(params.normtype,'nosatfit') && length(nosatidx)<2
    %detect if not enough non-saturated spectra for this normalization
    warning(['Fewer than two non-saturated spectra detected! Cannot use '...
        'scaling to normalize -- try the other normalization method instead'])

elseif strcmp(params.normtype,'ppmval') %normalize to one ppm signal value
    ppmnormind = find(abs(results.specppm-params.ppmnorm)==...
        min(abs(results.specppm-params.ppmnorm))); 
        %signal index for normalization
    if length(ppmnormind) > 1 %catch if >1 value found
        ppmnormind = ppmnormind(1);
    end
    normvec = abs(results.spec(1,ppmnormind)) ./ ...
        squeeze(abs(results.spec(:,ppmnormind)));
        %normalization factor for each B1 value
    normmat = repmat(normvec,[1 size(results.spec,2)]);
    results.spec = results.spec .* normmat;
end
end


% fitNoSatProfiles: Optimization function used to take two raw z-profiles 
% and determine scaling + offset parameters to minimize the residual 
% between them
function res=fitNoSatProfiles(x,ppm,ppmlim,spec1,spec2)
spec1vh=adjZProfile(spec1,x);
fitidx=(ppm>-ppmlim)&(ppm<ppmlim);
res=abs(spec2(fitidx)-spec1vh(fitidx));
end


% adjZProfile: Scales and offsets a spectral profile using the given
% parameters
function newspec=adjZProfile(spec,x)
specv=spec.*x(1);
bigpts=length(spec)*1E3; %larger number of points to bring things up to
specvbig=interp1(1:length(spec),specv,linspace(1,length(spec),bigpts)); 
    %interpolate many more points
delta_pts=round((1-x(2))*bigpts/2);
offset_pts=round(x(3)*bigpts);
if delta_pts<0
    specvbig([1:(abs(delta_pts)+offset_pts) (end+delta_pts-offset_pts+1):end])=[]; 
    %remove points at ends
else
    specvbig=[zeros(1,(delta_pts+offset_pts)) specvbig zeros(1,(delta_pts-offset_pts))];
end
newspec=interp1(1:length(specvbig),specvbig,linspace(1,length(specvbig),length(spec)));
    %back to original size
end
function [spectrum, yag_spec, residual] = compare_multi_spec(beam,yag_axis,yag_specs,gaussFilter)
%
% This function compares data with LiTrack simulations.
% 
% INPUTS: 
%
%   beam:       z = beam(:,1) [m]
%               e = beam(:,2) [GeV]
%   Nbin:       # of bins for histogram
%   yag_axis:   data Delta axis [%], must be same length as Nbin
%   yag_spec:   data Delta spectrum, must be same length as Nbin
%   smear:      (optional) a struct of smear parameters
%
% OUTPUTS:
%
%   z_ax:       data and LiTrack Z axis [microns]
%   lit_z_prof: LiTrack current profile binned on z_ax
%   e_ax:       data and LiTrack delta axis [%]
%   lit_e_spec: LiTrack Delta spectrum binned on e_ax
%   yag_spec:   Normalized YAG spectrum
%   residual:   Difference between LiTrack and data

% Convert LiTrack beam from [m,GeV] to [um,%]
lit_d = 100*(beam(:,2)-mean(beam(:,2)))/mean(beam(:,2));

% Histogram on shared axis, zero end bins, scale
[ne,eb] = hist(lit_d,yag_axis);
ne(1) = 0;
ne(end) = 0;

if nargin > 3
    
    ne = conv(ne, gaussFilter, 'same');
end

% lit_emax = max(ne);
% dat_emax = max(yag_spec);
% yag_spec = (lit_emax/dat_emax)*yag_spec;

% normalize to 1
ne = ne/sum(ne);

% Copy results
spectrum = ne';
spectra = repmat(spectrum,1,size(yag_specs,2));

[residual,yag_spec] = min(sum((yag_specs - spectra).^2));

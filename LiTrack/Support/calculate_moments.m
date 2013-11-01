function [SIGMA_Z,RMS_Z,I_PEAK,FWHM_D,RMS_D] = calculate_moments(profile,spectrum)

SIGMA_Z = mygaussfit(profile(:,1),profile(:,2));
RMS_Z   = calc_rms(profile(:,1),profile(:,2));
I_PEAK  = max(profile(:,2));

FWHM_D  = FWHM(spectrum(:,1),spectrum(:,2));
RMS_D   = calc_rms(spectrum(:,1),spectrum(:,2));


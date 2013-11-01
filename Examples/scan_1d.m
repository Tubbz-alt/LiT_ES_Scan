% load wake
global wake;
wake = load('slac.dat');

% load params
global PARAM;
param_5mm_R56;
init_beam = 1;

% scan parameter (see file Parameters/param_5mm_R56.m 
% for a full list of available parameters)
param_name = {'NRTL AMPL'};
low_limit  = 0.0400;
high_limit = 0.0410;
n_steps    = 5;

% create scan values
param_vals = linspace(low_limit,high_limit,n_steps);
nPar = 1;
par_str = cell(1,n_steps);
fwhm_str = cell(1,n_steps);
ip_str = cell(1,n_steps);

% color scaling
cmap = colormap;
cind = round(linspace(1,length(cmap),n_steps));

for i = 1:n_steps
    
    % Update parameters by calling SetPars and then executing beamline file
    SetPars(param_vals(i), param_name, nPar);
    FACET_NOTCH_bl;
    
    % Run LiTrack
    LiT_struct = LiT(beamline,init_beam,init_param,PARAM.SIMU.BIN,0);
    
    % Calculate beam moments (outside of LiTrack)
    [sigma_z,rms_z,I_peak,FWHM_d,RMS_d] = calculate_moments(profile,spectrum);
    
    % Generate phase spaces for plots
    [profile,spectrum,phase_space,beam_zd] = gen_ps(LiT_struct,PARAM.SIMU.BIN);
    
    % Plot phase space
    figure(1);
    subplot(2,2,1);
    plot(spectrum(:,2),spectrum(:,1),'color',cmap(cind(i),:),'linewidth',2);
    xlabel('N/bin','fontsize',14);
    ylabel('\delta [%]','fontsize',14);
    hold on;
    subplot(2,2,2);
    plot(beam_zd(:,1),beam_zd(:,2),'color',cmap(cind(i),:),'marker','*','linestyle','none');
    xlabel('Z [\mum]','fontsize',14);
    ylabel('\delta [%]','fontsize',14);
    hold on;
    subplot(2,2,4);
    plot(profile(:,1),profile(:,2),'color',cmap(cind(i),:),'linewidth',2);
    xlabel('Z [\mum]','fontsize',14);
    ylabel('I [kA]','fontsize',14);
    hold on;
    
    % legend strings
    par_str{i} = [char(param_name) ' = ' num2str(1000*param_vals(i),'%0.2f') ' MV'];
    fwhm_str{i} = ['FWHM = ' num2str(FWHM_d,'%0.2f') ' %'];
    ip_str{i} = ['I_{peak} = ' num2str(I_peak,'%0.2f') ' kA'];
end

% add legends
subplot(2,2,1);
legend(fwhm_str);
hold off;
subplot(2,2,2);
legend(par_str);
hold off;
subplot(2,2,4);
legend(ip_str);
hold off;
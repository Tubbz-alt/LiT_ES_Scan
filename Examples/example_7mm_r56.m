SI_consts;

% load wake
global wake;
wake = load('slac.dat');

% load params
global PARAM;
param_7mm_R56;

% Adjust NRTL compressor ampl (GV)
PARAM.NRTL.AMPL = 40.40E-3; % use 40.40E-3 for undercompressed, try 41.20E-3 for overcompressed

% Adjust 2-10 phase (degrees)
PARAM.LONE.PHAS = -21.9; % use -21.9 for undercompressed, try -19.2 for overcompressed

% Adjust notch and jaws (values are in delta)
PARAM.LI20.NLO   = -0.000; % low energy edge of notch
PARAM.LI20.NHI   = 0.000;  % high energy edge of notch
PARAM.LI20.ELO   = -0.025; % low energy (left) jaw
PARAM.LI20.EHI   = 0.025;  % high energy (right) jaw

% load beamline
FACET_NOTCH_bl;

% set beam init params
init_beam = 1;

% run LiT
out = LiT(beamline,init_beam,init_param,PARAM.SIMU.BIN,0);

% generate phase space and axis
[z_proj,e_proj,phase_space,beam_zd] = gen_ps(out,PARAM.SIMU.BIN);

% Plot phase space and profiles
figure(1);
subplot(2,2,1);
plot(e_proj(:,2),e_proj(:,1),'k'); axis tight; ylabel('\delta [%]','fontsize',16);
subplot(2,2,2);
pcolor(z_proj(:,1),e_proj(:,1),phase_space'); shading flat; ylabel('\delta [%]','fontsize',16); xlabel('Z [\mum]','fontsize',16);
subplot(2,2,4);
plot(z_proj(:,1),z_proj(:,2),'k'); axis tight; xlabel('Z [\mum]','fontsize',16);

% Plot phase space
figure(2);
pcolor(z_proj(:,1),e_proj(:,1),phase_space'); shading flat; ylabel('\delta [%]','fontsize',16); xlabel('Z [\mum]','fontsize',16);

% plot fit
figure(3);
[fit_object, error] = peakfit(z_proj,0,0,2);

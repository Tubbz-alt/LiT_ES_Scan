function LT_OUTPUT = LiT(beamline,init_beam,init_param,Nbin,show_all,save_img,save_dir,store_all,store_wake)
%function LT_OUTPUT = LiT(beamline,initialize_beam,init_params,show_all,save_img,save_dir)
%
%
%  INPUTS:
%	beamline:   A matrix of beamline elements and parameterms.		    
%   init_beam:  Boolean to determine whether or not to initialize beam
%   init_param: Set of beam initialization parameters
%
%  OUTPUTS:     
%	LT_OUPUT:	Struct containg beam distributions and moments
%
%==============================================================================================================

if nargin < 6
    save_img = 0;
    save_dir = '';
    store_all = 0;
    store_wake = 0;
end

global beam;
if init_beam
    beam = zeros(init_param.Nesim,2);
    beam(:,1) = init_Z(init_param.Nesim,init_param.sigz0,init_param.asym,init_param.cut,init_param.z0bar);	% generate asymmetric gaussian with cut at N sigma
    beam(:,2) = init_param.E0*(1+init_param.sigd0*randn(init_param.Nesim,1));               % generate gaussian energy distribution
    Nb = length(beam);
    Ne = init_param.Npart;
    Qp = Ne/Nb;
end

% n_el = number of beamline elements
[n_el,~] = size(beamline);
if nargin>=8 && store_all
    evolution = zeros(init_param.Nesim,2,n_el+1);
    evolution(:,:,1) = beam;
    NBs = zeros(1,n_el+1);
    NBs(1) = Nb;
    E0s = zeros(1,n_el+1);
    E0s(1) = mean(beam(:,2));
else
    evolution = [];
    NBs = [];
    E0s = [];
end

% store wake calcs
if nargin==9 && store_wake
    codes = beamline(:,1);
    n_wakes = sum(codes == 11);
    wakes = zeros(Nbin,2,n_wakes);
    wake_ind = 1;
    wake_inds = [];
    rfs = zeros(Nbin,2,n_wakes);
else
    wakes = [];
    rfs = [];
    wake_inds = [];
end

if show_all
    ps.BEAM = beam;
    ps.QP = Qp;
    plot_ps(ps,Nbin,1,1,0,0);
    pause;
    if save_img
        saveas(gcf,[save_dir 'step_00.eps'],'epsc');
    end
end

% Get orders from beamline
for j  = 1:n_el					

  code  = beamline(j,1);

  switch code
      
      case 6
          [beam, Nb] = bunch_compressor(beam,Nb,beamline(j,:));
          if store_all; evolution(1:Nb,:,j+1) = beam; NBs(j+1) = Nb; E0s(j+1) = mean(beam(:,2)); end;
          
      case 11
          [beam, Nb, wake_out, rf_out] = rf_acceleration(beam,Nb,Qp,Nbin,beamline(j,:));
          if store_all; evolution(1:Nb,:,j+1) = beam; NBs(j+1) = Nb; E0s(j+1) = mean(beam(:,2)); end;
          if store_wake 
              wakes(:,:,wake_ind) = wake_out; 
              rfs(:,:,wake_ind) = rf_out; 
              wake_ind = wake_ind+1; 
              wake_inds = [wake_inds; j+1];
          end
          
      case 22
          [beam, Nb] = ISR_energySpread(beam,Nb,beamline(j,:));
          if store_all; evolution(1:Nb,:,j+1) = beam; NBs(j+1) = Nb; E0s(j+1) = mean(beam(:,2)); end;
          
      case 26
          [beam, Nb] = energy_aperture(beam,Nb,beamline(j,:));
          if store_all; evolution(1:Nb,:,j+1) = beam; NBs(j+1) = Nb; E0s(j+1) = mean(beam(:,2)); end;
          
      case 28
          [beam, Nb] = notch_collimator(beam,Nb,beamline(j,:));
          if store_all; evolution(1:Nb,:,j+1) = beam; NBs(j+1) = Nb; E0s(j+1) = mean(beam(:,2)); end;
          
      case 37
          [beam, Nb] = cut_beam_tails(beam,Nb,beamline(j,:));
          if store_all; evolution(1:Nb,:,j+1) = beam; NBs(j+1) = Nb; E0s(j+1) = mean(beam(:,2)); end;
          
  end
      
  if show_all 
      ps.BEAM = beam;
      ps.QP = Qp;
      plot_ps(ps,Nbin,1,1,0,0);
      pause;
      if save_img
          saveas(gcf,[save_dir 'step_' num2str(j,'%02d') '.eps'],'epsc');
      end
  end
  
end % end loop over all beamline section

LT_OUTPUT.BEAM = beam;
LT_OUTPUT.QP = Qp;
LT_OUTPUT.EVO = evolution;
LT_OUTPUT.Nb = NBs;
LT_OUTPUT.E0 = E0s;
LT_OUTPUT.Wakes = wakes;
LT_OUTPUT.RFs = rfs;
LT_OUTPUT.Inds = wake_inds;
function LT_OUTPUT = LiT(beamline,init_beam,init_param,Nbin,show_all)
%function LT_OUTPUT = LiTrackOpt(beamline,initialize_beam,init_params)
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

global beam;
if init_beam
    beam = zeros(init_param.Nesim,2);
    beam(:,1) = init_Z(init_param.Nesim,init_param.sigz0,init_param.asym,init_param.cut);	% generate asymmetric gaussian with cut at N sigma
    beam(:,2) = init_param.E0*(1+init_param.sigd0*randn(init_param.Nesim,1));               % generate gaussian energy distribution
    Nb = length(beam);
    Ne = init_param.Npart;
    Qp = Ne/Nb;
end

% n_el = number of beamline elements
[n_el,~] = size(beamline);		            

if show_all; plot_ps(beam,Nbin,1,1); end;

% Get orders from beamline
for j  = 1:n_el					

  code  = beamline(j,1);

  switch code
      
      case 6
          [beam, Nb] = bunch_compressor(beam,Nb,beamline(j,:));
          
      case 11
          [beam, Nb] = rf_acceleration(beam,Nb,Qp,Nbin,beamline(j,:));
          
      case 22
          [beam, Nb] = ISR_energySpread(beam,Nb,beamline(j,:));
          
      case 26
          [beam, Nb] = energy_aperture(beam,Nb,beamline(j,:));
          
      case 28
          [beam, Nb] = notch_collimator(beam,Nb,beamline(j,:));
          
      case 37
          [beam, Nb] = cut_beam_tails(beam,Nb,beamline(j,:));
          
  end
      
  if show_all; plot_ps(beam,Nbin,1,1); end
  
end % end loop over all beamline section

if nargout == 1

      LT_OUTPUT.BEAM = beam;
      LT_OUTPUT.QP   = Qp;

%     LT_OUTPUT.Z.DIST = zposj;
%     LT_OUTPUT.Z.HIST = HIST_Z;
%     LT_OUTPUT.Z.AXIS = AXIS_Z;
%     %LT_OUTPUT.Z.AVG  = z_barj;
%     LT_OUTPUT.Z.FWHM = ZFWmmj;
%     if gzfit
%         LT_OUTPUT.Z.SIG  = sigzGj;
%         LT_OUTPUT.I.SIG  = I_pkfj;
%     end
% 
%     LT_OUTPUT.E.DIST = dE_Ej;
%     LT_OUTPUT.E.HIST = HIST_D;
%     LT_OUTPUT.E.AXIS = AXIS_D;
%     LT_OUTPUT.E.AVG  = Ebarj;
%     LT_OUTPUT.E.FWHM = dFWpctj;
%     if gdfit
%         LT_OUTPUT.E.SIG  = sigEGj;
%     end
% 
%     LT_OUTPUT.I.PART = fcutj;
%     LT_OUTPUT.I.PEAK = I_pkj;

end
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

if show_all; ps.BEAM = beam; ps.QP = Qp; plot_ps(ps,Nbin,0,1,1,0); end;

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
      
  if show_all; ps.BEAM = beam; ps.QP = Qp; plot_ps(ps,Nbin,0,1,1,0); end
  
end % end loop over all beamline section

LT_OUTPUT.BEAM = beam;
LT_OUTPUT.QP = Qp;
function LT_OUTPUT = LiT(beamline,initialize_beam,init_params)
%function LT_OUTPUT = LiTrackOpt(fn)
%
%	Function to do 2-D longitudinal phase space tracking of many
%	electrons through several accelerator sections including bunch
%	compression and the longitudinal wakefield of the accelerating
%  	structures.
%
%	The beam and the various sections are described in your M-file
%	called "*_lit.m" where "*" is a string contained in the above
%	input argument: fn.
%
%  	The initial zpos and dE/E particle coordinates can be 'LiTrack'-
%	generated Matlab gaussian distributions, uniform distributions, or
%	they can be input from a user's ASCII, 2-column file (Z/mm, dE/E/%).
%
%	If no output arguments are provided, this function generates plots
%	(see below).  If at least one output argument is provided, no plots
%	are generated and the longitudinal coordinates of the N particles
%	are returned as described below.
%
%  INPUTS:
%	fn:		    A string which describes the leading characters of the beam-
%			    line M-file with name "*_lit.m" where the fn string is represented
%			    here as "*" (SEE FOR EXAMPLE the file LCLS_LIT.M which is
%			    internall documented and would be run by setting fn='lcls').
%
%  OUTPUTS:     
%	LT_OUPUT:	Struct containg beam distributions and moments
%
%==============================================================================================================

global beam
if initialize_beam
    beam = zeros(init_params.Nesim,2);
    beam(:,1) = asym_gaussian(init_params.Nesim,init_params.sigz0,init_params.asym,init_params.cut);	% generate asymmetric gaussian with cuts & tails
    beam(:,2) = init_params.sigd0*randn(init_params.Nesim,1);	% always gaussian dE/E   
end

% nb=number of beamline sections in BL-file (e.g. accelerator, compressor, ...)
[n_el,n_spec] = size(beamline);		            


% Get orders from beamline
for j  = 1:n_el					
  code  = beamline(j,1);

  switch code
      
      case 6
          [beam, Nb] = bunch_compressor(beam,beamline(j,:));
          
      case 11
          [beam, Nb] = rf_acceleration(beam,beamline(j,:));
          
      case 13
          [beam, Nb] = energy_feedback(beam,beamline(j,:));
          
      case 22
          [beam, Nb] = ISR_energySpread(beam,beamline(j,:));
          
      case 26
          [beam, Nb] = energy_aperture(beam,beamline(j,:));
          
      case 28
          [beam, Nb] = notch_collimator(beam,beamline(j,:));
          
      case 37
          [beam, Nb] = cut_beam_tail(beam,beamline(j,:));
          
  end
      
  
  if (abs(cod)==11)	|| (abs(cod)==10)	% ACCELERATION SECTION (11 and 10)
    ecuts = 0;					% no dE/E cuts shown on plots
    zcuts = 0;					% no Z cuts shown on plots
    Eacc   = beamline(j,2);		% nominal acc (w/o wake and for phi=crest(=0)) [GeV]
    phi    = beamline(j,3);		% acc phase (crest=0, low-E head at phi < 0) [deg]
    lam    = beamline(j,4);		% RF wavelength [m]
    if lam<=0
      error('RF wavelength cannot be <= 0')
    end
    wakeon = beamline(j,5);		% wakeON=1,2, wakeOFF=0
    Lacc   = beamline(j,6);		% length of acc section (scales wake) [m]
    phir   = phi*pi/180;		% RF phase in radians
    if wakeon					% if wakes calc switched ON...
      iswake = 1;				% turns wake plot on
      nwake_fn = length(wake_fn(:,1));		        % count number of files provided
      if wakeon > nwake_fn
        error('Need multiple wake function file names when "wakeON/OFF" > 1')
      end
      wake_fn1 = strtrim(wake_fn(wakeon,:));			        % select proper wake function depending on wakeon (=1,2,...)
      disp(['Using wake function: ' wake_fn1])		% echo wake file being used
       %[dE_wake,zc_wake] = long_wake(z,Lacc,Ne1,...
       %                    Nbin,wake_fn1);	        % calculate dE_wake in MeV from Z-coordinates, acc length, N-particles, etc.
      [dE_wake,zc_wake] = fast_wake(z,Lacc,Ne1,Nbin,wake_fn1);
      dE_wakes = interp1(zc_wake,dE_wake,z,'*linear');	% inerpolate between wake calc points in bunch to evaluate dE for each e-
      dE_loss = 1E-3*mean(dE_wakes);	% wake loss [GeV]
    else                    			% if wake calc switched OFF...
      iswake = 0;				    	% no wake plot
      dE_wake = zeros(Nbin,1); 			% dE_wake is all zeros
      dE_loss = 0;						% wake loss = 0 without wakes ON
    end
    if abs(cod) == 10					% special case where Eacc is final energy, rather than acc-Voltage
      Eacc = -dE_loss + (Eacc - mean(E))/cos(phir);	% Eacc was final energy, now is acc-volts again [GeV]
    end
    Erf = E + Eacc*cos(phir + 2*pi*z/lam);	% energy of each particle from RF shape alone (no wake yet)
    if wakeon
      E = Erf + dE_wakes*1E-3;		% energy from RF phase and wake added [GeV]
    else
      E = Erf;						% energy from RF phase and NO wake added [GeV]
    end
    Ebar = mean(E);				    % mean particle energy [GeV]
    Ebarcuts = Ebar;				% mean energy after cuts - same as Ebar here [GeV]
    d     = (E - Ebar)/Ebar;		% relative energy deviation w.r.t. mean energy [ ]
  end
  if abs(cod) == 12                	% zero-out all energy deviations for diagnostics ONLY
    d = zeros(size(d));
    E = mean(E)*ones(size(E));
  end
  if abs(cod) == 13                	% energy feedback with two phase-opposed sections, each of eV0 volts @ crest
    ecuts = 0;					    % no dE/E cuts shown on plots
    zcuts = 0;					    % no Z cuts shown on plots
    Efin   = beamline(j,2);		    % Energy setpoint (goal) [GeV]
    eV0    = beamline(j,3);		    % acc. voltage available at crest for each of two fdbk sections [GeV]
    if eV0==0
      error('Feedback voltage of zero will not correct energy')
    end
    phi1r  = beamline(j,4)*pi/180;  % acc phase of 1st section (crest=0, low-E head at phi < 0) [deg]
    phi2r  = beamline(j,5)*pi/180;	% acc phase of 2nd section (crest=0, low-E head at phi < 0) [deg]
    lam    = beamline(j,6);		    % RF wavelength [m]
    if lam<=0
      error('RF wavelength cannot be <= 0')
    end
    iswake  = 0;				    	% no wake plot
    options = optimset;
    dphi    = fminsearch('fdbk_fun',0,options,phi1r,phi2r,(Efin-mean(E))/eV0);
    fprintf('Energy feedback phase set to %8.3f deg\n',dphi*180/pi)
    En      = mean(E) + eV0*cos(phi1r+dphi) + eV0*cos(phi2r-dphi);
    if abs(En-Efin)/Efin > 1E-4
      fprintf('Energy feedback phase maxed out at %8.3f deg\n',dphi*180/pi)
    end
    E = E  + eV0*cos(phi1r+dphi+2*pi*z/lam) + eV0*cos(phi2r-dphi+2*pi*z/lam);
    Ebar = mean(E);				    % mean particle energy [GeV]
    Ebarcuts = Ebar;				% mean energy after cuts - same as Ebar here [GeV]
    d     = (E - Ebar)/Ebar;		% relative energy deviation w.r.t. mean energy [ ]
  end                               % end code==13, energy feedback card
  
  if abs(cod) == 26                	% USER'S ENERGY CUTS (26) - doesn't change Ebar
    ecuts = 1;					    % show dE/E cuts on plots
    d1 = beamline(j,2);				% minimum dE/E to allow through [ ]
    d2 = beamline(j,3);				% maximum dE/E "   "     "      [ ]
    if d1 >= d2					    % bomb out if max<min (BT-file error)
      error(['Energy cuts (26) must have dE/E_min (col 2) < dE/E_max (col3) in ' fnfm])
    end
    i = find(d>d1 & d<d2);			% bomb out if cuts too tight
    if length(i) < 1
      error('Energy cuts (26) Emin=%7.4f %% and Emax=%7.4f %% threw out all particles',d1*100,d2*100)
    end
    Ni = length(i);				    % count particles left after cuts
    Ne1 = Ne1*Ni/Nesim;				% rescale N-particles to reflect cuts
    d = d(i);					    % reduce dE/E array inpose cuts
    z = z(i);					    % reduce Z array inpose cuts
    E = E(i);					    % reduce energy array inpose cuts
    Ebarcuts = mean(E);				% mean energy after cuts [GeV]
    disp([sprintf('E-cut (26): %6.3e',100*(1-Ni/Nesim)) '% of bunch'])  
    Nesim = Ni;					    % reduce number of simulation particles
  end
  if abs(cod) == 28                	% Notch collimator for M. Hogan
%    ecuts = 1;					    % show dE/E cuts on plots
    d1 = beamline(j,2);				% minimum dE/E for notch-collimator edge [ ]
    d2 = beamline(j,3);				% maximum dE/E for notch-collimator edge [ ]
    if d1 >= d2					    % bomb out if max<min (BT-file error)
        disp('No notch cut');
        %error(['Notch-collimator (28) must have dE/E_min (col 2) < dE/E_max (col3) in ' fnfm])
    else
        i = find(d<d1 | d>d2);			% bomb out if notch too wide
        if length(i) < 1
            error('Notch-collimator (28) Emin=%7.4f %% and Emax=%7.4f %% threw out all particles',d1*100,d2*100)
        end
        Ni = length(i);				    % count particles left after cuts
        Ne1 = Ne1*Ni/Nesim;				% rescale N-particles to reflect cuts
        d = d(i);					    % reduce dE/E array inpose cuts
        z = z(i);					    % reduce Z array inpose cuts
        E = E(i);					    % reduce energy array inpose cuts
        Ebarcuts = mean(E);				% mean energy after cuts [GeV]
        disp([sprintf('Notch-collimator (28) cut: %6.3e',100*(1-Ni/Nesim)) '% of bunch'])  
        Nesim = Ni;					    % reduce number of simulation particles
    end
  end

  
  if abs(cod) == 27                         % USER'S constant-dN/N dE/E cuts (27)
%    zcuts = 1;					            % show dE/E-cuts on plots
    dN_N = beamline(j,2);			        % fraction of max-dE/E-amplitude particles to cut [ ]
    no_charge_loss = beamline(j,3);			% if==1, no real charge cut intended, just better binning
    [dsort,idsort] = sort(abs(d-mean(d)));	% sort the absolute value of dE/E values (min to max)
    N1 = round(Nesim*dN_N);			        % throw out last N1 particles
    z(idsort((Nesim-N1):Nesim)) = [];		% now throw them out of zpos
    d(idsort((Nesim-N1):Nesim)) = [];		% now throw them out of dE/E
    E(idsort((Nesim-N1):Nesim)) = [];		% now throw them out of E
    Ni = length(z);				            % count particles left after cuts
    if no_charge_loss==0                    % if real charge cut intended...
      Ne1 = Ne1*Ni/Nesim;				        % ...rescale N-particles to reflect cuts
    end
    Ebarcuts = mean(E);				        % mean energy after cuts [GeV]
    disp([sprintf('Const-dN/N dE/E-cut (27): %6.3f',100*(1-Ni/Nesim)) '% of bunch'])  
    Nesim = Ni;					            % reduce number of simulation particles
  end

  if abs(cod) == 36                     	% USER'S Z-CUTS (36)
    zcuts = 1;					            % show Z-cuts on plots
    z1 = beamline(j,2);				        % minimum Z to allow through [m]
    z2 = beamline(j,3);				        % maximum Z "   "     "      [m]
    if z1 >= z2					            % bomb out if max<min (BT-file error)
      error(['Z-cuts (36) must have Z_min (col 2) < Z_max (col3) in ' fnfm])
    end
    i = find(z>z1 & z<z2);					% bomb out if cuts too tight
    if length(i) < 1
      error('Z-cuts (36) Zmin=%7.4f mm and Zmax=%7.4f mm threw out all particles',z1*1E3,z2*1E3)
    end
    Ni = length(i);							% count particles left after cuts
    Ne1 = Ne1*Ni/Nesim;						% rescale N-particles to reflect cuts
    d = d(i);								% reduce dE/E array inpose cuts
    z = z(i);								% reduce Z array inpose cuts
    E = E(i);								% reduce energy array inpose cuts
    Ebarcuts = mean(E);						% mean energy after cuts [GeV]
    disp([sprintf('Z-cut (36): %6.3f',100*(1-Ni/Nesim)) '% of bunch'])  
    Nesim = Ni;								% reduce number of simulation particles
  end
  if abs(cod) == 37                         % USER'S constant-dN/N z cuts (37)
%    zcuts = 1;					            % show Z-cuts on plots
    dN_N = beamline(j,2);			        % fraction of max-z-amplitude particles to cut [ ]
    no_charge_loss = beamline(j,3);			% if==1, no real charge cut intended, just better binning
    [zsort,izsort] = sort(abs(z-mean(z)));	% sort the absolute value of zpos values (min to max)
    N1 = round(Nesim*dN_N);			        % throw out last N1 particles
    z(izsort((Nesim-N1):Nesim)) = [];		% now throw them out of zpos
    d(izsort((Nesim-N1):Nesim)) = [];		% now throw them out of dE/E
    E(izsort((Nesim-N1):Nesim)) = [];		% now throw them out of E
    Ni = length(z);				            % count particles left after cuts
    if no_charge_loss==0                    % if real charge cut intended...
      Ne1 = Ne1*Ni/Nesim;				        % ...rescale N-particles to reflect cuts
    end
    Ebarcuts = mean(E);				        % mean energy after cuts [GeV]
    disp([sprintf('Const-dN/N Z-cut (37): %6.3f',100*(1-Ni/Nesim)) '% of bunch'])
    Nesim = Ni;					            % reduce number of simulation particles
  end


  if abs(cod) == 6 				% BUNCH COMPRESSION (R56/m, T566/m, E/GeV, U5666/m)
    ecuts = 0;					% no dE/E cuts shown on plots
    zcuts = 0;					% no Z cuts shown on plots
    iswake = 0;					% turn off induced voltage plot
    R56  = beamline(j,2);       % R56 value [m]
    T566 = beamline(j,3);       % T566 value [m] (=-3*R56/2 for non-quad system)
    U5666= beamline(j,5);       % U5666 value [m] (=2*R56 for non-quad system)
    E56  = beamline(j,4);       % Nominal energy of compressor [GeV]
    if E56 < 0.020				% need positive, reasonable nominal R56-energy [GeV]
      fprintf(['WARN: Compressor section (6) of R56=%7.4f m has nominal-energy too small (not ultra-relativistic) in ' fnfm '\n'],R56)
    end
    dd   = (E-E56)/E56;			% relative energy error w.r.t. nominal compressor energy
    z    = R56*dd + T566*dd.^2 + ...
           U5666*dd.^3 + z;		% compress or anti-compress bunch [m]
    sigz = std(z);				% re-calc bunch length for next possible pass through wake calculations
  end
  if abs(cod) == 7 				% BUNCH COMPRESSION CHICANE (R56/m, E/GeV [T566=-1.5*R56, U5666=2*R56])
    ecuts    = 0;					% no dE/E cuts shown on plots
    zcuts    = 0;					% no Z cuts shown on plots
    iswake   = 0;					% turn off induced voltage plot
    R56      = beamline(j,2);       % R56 value [m]
    dR56_R56 = beamline(j,4);       % relative R56 jitter 
    eps      = dR56_R56;            
    if R56>0
      error('R56 for chicane is always <0...  quitting')
    end
    if dR56_R56>0
      disp('Switched-on jitter on R56 in chicane')
    end
    T566 = -1.5*R56;            % T566 value [m]
    U5666=  2.0*R56;            % U5666 value [m]
    E56  = beamline(j,3);       % Nominal energy of compressor [GeV]
    if E56 < 0.020				% need positive, reasonable nominal R56-energy [GeV]
      error(['Chicane section (7) of R56=%7.4f m has nominal-energy too small (not ultra-relativistic) in ' fnfm],R56)
    end
    dd   = (E-E56)/E56;			% relative energy error w.r.t. nominal compressor energy
    z    = R56*(1-eps)*dd + T566*(1-eps)*dd.^2 + ...    % modified by P. Craievich 02/05/06 (relative R56 jitter)
           U5666*(1-eps)*dd.^3 + z ...
           + eps*R56/2;		    % compress or anti-compress bunch [m] 
    sigz = std(z);				% re-calc bunch length for next possible pass through wake calculations
  end
  
  if abs(cod) == 22			    % INCOHERENT ENERGY SPREAD ADDITION
    ecuts  = 0;				    % no dE/E cuts shown on plots
    zcuts  = 0;				    % no Z cuts shown on plots
    iswake = 0;				    % turn off induced voltage plot
    id_rms = beamline(j,2);     % rms incoherent relative energy spread to be added in quadrature [ ]
	d = d + id_rms*randn(length(d),1);	% incread dE/E by the incoherent addition [ ]
    E = Ebar*(1 + d);			% load energy array [GeV]
  end

  
  %ii = find(z);
  %z_bar = 1E3*mean(z(ii));	% mean z-pos AFTER CUTS [mm]
  
  if cod < 0 || cod == 99    % plot after each negative code point in beamline
                                                % end nargout < 1 (i.e., end plots section)
      jc = jc + 1;                                  % count output locations (where cod < 0 | cod == 99)
      dE_Ej(1:length(d),jc) = d(:);
      zposj(1:length(z),jc) = z(:);
      Ebarj(jc) = Ebar;
      %if nargout > 3				                % if FWHM output parameters wanted...
        zmm     = z*1E3;				            % convert to [mm]
        dpct    = d*100;				            % convert to [%]
		ii = find(zmm);
        [Nz1,Z] = hist(zmm(ii),Nbin);			    % bin the Z-distribution
        HIST_Z(1:Nbin,jc) = Nz1;                           % I'll take that too
        AXIS_Z(1:Nbin,jc) = Z;                             % With the axis
        ZFWmmj(jc) = FWHM(Z,Nz1,0.5);		        % calc. Z-FWHM [mm]
        dZ = mean(diff(Z));    	           	        % Z bin size [mm]
        I  = Nz1*(Ne1/1E10/Nesim)*elec*cspeed/dZ;	% convert N to peak current [kA]
        if gzfit
          [If,q,dq] = gauss_fit(Z,I,1E-3*ones(size(I)),0);
          sigzGj(jc) = q(4);				        % gaussian fit sigma_Z [mm]
          I_pkfj(jc) = max(If);
        end
 		ii = find(dpct);
        [Nd,D]  = hist(dpct(ii),Nbin);			    % bin the dE/E-distribution
        HIST_D(1:Nbin,jc) = Nd;                            % I'll take that too
        AXIS_D(1:Nbin,jc) = D;                             % With the axis
        dFWpctj(jc) = FWHM(D,Nd,0.5);		        % calc. dE/E-FWHM [%]
        %z_barj(jc) = z_bar;
        Ebarcutsj(jc) = Ebarcuts;
        I_pk1 = max(I);
        I_pkj(jc) = I_pk1;
        if gdfit
          [yf,q,dq] = gauss_fit(D,Nd,1E-3*ones(size(Nd)),0);
          sigEGj(jc) = q(4);				        % gaussian fit sigma_dE/E0 [%]
        end  
        %fcutj(jc)    = 1 - Ne1/Ne0;					% fraction of particles cut [ ]
        fcutj(jc) = Nesim;
      %end
  end                                               % end cod < 0 | cod == 99 stuff
  

end                                                 % end loop over all beamline section

if nargout == 1
    LT_OUTPUT.Z.DIST = zposj;
    LT_OUTPUT.Z.HIST = HIST_Z;
    LT_OUTPUT.Z.AXIS = AXIS_Z;
    %LT_OUTPUT.Z.AVG  = z_barj;
    LT_OUTPUT.Z.FWHM = ZFWmmj;
    if gzfit
        LT_OUTPUT.Z.SIG  = sigzGj;
        LT_OUTPUT.I.SIG  = I_pkfj;
    end

    LT_OUTPUT.E.DIST = dE_Ej;
    LT_OUTPUT.E.HIST = HIST_D;
    LT_OUTPUT.E.AXIS = AXIS_D;
    LT_OUTPUT.E.AVG  = Ebarj;
    LT_OUTPUT.E.FWHM = dFWpctj;
    if gdfit
        LT_OUTPUT.E.SIG  = sigEGj;
    end

    LT_OUTPUT.I.PART = fcutj;
    LT_OUTPUT.I.PEAK = I_pkj;

end

end_time = get_time;
disp(' ')
disp(['LiTrack ended:' end_time])

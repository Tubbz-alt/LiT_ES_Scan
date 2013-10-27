  global PARAM; % Initial beam and machine parameters  
  
  % FACET energy set points. We are pretty sure about these, I think. . .
  E0        = PARAM.ENRG.E0;    % GeV ... initial energy
  E1        = PARAM.ENRG.E1;    % GeV ... energy at LBCC
  E2        = PARAM.ENRG.E2;    % GeV ... energy at FACET
  
  % NRTL compressor klystron and R56 #s
  NRTL_ampl = PARAM.NRTL.AMPL;  % AMPL DR13 11 VDES
  NRTL_phas = PARAM.NRTL.PHAS;  % on the zero-crossing
  NRTL_leff = PARAM.NRTL.LEFF;  % cavity length
  NRTL_R56  = PARAM.NRTL.R56;   % This is design val
  NRTL_T566 = PARAM.NRTL.T566;  % Design val?
  NRTL_ELO  = PARAM.NRTL.ELO;   % NRTL low energy cut
  NRTL_EHI  = PARAM.NRTL.EHI;   % NRTL high energy cut
  
  % Phase and length of 02-10
  LONE_leff = PARAM.LONE.LEFF;  % Length of LI02-LI10 (m)
  LONE_phas = PARAM.LONE.PHAS;  % Chirp phase
  LONE_gain = PARAM.LONE.GAIN;  % energy gain in LI02-LI10 (GeV)
  LONE_ampl = PARAM.LONE.FBAM;  % feedback amplitude (GV)
  
  % S10 chcn #s
  LI10_R56  = PARAM.LI10.R56;   % Measured val?
  LI10_ISR  = PARAM.LI10.ISR;   %
  LI10_ELO  = PARAM.LI10.ELO;   % S20 low energy cut
  LI10_EHI  = PARAM.LI10.EHI;   % S20 high energy cut
  
  % Energy gain and length of 02-10
  LTWO_leff = PARAM.LTWO.LEFF;  % Length of LI02-LI10 (m)
  LTWO_phas = PARAM.LTWO.PHAS;  % Chirp phase
  LTWO_ampl = PARAM.LTWO.FBAM;  % feedback amplitude (GV)
  
  % S20 chcn R56 #s
  LI20_NLO  = PARAM.LI20.NLO;   % low energy notch edge
  LI20_NHI  = PARAM.LI20.NHI;   % high energy notch edge
  LI20_R56  = PARAM.LI20.R56;   % Measured val?
  LI20_T566 = PARAM.LI20.T566;  % Measured val?
  LI20_ISR  = PARAM.LI20.ISR;   %
  LI20_ELO  = PARAM.LI20.ELO;   % S20 low energy cut
  LI20_EHI  = PARAM.LI20.EHI;   % S20 high energy cut
  
  % 6mm bunches coming out of the ring teensy energy spread.
  sigz0     = PARAM.INIT.SIGZ0;	% rms bunch length used when inp=G or U above [m]
  sigd0     = PARAM.INIT.SIGD0;	% rms relative energy spread used when inp=G or U above [ ]
  z0_bar    = PARAM.INIT.Z0BAR; % axial offset of bunch [m] (used also with file input - mean of file removed first)
  d0_bar    = PARAM.INIT.D0BAR; % relative energy offset of bunch [ ]  (used also with file input - mean of file removed first)
  
  % 200K sim particles = 100K electrons per sim particle
  Nesim     = PARAM.INIT.NESIM;	% number of particles to generate for simulation when inp=G or U (reasonable: ~1000 to ~100000)
  Ne        = PARAM.INIT.NPART; % number of particles initially in bunch
  
  % The Holtzapple skew. Someday they'll name a skew after me. . .
  asym      = PARAM.INIT.ASYM;	% for inp='M' or 'G': sets rise/fall time width (-1<asym<1)
  cut       = PARAM.INIT.CUT;   % for inp='G': sets rise/fall time width (0.5<=cut<inf)

  fb_off = 0;
  fb_on  = 1;
  sband  = 1;
  NRTL_U5666 = 0;
  LI20_U5666 = 0;
  no_quad = 1;

beamline = [
       11		NRTL_ampl      NRTL_phas     sband      fb_off      NRTL_leff     % Compressor cavity AMPL DR13 13 VDES
       26	    NRTL_ELO       NRTL_EHI      0          0           0             % Approximate energy acceptance of NRTL
        6		NRTL_R56       NRTL_T566     NRTL_U5666 E0          0             % Design NRTL ~0.603, BDES to KMOD for E-164 gives 0.588
       11		E1             LONE_phas     sband      fb_on       LONE_leff     % 2-6, nominal 9GeV, no feedback
        6	    LI10_R56       0             0          E1          no_quad       % 2nd half of the chicane. Design was -0.0745, as built -0.076
       22       LI10_ISR       0             0          0           0             % Approximate SR growth in E-spread from chicane
       26       LI10_ELO       LI10_EHI      0          0           0             % Momentum Slits in FACET
       37		cut            0             0          0           0             % Clip any rediculously long tails
       11       E2             LTWO_phas     sband      fb_on       LTWO_leff     % Boost to 20.35
       28       LI20_NLO       LI20_NHI      0          0           0             % Notch collimate
       26       LI20_ELO       LI20_EHI      0          0           0             % Momentum Slits in FACET
       6		LI20_R56       LI20_T566     LI20_U5666 E2          0             % FACET 'dogleg' like chicane
       22       LI20_ISR       0             0          0           0             % Approximate SR growth in E-spread from dogleg
       37       cut            0             0          0           0             % Clip any rediculously long tails
       ];
   


% Sign conventions used:
% =====================
%
% phase = 0 is beam on accelerating peak of RF (crest)
% phase < 0 is beam ahead of crest (i.e. bunch-head sees lower RF voltage than tail)
% The bunch head is at smaller values of z than the tail (i.e. head toward z<0)
% With these conventions, the R56 of a chicane is < 0 (and R56>0 for a FODO-arc) - note T566>0 for both
% 
% * = Note, the Energy-window cut (25-code) floats the center of the given FW window in order center it
%     on the most dense region (i.e. maximum number of particles).
% **  1:=1st wake file (e.g. S-band) is used, 2:=2nd wake file (e.g. X-band)
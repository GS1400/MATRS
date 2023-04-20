%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% initTune.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initTune provides two data tuning struture for MATRS
%
function [ctune,itune] = initTune(ctune,itune,nc,ni)
% artificial lower and upper bounds to help prevent overflow
% initialize structure containing all tuning parameters 
%
if nc>0
    % ctune % structure containing all tuning parameters 
    %       % all parameters have a default that can be overwritten 
    %       % by specifying it as an input
    if ~exist('ctune'), ctune=[]; end
    % number of mutation points
    if ~isfield(ctune,'clambda'), ctune.clambda =max(6,nc); end
    % number of selected mutation points
    if ~isfield(ctune,'cmu'), ctune.cmu = 3+ceil(log(nc)); end
    % initial recombination step size
    if ~isfield(ctune,'sigmac'), ctune.sigmac = 1; end
    % parameter for trust region condition
    if ~isfield(ctune,'czeta'), ctune.czeta = 1e-20; end
    % parameter for updating line search step sizes
    if ~isfield(ctune,'cnu'), ctune.cnu = 2; end
    % parameter for updating trust region radii
    if ~isfield(ctune,'ctheta'), ctune.ctheta = 2; end
    % maximum value for sigma
    if ~isfield(ctune,'csigmamax'), ctune.csigmamax = 1e10;  end
    % minimum value for sigma
    if ~isfield(ctune,'csigmamin'), ctune.csigmamin = 1e-10;  end
    % parameter for cusequnce
    if ~isfield(ctune,'ckmax'), ctune.ckmax = 20;  end
    % factor for adjusting Y
    if ~isfield(ctune,'gammaX'), ctune.gammaX = 1e3; end
    % factor for adjusting gradient
    if ~isfield(ctune,'gammav'), ctune.gammav = 1e2; end
    % minimum threshold for Delta
    if ~isfield(ctune,'cDeltamin'), ctune.cDeltamin = 1e-3; end
    % maximum value for the initial Delta
    if ~isfield(ctune,'Delmax'), ctune.Delmax = 1;  end
    % factor for combination direction
    if ~isfield(ctune,'sc'), ctune.sc = 10;  end
    % parameter for restarting affine scaling matrix
    if ~isfield(ctune,'mmax'), ctune.mmax = 5;  end
    % maximum number of points wanted for cusequnce
    if ~isfield(ctune,'N'), ctune.N = 10000;  end
end
if ni>0
    % itune % structure containing all tuning parameters 
    %       % all parameters have a default that can be overwritten 
    %       % by specifying it as an input
    if ~exist('itune'), itune=[]; end
    % number of mutation points
    if ~isfield(itune,'ilambda'), itune.ilambda =max(6,ni); end
    % number of selected mutation points
    if ~isfield(itune,'imu'), itune.imu = 3+ceil(log(ni)); end
    % initial recombination step size
    if ~isfield(itune,'sigmai'), itune.sigmai = 1; end
    % parameter for the trust region condition
    if ~isfield(itune,'izeta'), itune.izeta = 1e-20; end
    % parameter for updating line search step sizes
    if ~isfield(itune,'inu'), itune.inu = 2; end
    % parameter for updating trust region radii
    if ~isfield(itune,'itheta'), itune.itheta = 2; end
    % maximum try to find feasible point in iTRS
    if ~isfield(itune,'stuckmax'), itune.stuckmax = 10;  end
    % upper bound on Delta to use Delta=Delta-1 to update radius in iTRS
    if ~isfield(itune,'iDeltabar'), itune.iDeltabar = 3;  end
    % maximum value for sigma
    if ~isfield(itune,'isigmamax'), itune.isigmamax = 1e2;  end
    % parameter for iusequence
    if ~isfield(itune,'ikmax'), itune.ikmax = 30;  end
    % maximum value for the initial Delta
    if ~isfield(itune,'Delmax'), itune.Delmax = 30;  end
    % minimum value for the initial Delta
    if ~isfield(itune,'Delmin'), itune.Delmin = 10;  end
    % parameter for restarting affine scaling matrix
    if ~isfield(itune,'mmax'), itune.mmax = 5;  end
    % factor for combination direction
    if ~isfield(itune,'sc'), itune.sc = 5;  end
    % factor for adjusting Y
    if ~isfield(itune,'gammaX'), itune.gammaX = 1e3; end
    % factor for adjusting gradient
    if ~isfield(itune,'gammav'), itune.gammav = 1e2; end
    % maximum number of points wanted for cusequnce
    if ~isfield(itune,'N'), itune.N = 10000;  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



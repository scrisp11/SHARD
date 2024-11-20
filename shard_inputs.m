%% Input for shard_main_noPlots
% clear all;      close all;
plots           = true;     set(0, 'DefaultAxesFontSize', 6.5);
plotParticleDists = true;
saveTag         = false;
animation       = false;

% bOff = 0;
% fName           = ['shard_buncher_opt_gamma12p_L5_noTaper_bOff_',num2str(floor(bOff*1e3))];
% fName           = ['shard_buncher_opt_gamma17p_L5_noTaper_bOff_',num2str(floor(bOff*1e3))];
% fName           = 'shard_Pond_gamma12_tapered_fig1';
% fName           = 'shard_mpPond_gamma12_untapered_fig1';

%% Structure
struc.c1        = 0.1;     % Ideal
struc.c2        = 0.1;
% struc.c1        = 0.025; % data
% struc.c2        = 0.025;
struc.zstop     = 5e-3;
% struc.zstop     = 1.28e-3;
struc.gap       = 400e-9;   % Ideal
% struc.gap       = 1200e-9;% data
struc.lambda    = 780e-9;


 %% Electrons
elec.kin          = (12 - 1)*511e3; % Exp: gamma 12, gamma 17.14
 
elec.emit         = 200e-9;  % Normalized
elec.deltagamma   = 0.0005;  % gam0*()
elec.sigmax       = 30e-6;  

elec.emit         = 1e-9;
elec.deltagamma   = 0*0.0005;  % gam0*()
elec.sigmax       = 1e-6;  

bunching        = false;     % Prebunch beam

elec            = calcElec(elec);
[struc,Ezfun]   = calcStruc(struc, elec); % (need beta0)
%% Simulation
elec.n          = 500;
nfreq           = 2^7;                                      % Number of Spatial Harmonics

opt             = true;                                   % Removes Particles as they hit the walls

nstep           = floor(struc.zstop/struc.lambda/25)*4;    % zstop/nstep should be an integer multiple of lambda
stepsize        = struc.zstop/(nstep-1);
z               = 0:stepsize:struc.zstop;
%% Inputs:
las.sigma       = 5e-3;

% Ponderomotive:
las.thI         = 0e-04; % angle incidence
las.E0          = 10e9;
taperOn         = true;
bOff            = 0;
pondTag         = 1;
constPsi        = 30; % resonant phase
apfTag          = 0;
phiSmooth       = false;
lambdap         = [linspace(400e-6,600e-6,sum(z<1e-3)),linspace(600e-6,860e-6,sum(z>=1e-3))];
buncher         = true;

% Flat:
% las.thI         = -3.1e-04; % angle incidence = mean(tap(z<2) - tap(1)) when taperOn = true 
% las.E0          = 2.4e9;
% taperOn   = false;
% bOff            = 0;
% pondTag         = 0;
% apfTag          = 0;
% constPsi        = 60; % resonant phase
% hardCodeBuncher = true;

% Tapered
% las.thI         = 0e-04; % angle incidence = mean(tap(z<2) - tap(1)) when taperOn = true 
% las.E0          = 1e9;
% taperOn   = false;
% % bOff            = 2*pi;
% pondTag         = 0;
% apfTag          = 0;
% constPsi        = 90; % resonant phase
% phiSmooth       = true;
% Ap              = 0; lambdap = 1;
% las.E0          = 1e9;
% 
% elec.kin          =(17.14-1)*511e3;
% st              = 0.64e-3;
% dr              = 0;
% elec.kin          =(12-1)*511e3;
% st              =  1.7e-3;
% dr              = 1.28e-3;

% APF:
% las.thI         = 0e-04; % angle incidence
% las.E0          = 10e9;
% taperOn         = true;
% bOff            = 0*pi/8;
% pondTag         = 0;
% apfTag          = true;
% apfFrac         = 3;
% ratio           = 0.5;
% constPsi        = 90 - apfTag*pi/apfFrac*180/pi; % resonant phase
% % % lambdap         = [ones(1,sum(z<1e-3)),linspace(600,600,sum(z>=1e-3))].*struc.lambda;
% % % lambdap         = [ones(1,sum(z<1e-3)),linspace(600,1400,sum(z>=1e-3))].*struc.lambda;
% % lambdap         = [ones(1,sum(z<1e-3)),linspace(600,1600,sum(z>=1e-3))].*struc.lambda;
% phiSmooth       = false;
% % lambdap         = [ones(1,sum(z<1e-3)),linspace(842,1775,sum(z>=1e-3))].*struc.lambda;
% buncher        = true;
% load(['LUT_APFstable_psiRes_',num2str(pi/apfFrac*1e4,'%3.0f')])
% LUTsubset = LUT(betagamma>=elec.gam0.*elec.beta0);
% betagammaSub = betagamma(betagamma>=elec.gam0.*elec.beta0);
% apfFit = polyfit(betagammaSub,LUTsubset,1);
% lambdap1        = elec.gam0.*elec.beta0*apfFit(1)+apfFit(2);
% lambdap2        = 16*apfFit(1)+apfFit(2);
% % lambdap         = [ones(1,sum(z<1e-3)),linspace(1200,1200,sum(z>=1e-3))].*struc.lambda;
% lambdap         = [ones(1,sum(z<1e-3)),linspace(lambdap1,lambdap2,sum(z>=1e-3))].*struc.lambda;

% Buncher
% las.thI         = -0*2.8e-04; % angle incidence
% las.E0          = 1e9;
% taperOn   = false;
% % bOff            = pi;
% pondTag         = 0;
% apfTag          = 0;
% constPsi        = 90; % resonant phase
% lambdap         = 1.5e-3*1/1.65;
% E0end           = 1e-3;
%Smoothing
sigma = 10e-6/stepsize;                             % pick sigma value for gaussian PSF of optical system

%% Laser
las.lambda      = 780e-9; las.delT = 100e-15;
las.k           = 2*pi/las.lambda;

% Phase
res             = (-pondTag)*ones(1,nstep);
psi_res         = -pi/180*constPsi*ones(1,nstep);

st = 0; dr = 0; % Default to no buncher
if apfTag % APF Based
    if buncher 
        st              = 0.06e-3; 
        dr              = 1.2e-3 - st;
    end
        psi_res(z<(dr+st)) = -pi;
        psi_res(z<(st))     = 0;
        psi_res(z>=(st+dr)&mod(z-dr-st,lambdap)./lambdap<ratio) = -pi - psi_res(z>=(st+dr)&mod(z-dr-st,lambdap)./lambdap<ratio);
elseif pondTag % Ponderomotive Based:
    if buncher
        st = 0.05e-3;     dr              = 1e-3 - st;
    end
    res(z<=st)             = 0;
    psi_res(z<=(st)) = 0;
elseif buncher
        st = 0.1; dr = 1e-3 - st;
end

LTot            = struc.zstop; % limits laser end

% Amplitude Ponderomotive Accel
captL           = dr + st; bL = captL;
Ap              = zeros(1,nstep);
Ap(z<captL)     = 0;
Ap(z>=captL)    = linspace(0.55,0.3,sum(z>=captL));
% Ap(z>=bL&z<=5e-3)    = linspace(0.05,0.35,sum(z>=bL&z<=5e-3));
% Ap(z>=5e-3)    = linspace(0.35,0.05,sum(z>=5e-3));

% Calculation
psi_res(res==-1)= psi_res(res==-1) - pi;          % resonant phase for ponderomotive case
las.phi     = psi_res + pondTag*Ap.*cos(2*pi.*(z- st-dr)./lambdap)-res.*2*pi.*(z- st-dr)./lambdap;      % Ponderomotive

las.phi(z>(st)) = las.phi(z>(st)) + bOff; phi0 = las.phi(1);
% las.phi         = unwrap(mod(las.phi,2*pi)); %make phase continuous
% las.phi         = las.phi - las.phi(1) + phi0;

% Amplitude
las.G_gauss     = las.E0*exp((-(z-struc.zstop/2).^2)/(2*las.sigma^2));
bR              = (z>st) & (z<(st+dr));
aR              = (z>=(st+dr));
las.G_gauss(bR) = las.G_gauss(bR)/100;
las.G_gauss(z>LTot) = 0;

%  Spatial resolution filtering
sigma = 10e-6/stepsize;                             % pick sigma value for gaussian PSF of optical system
gaussFilter = gausswin(round(6*sigma + 1))';
gaussFilter = gaussFilter / sum(gaussFilter);       % normalization
las.G_gauss1 = conv(las.G_gauss, gaussFilter, 'same');     % Apply filter
las.G_gauss1(1:10) = las.G_gauss(1:10);
las.G_gauss1(end-9:end) = las.G_gauss(end-9:end);
las.G_gauss =  las.G_gauss1;

if phiSmooth
    sigma = 10e-6/stepsize;                             % pick sigma value for gaussian PSF of optical system
gaussFilter = gausswin(round(6*sigma + 1))';
gaussFilter = gaussFilter / sum(gaussFilter);       % normalization
    phi1 = las.phi(end);
    las.phi = conv(las.phi, gaussFilter, 'same');
    las.phi(1:ceil(length(gaussFilter/2))) = phi0;
    las.phi(end - ceil(length(gaussFilter/2)) + 1:end) = phi1;
end
 
%%
plotGPhi(z,las.G_gauss*1e-6,las.phi,1,false);

%% Calculates phase to match accelerating electrons (or not(taperOn))
[taper, las, gam_res]  = calcTaper(elec,struc,las,z,taperOn);

%% Generate Particles
[yold,frac]     = genParts(elec,psi_res, bunching,struc);
% yold(elec.n+1:2*elec.n) = (pi + yold(elec.n+1:2*elec.n))/3 - 0.7;
drawnow;

%% Run Simulation
shard_main_noPlots; 

%% Output: gammap,thetap,x,xp, focusTrack
capt            = 0.511*(gammap(end,:) - elec.gam0)>(0.9*0.511*(gam_res(end)-elec.gam0));
focus           = logical(focus); focused = sum(focus);

capt            = capt&focus; accelFoc = sum(capt);
target          = sum(capt.*focus*(gam_res(end)-elec.gam0));
norm            = npart*(gam_res(end)-elec.gam0);
target          = target/norm;
if sum(capt)==0; capt=focus;end
inow            = nstep;
clear captured rejected cap_phase rej_phase

if ~opt
        capturedGm = gammap(inow,:);
        rejectedGm = gammap(inow,:);
        capturedGm(~focus) = NaN;
        rejectedGm(focus) = NaN;

        capturedTp = thetap(inow,:);
        rejectedTp = thetap(inow,:);
        capturedTp(~focus) = NaN;
        rejectedTp(focus) = NaN;
        
        cap_phase = xp(inow,:);
        rej_phase = xp(inow,:);
        cap_phase(~focus) = NaN;
        rej_phase(focus) = NaN;
    else
        capturedGm = gammap;
        rejectedGm = NaN;
        capturedTp = thetap;
        rejectedTp = NaN;
        cap_phase = xp;
        rej_phase = NaN;
    end
disp(['Focused Particles = ' num2str(sum(focus))])
disp(['Total captured particles =' num2str(sum(capt.*focus))])
disp(['Target cost function =' num2str(target)])
disp(['Phase Offset = ' num2str(bOff)]);
disp(['Resonant Energy Gain (in keV) = ' num2str((gam_res(end)-elec.gam0)*511)])
disp(['Maximum Energy Gain (in keV) = ' num2str(max(gammap(end,focus)-gammap(1,focus))*511)])
disp(['Maximum Energy Loss (in keV) = ' num2str(min(gammap(end,focus)-gammap(1,focus))*511)])
%%
if saveTag
    if ~opt
        save(fName,'struc','gam_res','DATA','focus','bOff','gammap','thetap','elec','las','taper','Ap','lambdap')
    else
        save(fName,'struc','gam_res','focus','bOff','gammap','thetap','elec','las','taper','Ap','lambdap')
    end
end
%%
if plots
    shard_quickPlots;
end
%%
function elec = calcElec(elec)
    % Must define elec.kin, elec.deltagamma
    elec.E0         = elec.kin+511e3;
    elec.gam0       = elec.E0/511e3;
    elec.beta0      = sqrt(1-1/elec.gam0^2);
    elec.deltagamma = elec.gam0*elec.deltagamma;
end

function [struc, Ezfun] = calcStruc(struc, elec)
    % Must define struc.lambda,struc.c1,c2, elec.beta0
    struc.Gamma     = 2*pi*sqrt(1/struc.lambda^2-(1/(struc.lambda/elec.beta0))^2);
    struc.yc        = log(abs((struc.c1)./(struc.c2)))./(2*struc.Gamma);
    struc.k         = 2*pi/struc.lambda;
    Ezfun           = @(c,x)((c(1)*exp(-struc.Gamma*x)+c(2)*exp(struc.Gamma*x)));
    struc.sf        = abs(Ezfun([struc.c1,struc.c2],struc.yc));
end

%%
function plotGPhi(z,G,phi,fig,holdTag)
    figure(fig);
    if holdTag; hold on;
    else; clf; end
    yyaxis left
    plot(z*1e3,G,'LineWidth',1.5);
    xlabel('z (mm)');
    ylabel('G (MeV/m)')
    yyaxis right
    plot(z*1e3,phi,'LineWidth',1.5)
    ylabel('\phi (rad)')
end
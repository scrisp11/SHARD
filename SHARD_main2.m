%% DLA Simulation for SLM Phase Profile SHARD v 2022
% clear all
% close all

% Shard Version 2022. Arbitrary on axis phase and amplitude profile
% are Fourier-analyzed to generate a set of spatial harmonics to describe
% the field in the DLA gap

%% Input Parameters
npart   = elec.n;             % Number of particles
E0_kin  = elec.E0_kin;       %5 MeV UCLA Pegasus
gam0    = elec.gam0;
beta0   = elec.beta0;
G0      = las.G0;           % Gradient
lambda  = las.lambda;
k       = las.k;
alpha   = G0 / 511e3 / k;

zstop   = las.lambda*5000;            %4mm total length

lambdap = params(1)*lambda;     %Periodicity Accel

La = params(3);
psi_bunch = 0;
Ld = (gam0)^2/(G0*La/E0_kin)*lambda/(2*pi);
bunchLength = La + Ld;
nstep = zstop/lambda/10;

stepsize = zstop/(nstep-1);
z = 0:stepsize:zstop;
sigmaBlur   = blurVar/stepsize;   % pick sigma value for gaussian PSF of optical system


%% Define Phase and Amplitude Mask 
psi_res = params(2);            % Choose resonant phase

phi=zeros(1,nstep)+psi_res;
% phi(mod(z+lambdap/2,lambdap)/lambdap<0.5) = -pi-phi(mod(z+lambdap/2,lambdap)/lambdap<0.5);  % APF-simple
phi(mod(z-lambdap/3,lambdap)/lambdap<0.5) = -pi-phi(mod(z-lambdap/3,lambdap)/lambdap<0.5);  % APF-simple

% 
% load('betaMax.mat')
% X = [ones(length(gamVec'),1) gamVec'];
% [minBetaMax,I] = min(beta_MaxFilt);
% b = X\NpVec(I)';
% npFunc0 = @(gamma1)(b(1)+b(2)*gamma1); % approx the optimal period for a given energy

npFunc = @(gamma1)((gamma1-gam0)*params(4)+params(1));
ptag = 1; phi = []; 
gamN = gam0;
nNext = floor(npFunc(gam0)/2*lambda*nstep/zstop);
xn = nNext*zstop/nstep;
while length(phi)<nstep
    if ptag < 0
        phi_next = zeros(1,nNext)-pi-psi_res;
    else
        phi_next = zeros(1,nNext)+psi_res;
    end
    phi = [phi phi_next];
    ptag = -ptag;
    gamN = gamN+0.9*G0/511e3*nNext*zstop/nstep;
    nNext = floor(npFunc(gamN)/2*lambda*nstep/zstop);
    xn = xn + nNext*zstop/nstep;
end

% phi = phi(1:length(phi));

phi_bunch = psi_bunch*ones(1,sum(z<1.0*bunchLength));
phi = [phi_bunch phi];
phi = phi(1:length(z));
% phi(z<bunchLength*10^-3) = psi_bunch;

% Amplitude
sigma = 10/sqrt(2);  % Field amplitude is modeled as a gaussian with sigma
G_gauss = G0*exp((-(z-zstop/2).^2)/(2*sigma^2));

G_gauss = G_gauss.*(1-fAmp(z,params(3))+fAmp(z,bunchLength));

G_crop      = [];
phi_crop    = [];
iBunch      = sum(z<=bunchLength);
if 1 == accel(1) %There is a buncher
    G_crop      = G_gauss(1,1:iBunch);
    phi_crop    = phi(1,1:iBunch);
end
if 1 == accel(2) %There is an accelerator
    G_crop      = [G_crop   G_gauss(1,iBunch+1:end)];
    phi_crop    = [phi_crop phi(1,iBunch+1:end)];
end
G_gauss = G_crop; phi = phi_crop;
nstep = length(G_gauss);
zstop = nstep*stepsize;
z = 0:zstop/(nstep-1):zstop;

% Phase and amplitude of mask on same plot
if plots
    figure(26)
    yyaxis left
    plot(z, abs(G_gauss.*1e-9),'LineWidth',1.5)
    ylabel('E, GV/m')
    yyaxis right
    plot(z,phi)
    ylabel('Phase, rad')
    xlabel('z, m')
    title('Electric Field Strength')
end


%% Acceleration Gradient and resonant phase tapering
gam_res=zeros(1,nstep);
beta_res=zeros(1,nstep);
tap=zeros(1,nstep);
phi_add = zeros(1,nstep);

gam_res(1)=gam0;
beta_res(1)=sqrt(1-1/gam0^2);
tap(1)=1/beta_res(1)-1;

for i=2:nstep
gam_res(i)=gam_res(i-1)-G_gauss(i)/511e3*stepsize*sin(phi(i));                               % APF case
%gam_res(i) = gam_res(i-1)-G_gauss(i)/511e3.*besselj(res,Ap).*stepsize.*sin(psi_res);        % Ponderomotive case
beta_res(i)=sqrt(1-1./gam_res(i).^2);
tap(i)=1/beta_res(i)-1;
phi_add(i) = trapz(z(2)-z(1),k.*tap);
end

gradient_GeVperm = (gam_res(end)-gam_res(1))*0.511/zstop*0.001;
maxenergy_MeV = gam_res(end)*0.511;
table(gradient_GeVperm,maxenergy_MeV)

phi = phi + phi_add;

%%  Spatial resolution filtering
if blur
    phi = blurF(sigmaBlur,phi);
end
%% Spatial Harmonics Fourier Analysis
E_axis = G_gauss.*exp(1i.*phi);
dwnsamplerate = max(floor(nstep/nfreq),1);
E_axis_dn = downsample(E_axis,dwnsamplerate);
Efft = fftshift(fft(E_axis_dn));
[~,idx] = max(abs(Efft));
nh = size(Efft,2);
kmax = pi/(zstop/nstep*dwnsamplerate);
%kmax=0.5/lambda/dwnsamplerate;   
deltak = kmax/(nh/2);  
kshift = (-nh/2:nh/2-1)*deltak; % zero-centered frequency range

%kntap = zeros(nh,1); betantap = kntap; gammantap = kntap;
for jh = 1:nh
kntap(jh) = k + kshift(jh);    
betantap(jh) = k ./ kntap(:,jh);
gammantap(jh) = sqrt(1./(1-betantap(:,jh).^2));
end

phi_tap.nh          = nh;
phi_tap.kntap       = kntap;
phi_tap.betantap    = betantap;
phi_tap.Efft        = Efft;
phi_tap.linearizedfields = linearizedfields;
phi_tap.aperture    = aperture;
phi_tap.G_gauss     = G_gauss;
phi_tap.phi         = phi-phi_add;
phi_tap.z           = z;

%% Particle generation
if genPtTag
    yold = genParts(las,elec,phi_tap,twissTag, R56Tag);
end
%% Tracking
focus = ones([1,npart]);
DATA(1,:)=yold;

tic
if opt
    for i=1:nstep
        [ynew, focus_new]=dla_push_RK4(z(i),yold,phi_tap,k,zstop/nstep,focus);
        filt0 = (1==(repmat(focus_new,1,4)));
        ynew = ynew(filt0~=0);
        yold=ynew;
        focus = ones([1,length(yold)/4]);
    end
else
    for i=1:nstep
        [ynew, focus_new]=dla_push_RK4(z(i),yold,phi_tap,k,zstop/nstep,focus);
        DATA(i+1,:)=ynew;
        yold=ynew;
        focus = focus_new;
    end
end
toc

%% Outputs
gammap = []; thetap = []; x = []; xp = [];
if opt
    npartN = length(yold)/4;
    gammap= real(ynew(1:npartN));
    thetap= real(ynew(npartN+1:2*npartN));
else
    DATA = real(DATA);
    for i = 1:nstep
        gammap(i,:)= DATA(i,1:npart);
        thetap(i,:)= DATA(i,npart+1:2*npart);
        x(i,:)= DATA(i,2*npart+1:3*npart);
        xp(i,:)= DATA(i,3*npart+1:4*npart);
    end
end

inow = nstep;
capt = gammap(end,:)>gam_res(end)*0.9;
focus = logical(focus); focused = sum(focus);
accelFoc = sum(capt.*focus);
if opt
    target = -focused;
else
    target = -sum(capt.*focus.*(gam_res(inow)-gam0));
end

table(gradient_GeVperm,maxenergy_MeV,target, focused, accelFoc)
%%
function phi = blurF(sigma,phi)
    gaussFilter = gausswin(round(6*sigma + 1))';
    gaussFilter = gaussFilter / sum(gaussFilter);       % normalization
    phi = conv(phi, gaussFilter, 'same');              % Apply filter
end
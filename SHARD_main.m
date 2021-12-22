%% DLA Simulation for SLM Phase Profile SHARD v 2022
clear all
close all

% Shard Version 2022. Arbitrary on axis phase and amplitude profile
% are Fourier-analyzed to generate a set of spatial harmonics to describe
% the field in the DLA gap


%% Input Parameters
plots = true; %Set false to remove plots (for optimization loops) 
bunching = false;
linearizedfields = false;

E0_kin = 5e6;       %5 MeV UCLA Pegasus
G0 = 2e9;           % Gradient
lambda = 0.80e-6;

npart = 1e3;             % Number of particles
zstop = lambda*6000;       % 4.8mm total structure length, lower for code tests
nstep = zstop/lambda/20;    % zstop/nstep should be an integer multiple of lambda
nfreq = 120;               % Number of spatial harmonics


% Ponderomotive case input data
lambdap = 360*lambda;
Ap = 0.12;          
res = -1;
% APF case input data
lambdap = 600*lambda;
Ap = 0;
res = 0;

k = 2*pi/lambda;
E0 =E0_kin+511e3;
gam0 = E0/511e3;
beta0 = sqrt(1-1/gam0^2);
alpha = G0 / 511e3 / k;

stepsize = zstop/(nstep-1);
z = 0:zstop/(nstep-1):zstop;

%% Define Phase and Amplitude Mask 
psi_res = -(pi/180)*45;            % Choose resonant phase
if (res == -1)
    psi_res = psi_res-pi;          % resonant phase for ponderomotive case
end

%phi = phi + Ap.*cos(2*pi.*z/lambdap)+pi/2-res*2*pi.*z/lambdap;      % Ponderomotive

phi=zeros(1,nstep)+psi_res;
%phi(mod(z,lambdap)/lambdap<0.25) = -pi/2;
%phi(mod(z,lambdap)/lambdap<0.5 & mod(z,lambdap)/lambdap>0.25) = -pi-psi_res;
%phi(mod(z,lambdap)/lambdap<0.75 & mod(z,lambdap)/lambdap>0.5) = -pi/2;
phi(mod(z,lambdap)/lambdap>0.5) = -pi-phi(mod(z,lambdap)/lambdap>0.5);  % APF-simple


% Amplitude
sigma = 10/sqrt(2);  % Field amplitude is modeled as a gaussian with sigma
G_gauss = G0*exp((-(z-zstop/2).^2)/(2*sigma^2));

%% Acceleration Gradient and resonant phase tapering
gam_res=zeros(1,nstep);
beta_res=zeros(1,nstep);
tap=zeros(1,nstep);

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
sigma = 10e-6/stepsize;                             % pick sigma value for gaussian PSF of optical system
gaussFilter = gausswin(round(6*sigma + 1))';
gaussFilter = gaussFilter / sum(gaussFilter);       % normalization
%phi = conv(phi, gaussFilter, 'same');              % Apply filter

%% Spatial Harmonics Fourire Analysis
E_axis = G_gauss.*exp(1i.*phi);
dwnsamplerate = max(floor(nstep/nfreq),1);
E_axis_dn = downsample(E_axis,dwnsamplerate);
Efft = fftshift(fft(E_axis_dn));
[~,idx] = max(abs(Efft));
nh = size(Efft,2)
kmax = pi/(zstop/nstep*dwnsamplerate);
%kmax=0.5/lambda/dwnsamplerate;   
deltak = kmax/(nh/2);  
kshift = (-nh/2:nh/2-1)*deltak; % zero-centered frequency range

for jh = 1:nh
kntap(jh) = k + kshift(jh);    
betantap(jh) = k ./ kntap(:,jh);
gammantap(jh) = sqrt(1./(1-betantap(:,jh).^2));
end

phi_tap.nh = nh;
phi_tap.kntap = kntap;
phi_tap.betantap = betantap;
phi_tap.Efft = Efft;
phi_tap.linearizedfields = linearizedfields;

%% Output 
if plots
    figure(31)
    plot(z,gam_res*0.511)
    ylabel('Energy (MeV)');
    xlabel('z [m]')
    
    figure(32)
    plot(z,mod(phi,2*pi))
    xlabel('z [m]')
    title('Electric Field Phase Profile')
    
    figure(33)
    plot(z, abs(G_gauss.*1e-9),'LineWidth',1.5)
    ylabel('E [GV/m]')
    xlabel('z [m]')
    title('Electric Field Strength')
    
    figure(34)
    plot(kshift,abs(Efft))
    title('Spatial Harmonics Content')
    drawnow
end

%% Particle generation
deltagamma = 0.001*gam0;
sigmax = 0.25e-6;
emit = 0.1e-9;

sigmaxp = emit/gam0/beta0/sigmax;
seed = 1;

theta0 = (hammersley(1,npart,seed)*2*pi-pi);
gamma0 = normrnd(gam0,deltagamma,[1,npart]);

% add bunching
if bunching
     Ab = 5*deltagamma;
     R56 = pi/Ab/2;
     gamma0 = gamma0 - Ab.*sin(theta0);
     theta0 = theta0+(gamma0-mean(gamma0))*R56+psi_res;
end

x0 = normrnd(0,sigmax,[1,npart]);
x0p = normrnd(0,sigmaxp,[1,npart]);

%f = 0.0005;                   % Add correlation at entrance of DLA
%x0p = x0p -1/f*x0;

% if plots
%     figure(40)
%     plot(theta0,gamma0,'.')
%     title('Initial Particle Distribution')
%     figure(41)
%     plot(x0,x0p,'.')
%     title('Initial Particle Distribution')
% end

%% Tracking
yold=[gamma0 theta0 x0 x0p];
focus = ones([1,npart]);
DATA(1,:)=yold;
if plots
    tic
end
for i=1:nstep
    [ynew, focus_new]=dla_push_RK4(z(i),yold,phi_tap,k,zstop/nstep,focus);
    DATA(i+1,:)=ynew;
    yold=ynew;
    focus = focus_new;
end
if plots
    toc
end
%% Outputs
DATA = real(DATA);
for i = 1:nstep
gammap(i,:)= DATA(i,1:npart);
thetap(i,:)= DATA(i,npart+1:2*npart);
x(i,:)= DATA(i,2*npart+1:3*npart);
xp(i,:)= DATA(i,3*npart+1:4*npart);
end

inow = nstep;
capt = gammap(end,:)>gam_res(end)*0.9;
focus = logical(focus);
target= sum(capt.*focus*(gam_res(inow)-gam0));
norm = npart*(gam_res(inow)-gam0);

target = target/norm;
clear captured rejected cap_phase rej_phase

if plots
    disp(['Resonant Energy Gain (in keV) = ' num2str((gam_res(end)-gam0)*511)])
    disp(['Accelerated Particles = ' num2str(sum(capt))])
    disp(['Focused Particles = ' num2str(sum(focus))])
    disp(['Total captured particles =' num2str(sum(capt.*focus))])
    disp(['Target cost function =' num2str(target)])
    
%    [f2,xi2] = ksdensity(gammap(inow,:));
    
    %f2=5*f2; %5 is just a plotting factor
    
    captured = gammap(inow,:);
    rejected = gammap(inow,:);
    captured(~focus) = NaN;
    rejected(focus) = NaN;
    
    cap_phase = xp(inow,:);
    rej_phase = xp(inow,:);
    cap_phase(~focus) = NaN;
    rej_phase(focus) = NaN;
    
    figure(20)
    scatter(x(inow,:),xp(inow,:),15,focus)
    xlim([-2e-6, 2e-6])
    ylim([-1e-3, 1e-3])
    
    figure(21)
    subplot(2,1,1)
    scatter(mod(thetap(inow,:),2*pi), gammap(inow,:),15,focus)
    subplot(2,1,2)
    scatter(mod(thetap(inow,:),2*pi), gammap(inow,:),15,capt)
    
    
    figure(1)
    plot(mod(thetap(1,:),2*pi),gammap(1,:),'.','MarkerSize',10)
    hold on
    plot(mod(thetap(inow,:),2*pi),rejected,'.','MarkerSize',10)
    plot(mod(thetap(inow,:),2*pi),captured,'.','MarkerSize',10)
%   plot(f2,xi2,'LineWidth',2)
    xlim([0,2*pi])
    hold off
    set(gca,'FontSize',18)
    xlabel('\theta')
    ylabel('\gamma')
    legend({'Initial','Rejected','Captured'}, 'Location',"eastoutside")
    figure(2)
    plot(x(1,:)*1e6,xp(1,:)*1e3,'o')
    hold on
    plot(x(inow,:)*1e6,rej_phase*1e3,'.','MarkerSize',10)
    plot(x(inow,:)*1e6,cap_phase*1e3,'.','MarkerSize',10)
    xline(-400e-3,'LineWidth',2)
    xline(400e-3,'LineWidth',2)
    hold off
    xlim([-2,2])
    set(gca,'FontSize',18)
    xlabel('y [\mum]')
    ylabel('y''[mrad]')
    ylim([-1,1])
    
%% Animation of the beam
%     F(nstep)=struct('cdata',[],'colormap',[]);
%     for i=1:nstep
%     figure(4)
%     plot(mod(thetap(i,:),2*pi),gammap(i,:),'.')
%     xlim([0 2*pi])
%     title(['z = ', num2str(i*zstop/nstep)])
%     F(i) = getframe;
%     end
      
end



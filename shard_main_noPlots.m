%% DLA Simulation for SLM Phase Profile SHARD ver 2024
% SHarD Spatial HARmonic DLA simulation
% Assumes single frequency for laser 

% Shard Version 2022. Arbitrary on axis phase and amplitude profile
% are Fourier-analyzed to generate a set of spatial harmonics to describe
% the field in the DLA gap

% Shard Version 2023 uses a better model for the fields in the DLA structures
% which are defined in terms of two complex constant c1 and c2
% This allows to simulate deflecting modes, dipole kicks as well as not symmetric focusing
% channels. 

% Update: this file no longer plots. 
% Update: this file no longer calculates taper, input taper

% Inputs: laser (phase + amplitude), elec (parameters), 
% y0 (initial electron distribution), opt (T if tracking
% all particles, F if throwing away particles as they hit the wall)
% Output: gammap,thetap,x,xp, focusTrack
%% Input Parameters
linearizedfields = true;
npart = elec.n;             % Number of particles

%% Spatial Harmonics Fourier Analysis
E_axis = las.G_gauss.*exp(1i.*las.phi);
dwnsamplerate = max(ceil(nstep/nfreq),1);
E_axis_dn = downsample(E_axis,dwnsamplerate);
dnL0       = length(E_axis_dn);
E_axis_dn = [E_axis_dn, zeros(1,mod(dnL0,2))]; % pad to even
Efft = fftshift(fft(E_axis_dn));
[~,idx] = max(abs(Efft));
nh = size(Efft,2)
kmax = pi/(struc.zstop/nstep*dwnsamplerate);
deltak = kmax/(nh/2);  
kshift = (-nh/2:nh/2-1)*deltak; % zero-centered frequency range

for jh = 1:nh
    kntap(jh) = struc.k + kshift(jh);    
    betantap(jh) = struc.k ./ kntap(:,jh);
    gammantap(jh) = sqrt(1./(1-betantap(:,jh).^2));
end

phi_tap.nh = nh;
phi_tap.kntap = kntap;
phi_tap.betantap = betantap;
phi_tap.Efft = Efft;
phi_tap.linearizedfields = linearizedfields;

%% Tracking
focus = ones([1,npart]);

tic
if opt
    i = 0;
    while (i<nstep)
        i = i+1;
        [ynew, focus_new]=dla_push_RK4(z(i),yold,phi_tap,struc,stepsize,focus);
        filt0 = (1==(repmat(focus_new,1,4)));
        ynew = ynew(filt0~=0);
        yold =ynew;
        focus = ones([1,length(yold)/4]);
    end
else
    focusTrack = ones([nstep,npart]);
    DATA(1,:)=yold;
    for i=1:nstep
        [ynew, focus_new]=dla_push_RK4(z(i),yold,phi_tap,struc,stepsize,focus);
        DATA(i+1,:)=ynew;
        yold=ynew;
        focus = focus_new;
        focusTrack(i,:) = focus;
    end
end
toc

%% Outputs
if 1==opt
    npartN = length(yold)/4;
    gammap= real(ynew(1:npartN));
    thetap= real(ynew(npartN+1:2*npartN));
    x = ynew(2*npartN+1:3*npartN);
    xp = ynew(3*npartN+1:4*npartN);
    focusTrack = focus;
else
    DATA = real(DATA);
    gammap = zeros(nstep,npart);thetap = zeros(nstep,npart);
    x = zeros(nstep,npart);     xp = zeros(nstep,npart);
    for i = 1:nstep
        gammap(i,:)= DATA(i,1:npart);
        thetap(i,:)= DATA(i,npart+1:2*npart);
        x(i,:)= DATA(i,2*npart+1:3*npart);
        xp(i,:)= DATA(i,3*npart+1:4*npart);
    end
end

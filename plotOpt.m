% Plots for save optimization data
% load('taperOpt_G1_3.mat') %(3)
% 
load('taperOpt3.mat') %(1/3)
accelFocMat0 = accelFocMat(1:length(np1s),1:length(np2s));
focusedMat0 = focusedMat(1:length(np1s),1:length(np2s));
np1s0 = np1s; np2s0 = np2s;
% load('taperOpt_G1_2.mat') %(2)
load('taperOpt4.mat') %(2/4)
accelFocMat = [accelFocMat; accelFocMat0];
focusedMat = [focusedMat; focusedMat0];
np1s = [np1s np1s0];

figure(100);
subplot(1,2,1);
surf(np2s,np1s,focusedMat);
view(2); 
xlabel('Slope Taper')
ylabel('Starting Freq')
title('Electrons Focused')
colorbar;
subplot(1,2,2);
surf(np2s,np1s,accelFocMat);
view(2); 
xlabel('Slope Taper')
ylabel('Starting Freq')
title('Electrons Focused and Accelerated')
colorbar;

%%
%  'taperOpt0.mat': 
% blurVar     = 10e-6;        %higher = blurrier

% %electron parameters
% elec.emit = 0.5e-9;
% 
% elec.E0_kin = 5e6;  %5 MeV UCLA Pegasus
% elec.n = 1e3; % Number of particles
% elec.gam0 = (elec.E0_kin+511e3)/511e3;
% elec.beta0 = sqrt(1-1/elec.gam0^2);
% %laser parameters
% las.G0      = 2.0e9; % Gradient
% las.lambda  = 0.8e-6;
% las.k       =  2*pi/las.lambda;
% las.psi_res = params(2);
% 
% nfreq = 120;               % Number of spatial harmonics
% aperture = 0.8e-6;
% 
% Ap = 0;
% res = 0;
% 
% R56, params = [np1 -pi/4 0.15e-3 np2];
%%'taperOpt_G12,3.mat': 
% blurVar     = 10e-6;        %higher = blurrier

% %electron parameters
% elec.emit = 0.5e-9;
% 
% elec.E0_kin = 5e6;  %5 MeV UCLA Pegasus
% elec.n = 1e3; % Number of particles
% elec.gam0 = (elec.E0_kin+511e3)/511e3;
% elec.beta0 = sqrt(1-1/elec.gam0^2);
% %laser parameters
% las.G0      = 1.0e9; % Gradient
% las.lambda  = 0.8e-6;
% las.k       =  2*pi/las.lambda;
% las.psi_res = params(2);
% 
% nfreq = 120;               % Number of spatial harmonics
% aperture = 0.8e-6;
% 
% Ap = 0;
% res = 0;
% 
% R56, params = [np1 -pi/4 0.15e-3 np2];
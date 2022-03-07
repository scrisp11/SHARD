% Input Parameters
opt         = false; %eliminates particles when they hit the walls, saving comp time
R56Tag      = true; %R56 bunch at the start of accel
twissTag    = false; %Use twiss matching parameters
linearizedfields = false;
blur        = true;
genPtTag    = true;

accel       = [0,1]; % [bunchingSection accelSection]
params      = [710 -pi/4 0.15e-3 80]; %[nlambdap psires lengthbunchAccel] (define in loop if looping)
%params(1) can be used as initial periodicity
%params(4) as slope of taper
blurVar     = 10e-6;        %higher = blurrier

%electron parameters
elec.emit = 10e-9;

elec.E0_kin = 5e6;  %5 MeV UCLA Pegasus
elec.n = 1e3; % Number of particles
elec.gam0 = (elec.E0_kin+511e3)/511e3;
elec.beta0 = sqrt(1-1/elec.gam0^2);
%laser parameters
las.G0      = 1.0e9; % Gradient
las.lambda  = 0.8e-6;
las.k       =  2*pi/las.lambda;
las.psi_res = params(2);

nfreq = 120;               % Number of spatial harmonics
aperture = 0.8e-6;

Ap = 0;
res = 0;
%%
% np1s = 500:25:800;
% np1s = 825:25:1000;
np2s = 0:25:200;
% np1s = 712.8682;
% np2s = 126.2239;
% np1s =    1.2038e+03;
np1s = 900:50:1050;
% np1s = 650:10:690;
% np2s = 0:20:200;
np1s = 1000;
np2s = 150;
accelFocMat = zeros(length(np1s),length(np2s));
focusedMat = zeros(length(np1s),length(np2s));
if opt
    if R56Tag % Use the same beam for every sim
        genPtTag    = false;
        yold00 = genParts(las,elec,100,0,1);
        figure;
        scatter(yold00(1001:2000),mod(yold00(1:1000),2*pi));
    end
    plots = false;
    % target = sum(focus); because for now we aren't worrying about exact
    % capture
    figure(44);clc; hold on;
    ylim([0,300]);
    iP = 1;
    for I_np1 = 1:length(np1s)
        np1 = np1s(I_np1)
        for I_np2 = 1:length(np2s)
            np2 = np2s(I_np2)
            params = [np1 -pi/4 0.15e-3 np2];
            yold = yold00;
            SHARD_main2;
            focusedMat(I_np1,I_np2) = focused;
            accelFocMat(I_np1,I_np2) = accelFoc;
            scatter(iP,focused,'rx');
            scatter(iP,accelFoc,'bo');
            iP = iP+1;
            pause(.05);
        end
    end
    hold off;
else
    plots = true;
    SHARD_main2;
end
%%
plotRoutines;
% save('taperOpt_G1_3.mat','focusedMat','accelFocMat','np1s','np2s')
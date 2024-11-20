%% SHarD Quick Plots

%%
E_axis_dn_reconstructed = ifft(ifftshift(Efft));

G_gauss_approx = abs(E_axis_dn_reconstructed);
phi_approx = unwrap(angle(E_axis_dn_reconstructed));
z_reconstructed = downsample(z,dwnsamplerate);

figure(1); clf;
subplot(1,2,1)
plot(kshift,abs(Efft));
title('Spatial Harmonics Content')

subplot(1,2,2)
yyaxis left
plot(z_reconstructed,G_gauss_approx(1:dnL0),'.'); hold on;
plot(z,las.G_gauss);
yyaxis right
plot(z_reconstructed,phi_approx(1:dnL0),'.');
hold on;
plot(z,las.phi)
title('Field Comparison')
L = legend('Reconstructed', 'Input');
L.Location = 'Southeast';
%%
if plotParticleDists
    figure(3); clf;
    plot(x(1,focus)*1e6,xp(1,focus)*1e3,'bo');hold on;
    plot(x(1,~focus)*1e6,xp(1,~focus)*1e3,'ro')
    if ~opt
        plot(x(end,:)*1e6,cap_phase*1e3,'b.','MarkerSize',10)
    end
    xline(-struc.gap/2*1e6,'LineWidth',2);    xline(struc.gap/2*1e6,'LineWidth',2);
    xlim([-2*struc.gap/2*1e6,2*struc.gap/2*1e6])
    legend('z = 0 (Transmitted)','z = 0 (Walls)','z = end (Transmitted)')
    
    xlabel('y [\mum]')
    ylabel('y''[mrad]')
    ylim([-1,1])
    
    if ~opt
        iPs = [1,nstep];
    else
        iPs = 1;
    end
    figure(4);clf;
    for ii = 1:length(iPs)
        iP = iPs(ii);
        subplot(1,length(iPs),ii)
        scatter((x(iP,~focus))*1e6,mod(thetap(iP,~focus)+taper(iP) - pi,2*pi) - pi,'ro'); hold on;
        plot((x(iP,focus))*1e6,mod(thetap(iP,focus)+taper(iP) - pi,2*pi) - pi,'go'); 
        plot((x(iP,capt))*1e6,mod(thetap(iP,capt)+taper(iP) - pi,2*pi) - pi,'bo'); hold off;
        xline(-struc.gap*1e6*.5,'LineWidth',2);
        xline(struc.gap*1e6*.5,'LineWidth',2);
        xlim([-2*struc.gap/2*1e6,2*struc.gap/2*1e6])
        ylim([-pi,pi])
        xlabel('y [\mum]')
        ylabel('phi [rad]')
        title(['z = ',num2str(z(iP)*1e3),' mm'])
    end
end
%%
binLims = .511*[0.95*min(gammap(end,focus)-elec.gam0),1.05*max(gammap(end,focus)-elec.gam0)];

if binLims(1)>binLims(2); binLims = flip(binLims);end

f = figure(5); clf; f.Units  = 'centimeters';
set(gcf,'renderer','Painters')
f.Position = [1,1,12,8];
histogram(.511*gammap(end,focus)-(elec.gam0*.511),'NumBins',50,'BinLimits',binLims,'FaceColor','r'); 
ylabel('Count')

%%
if ~opt
    figure(6);clf;
    if sum(focus)>1
        plot(z*1e3, .511*max(gammap(:,focus)-elec.gam0,[],2)); hold on;
    else
        plot(z*1e3, .511*max(gammap(:,:)-elec.gam0,[],2)); hold on;
    end
    plot(z*1e3,real((gam_res-elec.gam0)*0.511),'--'); hold off;
    legend('Simulation','Model')
    title('Maximum Energy Particle')
    ylabel('Energy (MeV)')
    xlabel('z (mm)')
end    

%%
[~,iP] = maxk(gammap(end,focus),10);
gpF = gammap(:,focus);xF = x(:,focus);

%% Change in phase over course of structure (Longitudinal Stability)
if ~opt
    fCapt = find(capt);
    xPl = 1e3*z'; xPlx = [xPl,xPl];
    thF = thetap(:,focus);

    figure(7); clf; 
    set(gcf,'renderer','Painters')

for iPlot = 1:min([10,sum(capt)])
    yPl = (thF(:,iP(iPlot)) + taper'); yPly = [yPl, yPl];
    zPlz = zeros(size(xPlx));
    cols = .511*(gpF(end,iP(iPlot)) - gpF(1,iP(iPlot)))*ones(size(yPly));
    hs = surf(xPlx,yPly,zPlz,cols,'EdgeColor','Interp'); hold on;
    hs.LineWidth = 1;
end
view(2)
a = colorbar;
ylabel(a,'\Delta E (MeV)','fontsize',16,'Rotation',90);
a.Label.Position(1) = 3.5;
xlabel('z (mm)')
ylabel('\phi (rad)')
grid off
box on
end

%% Energy, Transverse Position along structure
if ~opt
f = figure(8);clf;
set(gcf,'renderer','Painters')
f.Units  = 'centimeters';
f.Position = [5,5,14,10];
% hannah = [
%      21, 173, 173;
%      62, 138,  61;
%      172,  13,  27;
%     ]./256;
% nColors = 256*8;
% cmapSurf = zeros(nColors,3);
% for i = 1:3
%     cmapSurf(1:nColors/2, i) = linspace(hannah(1,i),hannah(2,i), nColors/2);
% end
% for i = 1:3
%     cmapSurf(nColors/2+1:nColors, i) = linspace(hannah(2,i),hannah(3,i), nColors/2);
% end
rainbowColormap;
for ii = 1:1:min([length(iP),50])
xx=[z*1e3; z*1e3];           %// create a 2D matrix based on "X" column
yy=[1e9*xF(:,iP(ii)), 1e9*xF(:,iP(ii))]';           %// same for Y
zz=zeros(size(xx)); %// everything in the Z=0 plane
cc =[.511*(gpF(:,iP(ii))-elec.gam0), .511*(gpF(:,iP(ii))-elec.gam0)]';         %// matrix for "CData"
hold on;
hs=surf(xx,yy,zz,cc,'EdgeColor','interp','FaceColor','none') ;
hs.LineWidth = 1.5;
view(2)
colormap(flip(rainbowCMap))
end
grid off
ylabel('y (nm)')
xlabel('z (mm)')
c = colorbar;
yl = ylabel(c,'\Delta{E} (MeV)','Rotation',90);
%     saveas(gcf,['C:\Users\sophi\Documents\UCLA Grad\ACHIP\2023_2024_Run\Paper\MatlabFigures\shard_Pond_traj_bunched'],'epsc')
%     saveas(gcf,['C:\Users\sophi\Documents\UCLA Grad\ACHIP\2023_2024_Run\Paper\MatlabFigures\shard_Pond_traj_bunched'],'pdf')
%         saveas(gcf,['C:\Users\sophi\Documents\UCLA Grad\ACHIP\2023_2024_Run\Paper\MatlabFigures\shard_Pond_traj_bunched'],'png')
end
%% Change in position over course of structure (Transverse Stability)
% if ~opt
% xPl = 1e3*z'; xPlx = [xPl,xPl];
% 
% figure(9); clf; 
% set(gcf,'renderer','Painters')
% 
% for iPlot = 1:1:sum(capt)
%     iCapt = fCapt(iPlot);
%     yPl = 1e6*x(:,iCapt); yPly = [yPl, yPl];
%     zPlz = zeros(size(xPlx));
%     cols = .511*(gammap(end,iCapt) - gammap(1,iCapt))*ones(size(yPly));
%     hs = surf(xPlx,yPly,zPlz,cols,'EdgeColor','Interp'); hold on;
%     hs.LineWidth = 1;
% end    
% 
% view(2)
% a = colorbar;
% ylabel(a,'\Delta E (MeV)','fontsize',16,'Rotation',90);
% a.Label.Position(1) = 3.5;
% xlabel('z (mm)')
% ylabel('y (\mu{m})')
% grid off
% box on
% end
%%
f = figure(10); clf; f.Units  = 'centimeters';
f.Position = [10,10,30,6];

iMM(1) =  find([1:nstep]*struc.zstop/nstep>=0e-3,1);
iMM(2) = find([1:nstep]*struc.zstop/nstep>=st,1);
% iMM(3) = find([1:nstep]*struc.zstop/nstep>=2.5e-3,1);
iMM(3) = find([1:nstep]*struc.zstop/nstep>=st + dr,1);

if iMM(2)==iMM(3)&dr~=0;    iMM(3) = iMM(3) + 1; end
iMM(4) = nstep;
if ~opt; iis = 1:4; else iis = 1; end
for i = iis
    fi = focusTrack(iMM(i),:)==1;
    subplot(1,max(iis),i)
    hold on
    plot(mod(thetap(iMM(i),fi) + taper(iMM(i)) - pi,2*pi) - pi,.511*(gammap(iMM(i),fi)-elec.gam0),'b.','MarkerSize',5);
    plot(mod(thetap(iMM(i),focus)+ taper(iMM(i)) - pi,2*pi) - pi,.511*(gammap(iMM(i),focus)-elec.gam0),'r.','MarkerSize',5);
    plot(mod(thetap(iMM(i),capt)+ taper(iMM(i)) - pi,2*pi) - pi,.511*(gammap(iMM(i),capt)-elec.gam0),'m.','MarkerSize',5);
    xlim([-pi,pi])
    hold off
    box on
    xlabel('\theta')
    title(['z = ',num2str(floor(1e6*iMM(i)*struc.zstop/nstep)/1e3,2),' mm'])
end

% subplot(1,max(iis),1)
ylabel('\Delta E (MeV)')
L = legend('Current','Transmitted','Accelerated');
L.ItemTokenSize(1) = 3;

%% Ponderomotive Stability
psi2        = pi*sign(psi_res) - psi_res;
km_eq       = @(m,kp) (struc.k + (m +1).*kp);
alpham_eq   = @(m,Ap,zi,kp) (struc.sf*las.G_gauss(zi).*besselj(m,Ap)./(km_eq(m,kp)*511e3));
eta_eq      = @(m,psi,gam_res,Ap,zi,kp,psi_res) (sqrt(2.*alpham_eq(m,Ap,zi,kp).*gam_res.*(cos(psi) - cos(psi2(zi)) + (psi - psi2(zi)).*sin(psi_res))));

psis        = linspace(0,2*pi,100);
%%
if animation
%     F(nstep)=struct('cdata',[],'colormap',[]);
        figure(40)
    for i=1:10:nstep
    clf;
    plot(mod(thetap(i,capt)+taper(i)-pi,2*pi) - pi,.511*gammap(i,capt),'.'); hold on;
    plot(mod(thetap(i,~capt)+taper(i)-pi,2*pi) - pi,.511*gammap(i,~capt),'.');
    xlim([-pi pi])
    if pondTag
        etas        = eta_eq(res(i),psis,gam_res(i),Ap(i),i,2*pi/lambdap(i),psi_res(i));
        gamsp        = .511*(etas+1)*gam_res(i);
        gamsm        = .511*(-etas+1)*gam_res(i);
        scatter(mod(psis,2*pi) - pi,gamsp,'.'); hold on;
        scatter(mod(psis,2*pi) - pi,gamsm,'.'); hold off;
    end
    ylim([min(.511*gammap(end,:)),1.2*max(.511*gammap(end,:))])
    title(['z = ', num2str(i*struc.zstop/nstep)])
%     F(i) = getframe;
    pause;
    end
        ylabel('Energy (MeV)')
    xlabel('Phase (Rad)')
%%
%     F(nstep)=struct('cdata',[],'colormap',[]);
    figure(40);clf;
% 
    for i=1:5:nstep
    subplot(1,2,1)
    plot((x(i,1==focusTrack(i,:)))*1e6,xp(i,1==focusTrack(i,:))*1e3,'.','MarkerSize',10);hold on;
    plot((x(i,0==focusTrack(i,:)))*1e6,xp(i,0==focusTrack(i,:))*1e3,'.','MarkerSize',10);hold off;
    title(['z = ', num2str(i*struc.zstop/nstep)])
    ylim([-1,1])
    xlim([-1e6*struc.gap/2*1.1,1e6*struc.gap/2*1.1])
    xlabel('y [\mum]')
    ylabel('y''[mrad]')
    xline(-struc.gap*1e6*.5,'LineWidth',2);
    xline(struc.gap*1e6*.5,'LineWidth',2);
    ylim([-1,1])
    subplot(1,2,2)
    subplot(1,2,2)
    plot(z*1e3,abs(las.phi-taper)); hold on;
    scatter(z(i)*1e3,abs(las.phi(i)-taper(i))); hold off;
    xlabel('z (mm)')
    ylabel('Apf Phi (rad)')
    pause(0.1)
%     F(i) = getframe;
    end
%%     
figure(45);clf;
% 
    for i=1:5:nstep
    subplot(1,2,1)
    plot((x(i,1==focusTrack(i,:)))*1e6,gammap(i,1==focusTrack(i,:))*.511,'.','MarkerSize',10);hold on;
    plot((x(i,0==focusTrack(i,:)))*1e6,gammap(i,0==focusTrack(i,:))*.511,'.','MarkerSize',10);hold off;
    title(['z = ', num2str(i*zstop/nstep)])
    ylim([(gam_res(i)-0.4)*.511,(gam_res(i)+0.4)*.511])
    xlim([-1,1])
    xlabel('y (\mum)')
    ylabel('Energy (MeV)')
    xline(-struc.gap*1e6*.5,'LineWidth',2);
    xline(struc.gap*1e6*.5,'LineWidth',2);
    subplot(1,2,2)
    subplot(1,2,2)
    plot(z*1e3,abs(las.phi-taper)); hold on;
    scatter(z(i)*1e3,abs(las.phi(i)-taper(i))); hold off;
    xlabel('z (mm)')
    ylabel('Apf Phi (rad)')
%         F(i) = getframe;
    end
end

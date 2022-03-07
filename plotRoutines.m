% All Plots
vids = false;
% plots = true;
%Initial distribution in both phase spaces
% if plots
%     figure(40)
%     plot(theta0,gamma0,'.')
%     title('Initial Particle Distribution')
%     figure(41)
%     plot(x0,x0p,'.')
%     title('Initial Particle Distribution')
% end

% Phase and amplitude of mask on same plot
if plots
    figure(33)
    yyaxis left
    plot(z, abs(G_gauss.*1e-9),'LineWidth',1.5)
    ylabel('E, GV/m')
    yyaxis right
    plot(z,phi)
    ylabel('Phase, rad')
    xlabel('z, m')
    title('Electric Field Strength')
end

%Energy particle over structure
if plots
    figure(31)
    plot(z,gam_res*0.511)
    ylabel('Energy (MeV)');
    xlabel('z [m]')
end

% Electric field (Actual)
if plots
%     figure(32)
%     plot(z,mod(phi,2*pi))
%     xlabel('z [m]')
%     title('Electric Field Phase Profile')
end

%Spatial Harmonics
if plots
    figure(34)
    plot(kshift,abs(Efft))
    title('Spatial Harmonics Content')
    drawnow
end


if plots
    clear captured rejected cap_phase rej_phase

    disp(['Resonant Energy Gain (in keV) = ' num2str((gam_res(end)-gam0)*511)])
    disp(['Accelerated Particles = ' num2str(sum(capt))])
    disp(['Focused Particles = ' num2str(sum(focus))])
    disp(['Total captured particles =' num2str(sum(capt.*focus))])
    disp(['Target cost function =' num2str(target)])
    
%    [f2,xi2] = ksdensity(gammap(inow,:));
    
    %f2=5*f2; %5 is just a plotting factor
    
    captured = gammap(inow,:);  % All energies
    rejected = gammap(inow,:); 
    focused = gammap(inow,:);
    captured(~focus) = NaN;     % filter out walls
    captured(~capt) = NaN;
    rejected(focus) = NaN;
    focused(~focus) = NaN;
    
    cap_phase = xp(inow,:);
    rej_phase = xp(inow,:);
    cap_phase(~focus) = NaN;
    rej_phase(focus) = NaN;
%     
%     figure(20);
%     scatter(x(inow,:),xp(inow,:),15,focus);
%     xlim([-2e-6, 2e-6])
%     ylim([-1e-3, 1e-3])
%     
%     figure(21);
%     subplot(2,1,1);
%     scatter(mod(thetap(inow,:),2*pi), gammap(inow,:),15,focus);
%     subplot(2,1,2);
%     scatter(mod(thetap(inow,:),2*pi), gammap(inow,:),15,capt);
%     
    
    figure(1);
    plot(mod(thetap(1,:),2*pi),gammap(1,:),'.','MarkerSize',10);
    hold on
    plot(mod(thetap(inow,:),2*pi),rejected,'r.','MarkerSize',10);
    plot(mod(thetap(inow,:),2*pi),focused,'b.','MarkerSize',10);
    plot(mod(thetap(inow,:),2*pi),captured,'g.','MarkerSize',10);
%   plot(f2,xi2,'LineWidth',2)
    xlim([0,2*pi])
    hold off
    set(gca,'FontSize',18)
    xlabel('\theta')
    ylabel('\gamma')
    title(['z = ', num2str(10^3*z(inow)),'mm'])
    legend({'Initial','Rejected','Focused','Captured'}, 'Location',"eastoutside")
%     figure(2)
%     plot(x(1,:)*1e6,xp(1,:)*1e3,'o');
%     hold on
%     plot(x(inow,:)*1e6,rej_phase*1e3,'.','MarkerSize',10)
%     plot(x(inow,:)*1e6,cap_phase*1e3,'.','MarkerSize',10);
%     xline(-400e-3,'LineWidth',2);
%     xline(400e-3,'LineWidth',2);
%     hold off
%     xlim([-2,2])
%     set(gca,'FontSize',18)
%     xlabel('y [\mum]')
%     ylabel('y''[mrad]')
%     ylim([-1,1])
    figure(6);
    nb =-2;
    plot(mod(thetap( sum(z<bunchLength)+nb,~focus),2*pi),gammap(sum(z<bunchLength)+nb,~focus),'r.','MarkerSize',10);
    hold on;
    plot(mod(thetap( sum(z<bunchLength)+nb,focus),2*pi),gammap(sum(z<bunchLength)+nb,focus),'b.','MarkerSize',10);
    plot(mod(thetap( sum(z<bunchLength)+nb,1==focus.*capt),2*pi),gammap(sum(z<bunchLength)+nb,1==focus.*capt),'g.','MarkerSize',10); hold off;
        title(['z = ', num2str(10^3*z(sum(z<bunchLength))),'mm'])
    xlabel('\theta')
    ylabel('\gamma')
    legend('Rejected','Focused','Focused and Accel')
    figure(7);
    plot(mod(thetap( 1,~focus),2*pi),gammap(1,~focus),'r.','MarkerSize',10);
    hold on;
    plot(mod(thetap( 1,focus),2*pi),gammap(1,focus),'b.','MarkerSize',10);
    plot(mod(thetap( 1,1==focus.*capt),2*pi),gammap(1,1==focus.*capt),'g.','MarkerSize',10); hold off;
        title(['z = ', num2str(10^3*z(1)),'mm'])
    xlabel('\theta')
    ylabel('\gamma')
end

if vids
    F(nstep)=struct('cdata',[],'colormap',[]);
    for i=1:nstep
    figure(4)
    plot(mod(thetap(i,:),2*pi),gammap(i,:),'.')
    xlim([0 2*pi])
    title(['z = ', num2str(i*zstop/nstep)])
    F(i) = getframe;
    pause(0.1)
    end
end
% clear all
close all

%% conditions for transverse stability 
lambda= 780e-9; k = 2*pi/lambda;

E1 = 2e9;
alpha = E1/(511e3*k);

j = 1;
%% APF beta
psir = pi/4;
for betagamma_exp = 1.0:0.03:2.4
    betagamma(j) = 0.1*10^betagamma_exp;
    i = 1;
    for beatperiod= 20:20:3000
        beatper(i) = beatperiod;
        gamma = sqrt(1+betagamma(j).^2);
        beta = sqrt(1 - 1./gamma.^2);
        K_0 = sqrt(alpha*k^2./ gamma.^3 ./ beta.^2 * cos(psir));
        fillfactor = 1;
        L = beatperiod*fillfactor*lambda;
        Ld = (1-fillfactor)*beatperiod*lambda/2;
        O = [1 Ld; 0 1];
        F = [cos(K_0.*L/4) 1./K_0.*sin(K_0.*L/4); -K_0.*sin(K_0.*L/4) cos(K_0.*L/4)];
        D = [cosh(K_0.*L/4) 1./K_0*sinh(K_0.*L/4); K_0.*sinh(K_0.*L/4) cosh(K_0.*L/4)];
        R = F*O*D*D*O*F;
        R2 = D*O*F*F*O*D;
        if(abs(trace(R))<2 & abs(trace(R2))<2)
            beta_hat(i,j) = abs(real(R(1,2))/sqrt(1-R(1,1)^2));
        else
            beta_hat(i,j) = NaN;
        end
        i = i+1;
    end
    j = j+1;
end

%% Plot beta max
figure(1);clf;

[bg,bp] = meshgrid(betagamma,beatper);

gap = 400e-9;
norm_emit = 1e-9;
sigma_max = sqrt(beta_hat./bg*norm_emit);

mask = ones(size(sigma_max));indV = [];
for ibp = 1:length(beatper)
    ind0 = find(isnan(sigma_max(ibp,:)),1,'last');
    if isempty(ind0); indV(ibp) = 1;else;indV(ibp) = ind0;end
    mask(ibp,1:indV(ibp)) = NaN;
end

[envelope,iMinSigmamax] = min(sigma_max.*mask);
contour(bg,bp,log10(sigma_max/gap),30)
zindex = 0;
hold on
contour(bg,bp,log10(sigma_max/gap),[0 0],'LineWidth',2);
plot(betagamma,beatper(iMinSigmamax));
plot(betagamma(indV),beatper,'m');

hold off
legend('log10(sigma/400nm)')
xlabel('relativistic gamma')
ylabel('modulation period')
title('Norm emittance 1 nm')
colorbar
%xlim([0.1 250])

ind = find(sigma_max./gap<1);
% bg(ind(end))
LUT = beatper(iMinSigmamax);

figure(2);clf;
plot(bg,envelope*1e9,'k');
xlabel('\beta\gamma')
ylabel('Beam Envelope (nm)')
%%
% save(['LUT_APFstable_psiRes_',num2str(psir*1e4,'%3.0f')],'betagamma','E1','beatper','lambda','LUT')
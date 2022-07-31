clear all
close all

%% conditions for transverse stability 
lambda= 800e-9;
E1 = 2e9;
alpha = 0.311*E1/1e6*lambda;
k = 2*pi/lambda;
j = 1;
%% APF beta
psir = pi/4;
for betagamma_exp = 1.0:0.03:4
    betagamma(j) = 0.1*10^betagamma_exp;
    i = 1;
    for beatperiod= 20:20:200000
        beatper(i) = beatperiod;
        gamma = sqrt(1+betagamma(j).^2);
        beta = sqrt(1 - 1./gamma.^2);
        K_0 = sqrt(alpha*k^2./ gamma.^3 ./ beta.^2 * cos(psir));
        fillfactor = 0.9;
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
[bg,bp] = meshgrid(betagamma,beatper);

gap = 1e-6;
norm_emit = 1e-9;
sigma_max = sqrt(beta_hat./bg*norm_emit);
contour(bg,bp,log10(sigma_max/gap),30)
zindex = 0;
hold on
contour(bg,bp,log10(sigma_max/gap),[0 0],'LineWidth',2)
hold off
legend('log10(sigma/1um)')
xlabel('relativistic gamma')
ylabel('modulation period')
title('Norm emittance 1 nm')
colorbar
%xlim([0.1 250])

ind = find(sigma_max./gap<1);
bg(ind(end))

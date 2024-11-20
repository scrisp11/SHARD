function dy = dla_push(z,y,phi_tap,struc,stepsize)
npart=size(y,2)/4;
k = struc.k;

%% Reshaping
gamma=y(1:npart);
tp=y(npart+1:npart*2);
x = y(2*npart+1:npart*3);
xp = y(3*npart+1:npart*4);
beta = sqrt(gamma.^2-1)./gamma;

%% Equations for particle motion
dgamma = zeros([1,npart]);
dxp = zeros([1,npart]);
k_perpM = sqrt(phi_tap.kntap.^2-(k)^2);
for jh = 1:phi_tap.nh
    k_perp = k_perpM(jh);
%     k_perp = struc.Gamma;
    Efield = phi_tap.Efft(jh)/phi_tap.nh;
    
    transv  = (struc.c1*exp(-k_perp.*x)+struc.c2*exp(k_perp.*x));
    transvd =   -struc.c1*exp(-k_perp.*x)+struc.c2*exp(k_perp.*x);

    if phi_tap.linearizedfields
        dgamma = dgamma - abs(struc.c1+struc.c2).*1/511e3*imag(Efield.*exp(1i*(tp+(phi_tap.kntap(jh)-k).*z)));
        dxp = dxp +  abs(struc.c1+struc.c2).*1/511e3.*x./gamma.*(1-phi_tap.betantap(jh).*beta)*phi_tap.kntap(jh).*...
            real(Efield.*exp(1i*(tp+(phi_tap.kntap(jh)-k).*z)));
    else
        dgamma = dgamma - 1/511e3*imag(Efield.*transv.*exp(1i*(tp+(phi_tap.kntap(jh)-k).*z)));
        if(k_perp ~= 0)
            dxp = dxp + 1/511e3.*(1-phi_tap.betantap(jh).*beta)./gamma*...
                phi_tap.kntap(jh).*...
                real(Efield.*transvd./k_perp.*exp(1i*(tp+(phi_tap.kntap(jh)-k).*z)));
        else
            dxp = dxp + 1/511e3.*x./gamma.*(1-phi_tap.betantap(jh).*beta)*...
                phi_tap.kntap(jh).*...
                real(Efield.*transv.*exp(1i*(tp+(phi_tap.kntap(jh)-k).*z)));
        end

        dgamma(abs(x*k_perp)>struc.gap*max(real(k_perpM))) = 0;
        dxp(abs(x*k_perp)>struc.gap*max(real(k_perpM))) = 0;
    end
    dtp = k*(1-1./beta);
    dx = xp;

%% Function output
dy = [dgamma dtp dx dxp];
end
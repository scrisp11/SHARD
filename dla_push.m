function dy = dla_push(z,y,phi_tap,k,stepsize)
npart=size(y,2)/4;

%% Reshaping
gamma=y(1:npart);
tp=y(npart+1:npart*2);
x = y(2*npart+1:npart*3);
xp = y(3*npart+1:npart*4);


%% Equations for particle motion
dgamma = zeros([1,npart]);
dxp = zeros([1,npart]);
for jh = 1:phi_tap.nh
    k_perp = sqrt(phi_tap.kntap(jh)^2-k^2);
    
    if phi_tap.linearizedfields
        dgamma = dgamma-1/511e3*abs(phi_tap.Efft(jh)/phi_tap.nh).*sin(tp+(phi_tap.kntap(jh)-k).*z+angle(phi_tap.Efft(jh))); %.*cosh(k_perp.*x);
        dxp = dxp + 1/511e3.*x./gamma.*(1-phi_tap.betantap(jh).*sqrt(gamma.^2-1)./gamma)*phi_tap.kntap(jh).*...
            abs(phi_tap.Efft(jh)/phi_tap.nh).*cos(tp+(phi_tap.kntap(jh)-k).*z+angle(phi_tap.Efft(jh)));
        
    else 
        dgamma = dgamma-1/511e3*abs(phi_tap.Efft(jh)/phi_tap.nh).*sin(tp+(phi_tap.kntap(jh)-k).*z+angle(phi_tap.Efft(jh))).*cosh(k_perp.*x);
        if(k_perp ~= 0)
            dxp = dxp + 1/511e3.*abs(sinh(k_perp.*x)./k_perp)./gamma.*(1-phi_tap.betantap(jh).*sqrt(gamma.^2-1)./gamma)*phi_tap.kntap(jh).*...
                abs(phi_tap.Efft(jh)/phi_tap.nh).*cos(tp+(phi_tap.kntap(jh)-k).*z+angle(phi_tap.Efft(jh)));
        else
            dxp = dxp + 1/511e3.*x./gamma.*(1-phi_tap.betantap(jh).*sqrt(gamma.^2-1)./gamma)*phi_tap.kntap(jh).*...
                abs(phi_tap.Efft(jh)/phi_tap.nh).*cos(tp+(phi_tap.kntap(jh)-k).*z+angle(phi_tap.Efft(jh)));
        end
        dgamma(abs(k_perp*x)>1) = 0;
        dxp(abs(k_perp*x)>1) = 0;
    end
dtp = k*(1-gamma./sqrt(gamma.^2-1));
dx = xp;

%% Function output
dy = [dgamma dtp dx dxp];
end









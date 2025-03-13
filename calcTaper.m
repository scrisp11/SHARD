function [taper, las, gam_res,tap] = calcTaper(elec,struc,las,z,taperOn)
%    
    nstep = length(z); stepsize        = struc.zstop/(nstep-1);
    gam_res=zeros(1,nstep);         beta_res=zeros(1,nstep);    
    phi_add = zeros(1,nstep);       tap=zeros(1,nstep);
    Ezfun           = @(c,x)((c(1)*exp(-struc.Gamma*x)+c(2)*exp(struc.Gamma*x)));

    gam_res(1)=elec.gam0;
    beta_res(1)=sqrt(1-1/gam_res(1)^2);
    tap(1)=1/beta_res(1)-1+sin(las.thI)*las.k/struc.k;

    phi_add(1) = struc.k.*tap(1)*(z(2)-z(1))/2;
    % phi_add(1) = 0;
    for i=2:nstep
        gam_res(i)=gam_res(i-1)-las.G_gauss(i).*abs(Ezfun([struc.c1,struc.c2],struc.yc))/511e3*stepsize*sin(las.phi(i));
        beta_res(i)=sqrt(1-1./gam_res(i).^2);
        tap(i)=1/beta_res(max([i*taperOn,1]))-1+sin(las.thI)*las.k/struc.k;
        phi_add(i) = trapz(z(2)-z(1),struc.k.*tap);
    end

    gradient_GeVperm = (gam_res(end)-gam_res(1))*0.511/struc.zstop*0.001;
    maxenergy_MeV = gam_res(end)*0.511;
    table(gradient_GeVperm,maxenergy_MeV)
    las.phi     = las.phi + phi_add;
    taper       = phi_add;
end
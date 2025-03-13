function [y,frac] = genParts(elec,psi_res, bunching,struc)
    sigmaxp = elec.emit/elec.gam0/elec.beta0/elec.sigmax;
    seed = 1;
    rng(135);
   
    nK      = 0;
    yK      = [];
    nTot    = 0;
    while nK < elec.n
        theta0 = hammersley(1,elec.n,seed)*2*pi-pi;
        gamma0 = normrnd(elec.gam0,elec.deltagamma,[1,elec.n]);
        % add bunching
        if bunching
             Ab = 5*elec.deltagamma;
             R56 = pi/Ab/2;
             gamma0 = gamma0 - Ab.*sin(theta0);
             theta0 = theta0+(gamma0-mean(gamma0)).*R56-psi_res(1)-pi;
        else
            theta0 = hammersley(1,elec.n,seed)*2*pi-pi;
        end

        x0 = normrnd(0,elec.sigmax,[1,elec.n]);
        x0p = normrnd(0,sigmaxp,[1,elec.n]);        

        filt = (abs(x0)<struc.gap/2);
        gamma0 = gamma0(filt);
        x0 = x0(filt);
        x0p = x0p(filt);
        hold on;    
        yK0 = [gamma0; x0; x0p];
        yK  = [yK yK0];
        nK = length(yK);
        nTot = nTot + elec.n;
    end
    frac = nK/nTot;
    y = yK(:,1:elec.n);
    y = [y(1,:) theta0(1:elec.n) y(2,:) y(3,:)];
end



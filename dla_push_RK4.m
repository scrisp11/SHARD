function [ynew,focus_new] = dla_push_RK4(zold,yold,phi_tap,struc,stepsize,focus)
focus_new = aperture(focus,yold,struc);

k1=stepsize*dla_push(zold,yold,phi_tap,struc,stepsize);
k2=stepsize*dla_push(zold+stepsize/2,yold+k1/2,phi_tap,struc,stepsize);
k3=stepsize*dla_push(zold+stepsize/2,yold+k2/2,phi_tap,struc,stepsize);
k4=stepsize*dla_push(zold+stepsize,yold+k3,phi_tap,struc,stepsize);

ynew=yold+1/6*(k1+2*k2+2*k3+k4);
end


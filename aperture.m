function focus_new = aperture(focus_old,y,struc)
npart=size(y,2)/4;
x = y(2*npart+1:npart*3);

%% Aperture
gapsize = struc.gap;
focus_new = focus_old.*(abs(x(1,:))<gapsize/2);

end
load ts.dat;
load nf.dat;
load qf.dat
load Tef.dat
load Tif.dat
load Vthf.dat
load Vzf.dat
load Erf.dat

ind = find (nf(:,1)==ts(:,1))

qf(ind,2)-ts(:,2)



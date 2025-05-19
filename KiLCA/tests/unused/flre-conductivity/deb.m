dta = load('fort.100');

figure;
plot(dta(:,1), dta(:,2),'r.', dta(:,1), dta(:,3),'b.');

sz = size(dta);

len = (sz(2) - 3) / 2; % number of complex eqs

figure;
plot(dta(:,1), dta(:,4:3+len), 'r', dta(:,1), dta(:,4+len:4+2*len-1), 'b');

%  dta2 = load('fort.200');
%
%  xx = dta2(1:100:10000,1);
%  yy = dta2(1:100,2);
%
%  ref = reshape(dta2(:,3), [100, 100]);
%  imf = reshape(dta2(:,4), [100, 100]);
%
%  %fun = log(ref + i*imf);
%  fun = ref + i*imf;
%
%  figure;
%  subplot(1,2,1);
%  mesh(xx, yy, real(fun));
%  %plot3(0,0,0,'*');
%
%  subplot(1,2,2);
%  mesh(xx, yy, imag(fun));
%  %plot3(0,0,0,'*');

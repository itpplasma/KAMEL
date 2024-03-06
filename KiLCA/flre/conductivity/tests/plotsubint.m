function w = plotsubint()

tmax = 3.0e2

m = 0;
n = 0
l = 20;

x1 = -0.15548575383155275
x2 = -1.3981937829735058
x3 = -2.0838042516747982E-002
x4 = -2.6022827172207063E-002

  x1 = -189.02963851415251
  x2 = -48.812002536152967
  x3 = -0.54539986769432813
  x4 =  0.56254880990715950

  x1 =  -54.373205158510601
  x2 = -913179.34786194994
  x3 = 2.4474542564594847E-003
  x4 = -9.7688538680105260E-004

 x1 = -189.02963851415251
 x2 = -48.812002536152967
 x3 = -0.54539986769432813
 x4 =  0.56254880990715950

%  x1 = -0.15548575383155275
%  x2 = -1.3981937829735058
%  x3 = -2.0838042516747982E-002
%  x4 = -2.6022827172207063E-002

alpha = sqrt(0.25e0 - x3*i) - 0.5e0;

gamma = alpha / (1.0 + 2.0e0*alpha);

x1 = x1 / (1.0e0 + 2.0e0*alpha)^1.5;

x2 = x2 / (1.0e0 + 2.0e0*alpha) + gamma * i;

x4 = x4 / (1.0e0 + 2.0e0*alpha);

x = 0.0;

y = 0.0;

tau = linspace(0.0, tmax, 10000);

arg = (-x1*x1 + i*x2)

taumax = -40.0 / real(arg)

phi = - atan(imag(arg)/real(arg))

si = subint(tau, l, x, y, x1, x2, x3, x4);
sil = subintlim(tau, l, x, y, x1, x2, x3, x4);
sia = subintarg(tau, l, x, y, x1, x2, x3, x4);

figure('Name', 'si');

subplot(1,2,1)
plot(tau, real(si), 'r', tau, imag(si), 'b', taumax, 0, 'r*');

subplot(1,2,2)
plot(real(si), imag(si));

figure('Name', 'si arg');

plot(tau, real(sia), 'r', tau, imag(sia), 'b');

[q, errbnd] = quadgk(@(z)subint(z, l, x, y, x1, x2, x3, x4), 0.0, inf, 'RelTol', 1e-12, 'AbsTol', 1e-12);

q * 1.0e0 /(1.0e0 + 2.0e0*alpha)^((3 + m + n)/2)

errbnd

kmax = 100;
tau_z = linspace(-10,10,kmax);
siarg = zeros(kmax,kmax);

for k1=1:kmax
    for k2=1:kmax
    
        siarg(k1,k2) = subintarg(tau_z(k1)+tau_z(k2)*i, l, x, y, x1, x2, x3, x4);

    end
end

figure;
subplot(1,2,1); hold on;
mesh(tau_z, tau_z, real(siarg'));
plot3(0,0,0,'*');

subplot(1,2,2);
mesh(tau_z, tau_z, imag(siarg'));

phimin = - pi/2;
phimax =   pi/2;
phimin = - pi/4;
phimax =   pi/4;

phiopt =   phi;

if(real(arg)*cos(phimin)-imag(arg)*sin(phimin) < real(arg)*cos(phimax)-imag(arg)*sin(phimax))
    phi = phimin;
else
    phi = phimax;
end

if(real(arg)*cos(phiopt)-imag(arg)*sin(phiopt) < real(arg)*cos(phi)-imag(arg)*sin(phi))
    phi = phiopt;
end

phi = -1.570796326794896

phi

% contour:
ll = linspace(0,tmax,10000);

zz = ll*cos(phi)+ i*ll*sin(phi);

sia_c = subintarg(zz, l, x, y, x1, x2, x3, x4);

tt = linspace(0,-1d9,1000);

cut = -0.5 * log(((gamma-1)^2 - tt)/gamma^2);

llmax = - 40.0 / real(arg*(cos(phi)+ i*sin(phi)))

zzmax = llmax*(cos(phi)+ i*sin(phi));

figure('Name', 'contour');

subplot(1,2,1);
plot(real(zz), imag(zz), 'b', real(cut), imag(cut), 'r', real(zzmax), imag(zzmax), 'ro');

subplot(1,2,2); hold on;
plot(ll, real(sia_c), 'r.', ll, imag(sia_c), 'b.', tau, real(sia), 'r', tau, imag(sia), 'b', llmax, 0, 'ro');

zzmin = ll*cos(phimin)+ i*ll*sin(phimin);
siamin = subintarg(zzmin, l, x, y, x1, x2, x3, x4);

zzmax = ll*cos(phimax)+ i*ll*sin(phimax);
siamax = subintarg(zzmax, l, x, y, x1, x2, x3, x4);

zzopt = ll*cos(phiopt)+ i*ll*sin(phiopt);
siaopt = subintarg(zzopt, l, x, y, x1, x2, x3, x4);

zzfin = ll*cos(phi)+ i*ll*sin(phi);
siafin = subintarg(zzfin, l, x, y, x1, x2, x3, x4);

figure('Name', 'contours');
plot(ll, real(siamin), 'r', ll, imag(siamin), 'b', ll, real(siamax), 'r--', ll, real(siaopt), 'r:', ll, imag(siaopt), 'b:', ll, real(siafin), 'r.', ll, imag(siafin), 'b.');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tmax = 3.0e1

l = 0
m = 0
n = 0
j1 = 0
j2 = 0

phi =   0.39269908169872414 


% contour:
ll = linspace(0,tmax,10000);

zz = ll*cos(phi)+ i*ll*sin(phi);

sia_gen = subintarggen(zz, l, m, n, j1, j2, x1, x2, x3, x4);

llmax = - 40.0 / real(arg*(cos(phi)+ i*sin(phi)))

zzmax = llmax*(cos(phi)+ i*sin(phi));

figure('Name', 'contour gen');

subplot(1,2,1);
plot(real(zz), imag(zz), 'b', real(zzmax), imag(zzmax), 'ro');

subplot(1,2,2); hold on;
plot(ll, real(sia_gen), 'r', ll, imag(sia_gen), 'b', llmax, 0, 'ro');

siarggen = zeros(kmax,kmax);
for k1=1:kmax
    for k2=1:kmax

        siarggen(k1,k2) = subintarggen(tau_z(k1)+tau_z(k2)*i, l, m, n, j1, j2, x1, x2, x3, x4);

    end
end

figure;
subplot(1,2,1); hold on;
mesh(tau_z, tau_z, real(siarggen'));
%plot3(0,0,0,'*');
view(3);

subplot(1,2,2); hold on;
mesh(tau_z, tau_z, imag(siarggen'));
view(3);

end

function w = subint(tau, l, x, y, x1, x2, x3, x4)

alpha = sqrt(0.25e0 - x3*i) - 0.5e0;

gamma = alpha / (1.0 + 2.0e0*alpha);

w = -sqrt(2.0*pi) ./ (1.0 + i*x4*tau).^(l + 1) ./ sqrt((gamma - 1.0)^2 - gamma^2 * exp(-2.0*tau)) .* ...
     exp((-x1*x1 + i*x2)*tau + (1.0 - exp(-tau))./(gamma - 1.0 + gamma*exp(-tau)) * x1 * ...
     ( (2.0*gamma - 1.0)*x1 + (x + y)*i) + (gamma*(x*x + y*y)*exp(-2.0*tau) + 2*x*y*exp(-tau) - (gamma - 1.0)*(x*x + y*y)) ./ ...
     2.0 ./ ((gamma - 1.0)^2 - gamma^2*exp(-2.0*tau)));

end

function w = subintlim(tau, l, x, y, x1, x2, x3, x4)

alpha = sqrt(0.25e0 - x3*i) - 0.5e0;

gamma = alpha / (1.0 + 2.0e0*alpha);

w = -sqrt(2.0*pi) ./ (1.0 + i*x4*tau).^(l + 1) ./ sqrt((gamma - 1.0)^2) .* ...
     exp((-x1*x1 + i*x2)*tau + (1.0)./(gamma - 1.0) * x1 * ...
     ( (2.0*gamma - 1.0)*x1 + (x + y)*i) - (x*x + y*y) ./ 2.0 ./ (gamma - 1.0));

end

function w = subintarg(tau, l, x, y, x1, x2, x3, x4)

alpha = sqrt(0.25e0 - x3*i) - 0.5e0;

gamma = alpha / (1.0 + 2.0e0*alpha);

w = -(l+1)*log(1.0 + i*x4*tau) - 0.5*log((gamma - 1.0)^2 - gamma^2 * exp(-2.0*tau)) + ...
     (-x1*x1 + i*x2)*tau + (1.0 - exp(-tau))./(gamma - 1.0 + gamma*exp(-tau)) * x1 * ...
     ( (2.0*gamma - 1.0)*x1 + (x + y)*i) + (gamma*(x*x + y*y)*exp(-2.0*tau) + 2*x*y*exp(-tau) - (gamma - 1.0)*(x*x + y*y)) ./ ...
     2.0 ./ ((gamma - 1.0)^2 - gamma^2*exp(-2.0*tau));

end


function arg = subintarggen(tau, ll, mm, nn, j1, j2, x1, x2, x3, x4)

alpha = sqrt(0.25e0 - x3*i) - 0.5e0;

gamma = alpha / (1.0 + 2.0e0*alpha);

gr = -x1*x1 + i*x2;

expon = exp(-tau);

onepluexp = 1.0e0 + expon;
oneminexp = 1.0e0 - expon;

denplu = 1.0e0 - gamma * onepluexp;
denmin = 1.0e0 - gamma * oneminexp;

fac1 = onepluexp ./ denplu;
fac2 = oneminexp ./ denmin;
fac3 = oneminexp ./ denplu;

arg = - (ll+1) * log(1.0e0 + x4 * tau * i) - 0.5 * log(denmin .* denplu) + ...
        j1 * log(fac1) + j2 * log(fac2) + (mm + nn - 2 * (j1 + j2)) * log(fac3) + ...
        gr * tau + fac3 * (1.0e0 - 2.0 * gamma) * x1 * x1;

end


function ans = subintgen(tau, ll, mm, nn, j1, j2, x1, x2, x3, x4)

alpha = sqrt(0.25e0 - x3*i) - 0.5e0;

gamma = alpha / (1.0 + 2.0e0*alpha);

gr = -x1*x1 + i*x2;

expon = exp(-tau);

onepluexp = 1.0e0 + expon;
oneminexp = 1.0e0 - expon;

denplu = 1.0e0 - gamma * onepluexp;
denmin = 1.0e0 - gamma * oneminexp;

fac1 = onepluexp ./ denplu;
fac2 = oneminexp ./ denmin;
fac3 = oneminexp ./ denplu;

ans = 1.0e0 ./ (1.0e0 + x4 * tau * i).^(ll+1) ./ sqrt(denmin .* denplu) .* ...
      fac1.^j1 .* fac2.^j2 .* fac3.^(mm + nn - 2 * (j1 + j2)) .* ...
      exp( gr * tau + fac3 * (1.0e0 - 2.0 * gamma) * x1 * x1 );

end

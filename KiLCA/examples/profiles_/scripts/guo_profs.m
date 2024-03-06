global a alpha theta_0 mu_0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mu_0 = 4.0*pi*1.0E-7;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%path_out = 'profiles_guo_1/';  %corresponds to p4
path_out = 'profiles_guo_2/';  %corresponds to p4

%path_out = 'profiles_guo_p5/';  %rtor/a parameter

par = 4;

%case settings in SI units:

a = 0.5;

%guo_1:
%  alpha = 9.66;  %5.5;
%  theta_0 = 1.5; %1.6;

%guo_2:
alpha = 5.45;  %5.5;
theta_0 = 1.6; %1.6;

%rtor = 4*a;    %???
rtor = par*a;   %???

Bth_0  = 0.0;
Bz_0 = 0.5;    %???
p_0 = 0.0;

%below constants are in Gauss units:
m_i = 2.0*1.67262158E-24;

n_0 = 1.0e13; %???
Ti_0 = 0.0;
Te_0 = 0.0;

Vth_0 = 0.0;
%Vz_0 = 0.0;

Bth_a = 1.878101216715334E3; %Bth at the edge
Vth_a = Bth_a/sqrt(4.0*pi*m_i*n_0);
Vz_0 = 1.0*Vth_a;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mkdir (path_out);

rb = 0.0:0.001:a;

len = length(rb);

%computation of u = {Bth, Bz, p}:
options = odeset ('RelTol', 1e-12);
[rb, u_ode] = ode113 (@u_rhs, rb, [rb(1)*Bth_0 Bz_0 p_0], options);

Bth = u_ode(:,1)./rb; Bth(1) = Bth_0;
Bz = u_ode(:,2);
p = u_ode(:,3);

fB = figure('Name', 'B'); hold on;
plot (rb, Bth, 'r', rb, Bz, 'g');

fp = figure('Name', 'press'); hold on;
plot (rb, p);

q = (rb/rtor).*(Bz./Bth); q(1) = a/rtor/theta_0;

fp = figure('Name', 'q'); hold on;
plot (rb, q, 'r', rb, q*rtor/a, 'g');

%current densities:
mu = 2.0/a * theta_0 * (1.0 - (rb/a).^alpha);
Jth = mu/mu_0 .* Bth;
Jz  = mu/mu_0 .* Bz;

fp = figure('Name', 'J'); hold on;
plot (rb, Jth, 'r', rb, Jz, 'g');

%rescaling to SGS:
ind = 2:length(rb);

%q
Q = [100*rb(ind) q(ind)];
fid = fopen([path_out 'q.dat'],'wb');
fprintf(fid,'%.15e\t%.15e\n', Q');
fclose(fid);
figure('Name', 'q'); hold on;
plot (Q(:,1), Q(:,2), 'r');

%n
Q = [100*rb(ind) n_0*ones(length(ind),1)];
fid = fopen([path_out 'n.dat'],'wb');
fprintf(fid,'%.15e\t%.15e\n', Q');
fclose(fid);
figure('Name', 'n'); hold on;
plot (Q(:,1), Q(:,2), 'r');

%Ti
Q = [100*rb(ind) Ti_0*ones(length(ind),1)];
fid = fopen([path_out 'Ti.dat'],'wb');
fprintf(fid,'%.15e\t%.15e\n', Q');
fclose(fid);
figure('Name', 'Ti Te'); hold on;
plot (Q(:,1), Q(:,2), 'r');

%Te
Q = [100*rb(ind) Te_0*ones(length(ind),1)];
fid = fopen([path_out 'Te.dat'],'wb');
fprintf(fid,'%.15e\t%.15e\n', Q');
fclose(fid);
plot (Q(:,1), Q(:,2), 'b');

%Vth
Q = [100*rb(ind) Vth_0*ones(length(ind),1)];
fid = fopen([path_out 'Vth.dat'],'wb');
fprintf(fid,'%.15e\t%.15e\n', Q');
fclose(fid);
figure('Name', 'Vth Vz'); hold on;
plot (Q(:,1), Q(:,2), 'r');

%Vz
Q = [100*rb(ind) Vz_0*ones(length(ind),1)];
fid = fopen([path_out 'Vz.dat'],'wb');
fprintf(fid,'%.15e\t%.15e\n', Q');
fclose(fid);
plot (Q(:,1), Q(:,2), 'b');

%Er
Q = [100*rb(ind) 0.0*ones(length(ind),1)];
fid = fopen([path_out 'Er.dat'],'wb');
fprintf(fid,'%.15e\t%.15e\n', Q');
fclose(fid);
figure('Name', 'Er'); hold on;
plot (Q(:,1), Q(:,2), 'r');

%Bth
Q = [100*rb(ind) 1.0E4*Bth(ind)];
fid = fopen([path_out 'Bth.dat'],'wb');
fprintf(fid,'%.15e\t%.15e\n', Q');
fclose(fid);
figure('Name', 'Bth Bz'); hold on;
plot (Q(:,1), Q(:,2), 'r');

%Bz
Q = [100*rb(ind) 1.0E4*Bz(ind)];
fid = fopen([path_out 'Bz.dat'],'wb');
fprintf(fid,'%.15e\t%.15e\n', Q');
fclose(fid);
plot (Q(:,1), Q(:,2), 'b');

%Jth
Q = [100*rb(ind) 3.0E5*Jth(ind)];
fid = fopen([path_out 'Jth.dat'],'wb');
fprintf(fid,'%.15e\t%.15e\n', Q');
fclose(fid);
figure('Name', 'Jth Jz'); hold on;
plot (Q(:,1), Q(:,2), 'r');

%Jz
Q = [100*rb(ind) 3.0E5*Jz(ind)];
fid = fopen([path_out 'Jz.dat'],'wb');
fprintf(fid,'%.15e\t%.15e\n', Q');
fclose(fid);
plot (Q(:,1), Q(:,2), 'b');

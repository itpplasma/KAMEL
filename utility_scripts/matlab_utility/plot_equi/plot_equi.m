% plot_equi.m
% ----------------------------------------------------------------------------
% This script is intended to plot the equilibrium flux surface of a given 
% shot/time pair.
% ----------------------------------------------------------------------------

% author: Markus Markl, with usage of code written by Philipp Ulbl
% created: 21.07.2021

% ----------------------------------------------------------------------------

PRINT = true;

addpath('~/Dokumente/plasma/code/libneo/matlab/EFIT');

equi_path = '/temp/markl_m/BALANCE_4CASES/DATA/EQUI/33353/g33353.3200_EQH';

e = efit(equi_path, [], []);
e.read();

phi = linspace(0, 1.5*pi, 100);
r = e.rbbbs .* 100;
z = e.zbbbs .* 100;

col = 0.7 .* [0, 1, 1];

[R, PHI] = meshgrid(r, phi);
XB = R .* cos(PHI);
YB = R .* sin(PHI);
ZB = repmat(z, numel(phi), 1);

fh = figure('units', 'normalized', 'outerposition', [0, 0, 1.25, 1.0]);

subplot(1,2,1)
daspect([1 1 1])
sb = surf(XB, YB, ZB, 'FaceColor', col, 'EdgeColor', 'none', 'FaceAlpha', ...
	0.4, 'DisplayName', 'Plasma boundary');
patch(reducepatch(surf2patch(sb), 1e-6), 'EdgeColor', 'none', ...
	'HandleVisibility', 'Off');
set(sb, 'FaceLighting', 'gouraud');
light

% View of 3D plot
az = 24.2848;
el = 50.0;%60.7185;
view(az, el);

xl = xlim;
yl = ylim;
zl = zlim;

%axis off

% Torus
subplot(1,2,2)
daspect([1 1 1])

major_r = 1.65*100;
minor_r = 0.5*100; %centimeters


tor_p = linspace(0, 2*pi, 100);
tor_t = linspace(0, 1.5*pi, 100);

[PP, TT] = meshgrid(tor_p, tor_t);

tor_x = (major_r + minor_r .* cos(PP)) .* cos(TT);
tor_y = (major_r + minor_r .* cos(PP)) .* sin(TT);
tor_z = minor_r * sin(PP);

tor_s = surf(tor_x, tor_y, tor_z, 'FaceColor', col, 'EdgeColor', 'none', 'FaceAlpha', ...
	0.4, 'DisplayName', 'Plasma boundary');
patch(reducepatch(surf2patch(tor_s), 1e-6), 'EdgeColor', 'none', ...
	'HandleVisibility', 'Off');
set(tor_s, 'FaceLighting', 'gouraud');
light

view(az, 30.0);

%axis off
xlim(xl);
ylim(yl);
zlim(zl);



%tor_sb = surf(, YB, ZB, 'FaceColor', col, 'EdgeColor', 'none', 'FaceAlpha', ...


if PRINT == true
	system(['mkdir -p out'])
	print('./out/equitorus.png', '-dpng', '-r200');
end

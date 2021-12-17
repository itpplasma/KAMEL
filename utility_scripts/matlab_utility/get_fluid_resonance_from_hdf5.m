%##############################################################################
% get_fluid_resonance_from_hdf5.m
%##############################################################################
% description:
%------------------------------------------------------------------------------
% This function is used to determine the position of the fluid resonance
% in relation to the resonant surfaces.
%##############################################################################

% author: Markus Markl
% created: 24.11.2021

function [pos, res_surf, time] = get_fluid_resonance_from_hdf5(shot, time, datapath)
	% pos ... radial location of the fluid resonance, is a vector if time is a vector
	% res_surf ... radial location of the resonant surfaces for m=5,6,7 (n=2)

	echarge = 4.8e-10;
	kB = 1.3807e-16;
	eVK = 1.1604e4;

	if (nargin<3 || isempty(datapath))
		datapath = ['/temp/markl_m/BALANCE_4CASES/PRERUNDATA/', num2str(shot), '/'];
	end

	pos = zeros(numel(time), 1);
	res_surf = zeros(3, numel(time));
	

%	for ti=1:numel(time)
	ti = 1;
	while true
		if ti > numel(time)
			break
		end
		fname= [datapath, num2str(shot), '_', num2str(time(ti)),'_mi_2.hdf5'];

		try
			r_v = h5read(fname, '/KiLCA_flre/output/postprocessor1/r');
		catch
			disp(['no data for ', num2str(time(ti))]);
			time(ti) = []; % delete time index, since no data is available
			res_surf(:,ti) = [];
			pos(ti) = [];
			continue
		end
		v_ed = h5read(fname, '/KiLCA_flre/output/postprocessor1/ved');
		v_ExB = h5read(fname, '/KiLCA_flre/output/postprocessor1/vExB');

		r_prof = h5read(fname, '/preprocprof/r_out');
		temp = interp1(r_prof, h5read(fname, '/preprocprof/Te'), r_v) .*kB .*eVK;
		grad_temp = gradient(temp,r_v);
		dens = interp1(r_prof, h5read(fname, '/preprocprof/n'), r_v);
		grad_dens = gradient(dens, r_v);

		r_b = h5read(fname, '/KiLCA_flre/output/background/R');
		b0 = interp1(r_b, h5read(fname, '/KiLCA_flre/output/background/b0'), r_v);
		r_res = h5read(fname, '/output/r_res');
		res_surf(:,ti) = r_res;
		
		v_fluid_wo_gradTe = v_ExB + v_ed - (1/echarge) .* (grad_temp) ./ b0;
		%zero_crossing = @(v) find(diff(sign(v)));
		%r_fluid_resonance = r_v(zero_crossing(v_fluid_wo_gradTe));
		r_fluid_resonance = find_resonance(v_fluid_wo_gradTe, r_v, r_res);

		%if numel(r_fluid_resonance) > 1
		%	[a, ind1] = min(abs(r_res - r_fluid_resonance), [], 2);
		%	[b, ind2] = min(a, [], 1);
		%	r_fluid_resonance = r_fluid_resonance(ind2);
		%end

		pos(ti) = r_fluid_resonance;
		ti = ti+1;
	end

end % get_fluid_resonance_from_hdf5


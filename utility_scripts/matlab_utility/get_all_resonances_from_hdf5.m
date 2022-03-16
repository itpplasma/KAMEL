%##############################################################################
% get_all_resonances_from_hdf5.m
%##############################################################################
% description:
%------------------------------------------------------------------------------
% This function is used to determine the position of the fluid resonance
% in relation to the resonant surfaces.
%##############################################################################

% author: Markus Markl
% created: 30.11.2021

function [fl_wo_dTe, fl_w_dTe, el_res, res_surf, time] = get_all_resonances_from_hdf5(shot, time, outpath,datapath)
	% fl_wo_dTe ... radial location of the fluid resonance without electron temperature gradient, is a vector if time is a vector
	% fl_w_dTe ... radial location of the fluid resonance with electron temperature gradient, is a vector if time is a vector
	% res_surf ... radial location of the resonant surfaces for m=5,6,7 (n=2)

	echarge = 4.8e-10;
	kB = 1.3807e-16;
	eVK = 1.1604e4;

	if (nargin<4 || isempty(datapath))
		datapath = ['/temp/markl_m/BALANCE_4CASES/PRERUNDATA/', num2str(shot), '/'];
	end

	fl_wo_dTe = zeros(numel(time), 1);
	fl_w_dTe = zeros(numel(time), 1);
	el_res = zeros(numel(time), 1);
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
			fl_w_dTe(ti) = [];
			fl_wo_dTe(ti) = [];
			el_res = [];
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
		[r_fluid_resonance,more] = find_resonance(v_fluid_wo_gradTe, r_v, r_res);
		disp([num2str(shot), '_', num2str(time(ti)), ' fluid more: ', num2str(more)])
		fileID = fopen([outpath, num2str(shot), '_fluid_resonance.dat'], 'a');
		write_resonance_to_file(time(ti), r_fluid_resonance, fileID);
		%fl_wo_dTe(ti) = r_fluid_resonance;

		%if numel(r_fluid_resonance) > 1
		%	[a, ind1] = min(abs(r_res - r_fluid_resonance), [], 2);
		%	[b, ind2] = min(a, [], 1);
		%	r_fluid_resonance = r_fluid_resonance(ind2);
		%end

%		fl_wo_dTe(ti) = r_fluid_resonance;


		[res_val ,more] = find_resonance(v_ExB + v_ed, r_v, r_res);
		fileID = fopen([outpath, num2str(shot),'_fluid_resonance_wdTe.dat'], 'a');
		write_resonance_to_file(time(ti), res_val, fileID);
		%fl_w_dTe(ti) = res_val;


		disp([num2str(shot), '_', num2str(time(ti)), ' fluid wdTe more: ', num2str(more)])

		[res_val, more] = find_resonance(v_ExB, r_v, r_res);
		fileID = fopen([outpath, num2str(shot), '_electric_resonance.dat'], 'a');
		write_resonance_to_file(time(ti), res_val, fileID);
		%el_res(ti) = res_val;

		disp([num2str(shot), '_', num2str(time(ti)), ' electric more: ', num2str(more)])

		fileID = fopen([outpath, num2str(shot), '_resonant_surfaces.dat'], 'a');
		write_to_file(time(ti), r_res, fileID);

		ti = ti+1;
	end

end % get_fluid_resonance_from_hdf5


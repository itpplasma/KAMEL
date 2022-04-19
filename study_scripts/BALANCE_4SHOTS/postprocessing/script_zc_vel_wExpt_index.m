%##########################################################################
% script_zc_vel_wExpt_index.m
%##########################################################################
% description:
%--------------------------------------------------------------------------
% This script plots the location of the fluid resonance relative to the 
% experimental diagnostics. This is done by employing an index function
% that creates a .index or .indexh5 file which lists the available data.
%##########################################################################
% Note: due to the change of the ql-balance code to the hdf5 file format for
% reading and writing, two different hdf5 files are invoked.
% If functions are used that say something like "*_from_hdf5", data from
% a new version run is used.
     
%author:   Markus Markl
%created:  24.11.2021

PRINT = false;
COPY = false;

addpath('~/Dokumente/plasma/code/libneo_old/matlab/Utility/xxaxis');
addpath('/temp/markl_m/GITHUB/Balance/utility_scripts/matlab_utility/');
addpath('./Plot');

shot = 33353;

prerundatapath = ['/temp/markl_m/BALANCE_4CASES/PRERUNDATA/',num2str(shot),'/'];
outpath = './result/';

system(['rm ', outpath, num2str(shot), '_fluid_resonance.dat']);
system(['rm ', outpath, num2str(shot), '_fluid_resonance_wdTe.dat']);
system(['rm ', outpath, num2str(shot), '_electric_resonance.dat']);
system(['rm ', outpath, num2str(shot), '_resonant_surfaces.dat']);



%create_index(shot);
create_index_from_hdf5(shot, prerundatapath);
%[time, studyname] = get_shot_times(shot);
timeh5 = get_shot_times_from_hdf5(shot);

%[fluid_resonance_h5, res_surf_h5, timeh5] = get_fluid_resonance_from_hdf5(shot, timeh5);
[fluid_resonance_h5, fluid_resonance_wdTe_h5, electric_resonance_h5, res_surf_h5, timeh5] = get_all_resonances_from_hdf5(shot, timeh5, outpath);

zcvec = zeros(numel(timeh5), 1); % fluid resonance without Te gradient
%zcvec = zeros(numel(time)+numel(timeh5), 1); % fluid resonance without Te gradient
electric_res = zeros(numel(timeh5), 1);% electric resonance v_ExB = 0
%electric_res = zeros(numel(time)+numel(timeh5), 1);% electric resonance v_ExB = 0
fl_res_wdTe = zeros(numel(timeh5), 1); % fluid resonance with Te gradient
%fl_res_wdTe = zeros(numel(time)+numel(timeh5), 1); % fluid resonance with Te gradient
ressurf = zeros(3,numel(timeh5));
%ressurf = zeros(3,numel(time)+numel(timeh5));

%zcvec(numel(time)+1:numel(time)+numel(timeh5)) = fluid_resonance_h5;
%electric_res(numel(time)+1:numel(time)+numel(timeh5)) = electric_resonance_h5;
%fl_res_wdTe(numel(time)+1:numel(time)+numel(timeh5)) = fluid_resonance_wdTe_h5;
%ressurf(:,numel(time)+1:numel(time)+numel(timeh5)) = res_surf_h5;

runname = [num2str(shot)];

%datapath = ['/temp/ulbl_p/AUG/SHOTS/', num2str(shot), '/'];
%outpath = ['/temp/markl_m/'];

%Ipolsoli = [num2str(shot), '_MAC_Ipolsoli.dat'];
%IBl6 = [num2str(shot), '_MAW_IBl6.dat'];
%ELMi = [num2str(shot), '_POT_ELMi-Han.dat'];

%Ipolsoli_data = importdata([datapath, Ipolsoli]);
%IBl6_data = importdata([datapath, IBl6]);
%ELMi_data = importdata([datapath, ELMi]);

%I_ymax = 5e4;
%Da_ymax = 1;

%outpath = ['/temp/markl_m/BALANCE_4CASES/POST/velscans/', num2str(shot), '/'];


% go through old version of balance code
for i = 1:numel(time)
	%disp(num2str(time(i)))
	runname = [num2str(shot), '_', num2str(time(i))];

	fname = ['/temp/markl_m/BALANCE_4CASES/POST/', studyname{i}, '/', runname, '/', studyname{i}, '_NRef_', runname, '.hdf5'];


	echarge = 4.8e-10;   % electron charge in cgs
	kB = 1.3807e-16;     % Boltzmann constant
	eVK = 1.1604e4;      % eV -> deg(K)

	try
		r_v = h5read(fname, '/KiLCA_flre/output/postprocessor1/r');
	catch
		disp(['no data for ', num2str(time(i))]);
		continue
	end

	v_ed = h5read(fname, '/KiLCA_flre/output/postprocessor1/ved');

	v_ExB = h5read(fname, '/KiLCA_flre/output/postprocessor1/vExB');

	r_prof = h5read(fname, '/profiles/r_out');
	temp = interp1(r_prof, h5read(fname, '/profiles/Te'), r_v) .* kB .* eVK;
	grad_temp = gradient(temp, r_v);
	dens = interp1(r_prof, h5read(fname, '/profiles/n'), r_v);
	grad_dens = gradient(dens, r_v);

	r_b = h5read(fname, '/KiLCA_flre/output/background/R');
	b0 = interp1(r_b, h5read(fname, '/KiLCA_flre/output/background/b0'), r_v);
	r_res = h5read(fname, '/output/r_res');
	ressurf(:, i) = r_res;

	v_ed_nodTedr = v_ed - (1 / echarge) .* (grad_temp) ./ b0;
	v_ed_nodndr = v_ed - (1 / echarge) .* (temp .* grad_dens ./ dens) ./ b0;

	v_final = v_ExB + v_ed - (1/echarge) .* (grad_temp) ./b0;
	%zci = @(v) find(diff(sign(v)));
	%r_zc = r_v(zci(v_final));
	[r_zc, more] = find_resonance(v_final, r_v, r_res);
	disp([runname, ' more: ', num2str(more)]);

	%if numel(r_zc) > 1
		%[a, ind1] = min(abs(r_res - r_zc), [], 2);
		%[b, ind2] = min(a,[],1);
		%r_zc = r_zc(ind2);
	%end

	%zcvec(i) = r_zc;
	fileID = fopen([outpath, 'fluid_resonance.dat'],'a');
	write_resonance_to_file(time(i), r_zc, fileID);

	[r_zc, more] = find_resonance(v_ExB, r_v, r_res);
	fileID = fopen([outpath, 'electric_resonance.dat'],'a');
	write_resonance_to_file(time(i), r_zc, fileID);
	%electric_res(i) = r_zc;

	[r_zc, more] = find_resonance(v_ExB + v_ed, r_v, r_res);
	%fl_res_wdTe(i) = r_zc;
	fileID = fopen([outpath, 'fluid_resonance_wdTe.dat'],'a');
	write_resonance_to_file(time(i), r_zc, fileID);

	
	fileID = fopen([outpath, 'resonant_surfaces.dat'],'a');
	write_to_file(time(i), r_res, fileID);


end % end of loop over time slices

%[time, sortI] = sort([time;timeh5]);
%zcvec = zcvec(sortI);
%ressurf = ressurf(:,sortI);

%fileID = fopen([outpath, 'resonant_surfaces.dat'],'w');
%fprintf(fileID, '%u %u %u %u\n', [time'; ressurf]);
%fclose(fileID);





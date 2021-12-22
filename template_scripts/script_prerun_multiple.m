%#########################################################################
% script_prerun_multiple.m
%##########################################################################
% description:
%--------------------------------------------------------------------------
% This script runs the prerun to the balance code for multiple time slices 
% for one shot.
% The prerun comprises the calculation of the Fourier modes, the profile 
% preprocessor, the KiLCA vacuum and flre run and the GPEC run.
%
%##########################################################################

%author:   Markus Markl
%created:  22.12.2021

ion_mass = 2; % mass is ion_mass * proton mass

% path to KiLCA interface and balance class. Needs to be changed individually.
libBalance = '/temp/markl_m/GITHUB/BalanceSuite/matlab/balance';

addpath(genpath(libBalance))

mpath = pwd();

%Runs to make
shot = 33353; 
time = [1250, 1300, 1350, 1400, 1450, 1500, 1550, 1600, 1650, 1700, 1750, 1800, 1850,1900,1950,2000, 2100, 2150, 2200, 2250,2300, 2350, 2400, 2450,2500,2550,2600,2650,2700,2750,2770,2850,2900,2950,3000,3050,3100,3150,3200];  
% if profiles have an offset in time in the file names, use this:
profoffset = [-2, -9, -3, 0, 1, -1, 0, 0, 1, 0, 1, 4, 0,4,1,0, 3,-3,0, 3, -4,-6,-1,-1,2,0,-4,3,1,2,1,0,0,0,0,0,0,0,0];

% if the g-file has an offset in time, use this:
equioffset = 0;

proftype = '_MMARKL_rho_pol';
studyname = ['BALANCE_4CASES', proftype];
% TODO: change the basepath to your needs
basepath = '/temp/markl_m/BALANCE_4CASES/';
datapath = [basepath, 'DATA/'];

runname = ['prerun'];
m = [5,6,7];
n = 2 .* ones(size(m));

% iterate over all time slices
for i=1:numel(time)
	tstart = tic;
	shotname = [num2str(shot), '_', num2str(time(i))];

	path2preh5 = [basepath,'PRERUNDATA/',num2str(shot), '/', num2str(shot), '_', num2str(time(i)), '_mi_', num2str(ion_mass), '.hdf5'];

	runpath = [basepath, 'RUNS_PRE/', studyname,'/', num2str(shot), '_', num2str(time(i)),'/',runname,'/'];

	%input
	gfile  = [datapath, 'EQUI/', num2str(shot),'/g',num2str(shot),'.',num2str(time(i)+equioffset),'_EQH'];
	%gfile  = [datapath, 'DATA/EQUI/', num2str(shot),'/g',num2str(shot),'.',num2str(time),'_COCOS3'];
	cfile  = [datapath, 'COIL/', num2str(shot),'/',num2str(shot),'.',num2str(time(i)),'_coil.dat'];
	%dapath = '/temp/markl_m/DA/ASTRA/';
	dapath = ['/temp/markl_m/AUG/SHOTS/', num2str(shot), '/'];

	filehead = [datapath, 'PROF/', num2str(shot),'/',num2str(shot),'.',num2str(time(i)+profoffset(i))];
	neprof = [filehead,'_ne_PED',proftype,'.dat'];
	Teprof = [filehead,'_Te_PED',proftype,'.dat'];
	Tiprof = [filehead,'_Ti_PED',proftype,'.dat'];
	vtprof = [filehead,'_vt_PED',proftype,'.dat'];
	copypath = ['/temp/markl_m/BALANCE_4CASES/RUNS_PRE/BALANCE_4CASES_MMARKL_rho_pol/', num2str(shot), '_', num2str(time(i)),'/prerun/profiles/'];


%TODO: change accordingly
	%location of field.dat
	pfile = ['/temp/markl_m/MESH3D/',num2str(shot),'/field.dat']; %will be calculated if not present
	fluxdatapath = ['/temp/markl_m/FLUXDATA/',num2str(shot),'/',num2str(time(i)),'/']; %will be calculated if not present
	gpecpath = ['/temp/markl_m/GPEC/', studyname, '/', num2str(shot), '_', num2str(time(i)), '/'];

	% Create balance object
	bal = Balance(runpath, shot, time(i), [studyname, '_', runname], path2preh5);
	bal.setModes(m, n);
	try
		bal.setCoil(cfile, pfile);
	catch
		fileID = fopen('log.txt','a');
		fprintf(fileID, ['Problem in coil for ', num2str(time(i)),'; ', datestr(datetime('now')), ' \n']);
		fclose(fileID);
		continue
	end
	try
		bal.setEqui(gfile, fluxdatapath);
	catch
		fileID = fopen('log.txt','a');
		fprintf(fileID, ['Problem in equi for ', num2str(time(i)), '; ', datestr(datetime('now')),' \n']);
		fclose(fileID);
		continue
	end
	try
		bal.setTMHDCode('GPEC', gpecpath);
	catch
		fileID = fopen('log.txt','a');
		fprintf(fileID, ['Problem in GPEC for ', num2str(time(i)), '; ', datestr(datetime('now')),' \n']);
		fclose(fileID);
		continue
	end
	try
		bal.setProfiles(neprof, Teprof, Tiprof, vtprof, copypath);
	catch
		fileID = fopen('log.txt','a');
		fprintf(fileID, ['Problem in profiles for ', num2str(time(i)), '; ', datestr(datetime('now')),' \n']);
		fclose(fileID);
		continue
	end

	bal.setKiLCA(ion_mass);
	bal.setDaEstimation(dapath);
	bal.write();

	system(['mkdir -p ', basepath, 'POST/', studyname, '/', shotname, '/']);
%bal.export2HDF5([basepath, 'POST/', studyname, '/', shotname, '/'], [studyname, '_', runname]);

	cd(mpath);

	% link the resulting hdf5 files, if you want
	%system(['ln -sf ', path2preh5, ' ', '/temp/markl_m/inputdata_h5/', num2str(shot),'/',...
	%num2str(shot), '_', num2str(time(i)), '.hdf5']);

	disp(['Finished! :), shot: ', num2str(shot), ' at time: ', num2str(time(i)), ' with ion mass= ', num2str(ion_mass)]);

	tend = toc(tstart);

	fileID = fopen('log.txt','a');
	fprintf(fileID, '------------------------------------------------------------------------------\n');
	fprintf(fileID, ['Finished ', num2str(time(i)), ' successfully at ', datestr(datetime('now')),' \n']);
	fprintf(fileID, ['It took ', num2str(tend), ' seconds\n']);
	fprintf(fileID, '------------------------------------------------------------------------------\n\n');
	fclose(fileID);

end

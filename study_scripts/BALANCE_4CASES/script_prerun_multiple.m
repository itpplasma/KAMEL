%#########################################################################
% script_prerun.m
%##########################################################################
% description:
%--------------------------------------------------------------------------
% This script runs the prerun to the balance code for given shot/time pair.
% The prerun comprises the calculation of the Fourier modes, the profile 
% preprocessor, the KiLCA vacuum and flre run and the GPEC run.
%
%##########################################################################

%author:   Markus Markl
%created:  28.06.2021

ion_mass = 2; % mass is ion_mass * proton mass

% path to KiLCA interface and balance class. Needs to be changed individually.
libBalance = '/temp/markl_m/GITHUB/Balance/matlab/balance';

addpath(genpath(libBalance))

mpath = [pwd(), '/'];

%Runs to make
shot = 33353; 

timelist = importdata(['/temp/markl_m/BALANCE_4CASES/DATA/',num2str(shot),'_time_list.txt'])

time = timelist(:,1);%[1250, 1300, 1350, 1400, 1450, 1500, 1550, 1600, 1650, 1700, 1750, 1800, 1850,1900,1950,2000, 2100, 2150, 2200, 2250,2300, 2350, 2400, 2450,2500,2550,2600,2650,2700,2750,2770,2850,2900,2950,3000,3050,3100,3150,3200];  
tcut = 3000;
profoffset = timelist(:,2);%[-2, -9, -3, 0, 1, -1, 0, 0, 1, 0, 1, 4, 0,4,1,0, 3,-3,0, 3, -4,-6,-1,-1,2,0,-4,3,1,2,1,0,0,0,0,0,0,0,0];
profoffset = profoffset(time > tcut);
equioffset = timelist(:,3);
equioffset = equioffset(time > tcut);

time = time(time >tcut);

proftype = '_MMARKL_rho_pol';
studyname = ['BALANCE_4CASES', proftype];
basepath = '/temp/markl_m/BALANCE_4CASES/';
datapath = [basepath, 'DATA/'];

runname = ['prerun'];
m = [5,6,7];
n = 2 .* ones(size(m));

for i=1:numel(time)
	disp('==== ==== ==== ==== ==== ==== ==== ==== ====')
	disp(['Processing now: shot ', num2str(shot), ' at time ', num2str(time(i))])
	disp(' ')
	tstart = tic;
	shotname = [num2str(shot), '_', num2str(time(i))];

	path2preh5 = [basepath,'PRERUNDATA/',num2str(shot), '/', num2str(shot), '_', num2str(time(i)), '_mi_', num2str(ion_mass), '.hdf5'];

	runpath = [basepath, 'RUNS_PRE/', studyname,'/', num2str(shot), '_', num2str(time(i)),'/',runname,'/'];

	%input
	gfile  = [datapath, 'EQUI/', num2str(shot),'/g',num2str(shot),'.',num2str(time(i)+equioffset(i)),'_EQH'];
	%gfile  = [datapath, 'DATA/EQUI/', num2str(shot),'/g',num2str(shot),'.',num2str(time),'_COCOS3'];
	cfile  = [datapath, 'COIL/', num2str(shot),'/',num2str(shot),'.',num2str(time(i)),'_coil.dat'];
	%dapath = '/temp/markl_m/DA/ASTRA/';
	dapath = ['/temp/markl_m/AUG/SHOTS/', num2str(shot), '/'];

	filehead = [datapath, 'PROF/', num2str(shot),'/',num2str(shot),'.',num2str(time(i)+profoffset(i))];
	neprof = [filehead,'_ne_PED',proftype,'.dat'];
	Teprof = [filehead,'_Te_PED',proftype,'.dat'];
	Tiprof = [filehead,'_Ti_PED',proftype,'.dat'];
	vtprof = [filehead,'_vt_PED',proftype,'.dat'];
	%copypath = ['/temp/markl_m/BALANCE_4CASES/RUNS_PRE/BALANCE_4CASES_MMARKL_rho_pol/', num2str(shot), '_', num2str(time(i)),'/prerun/profiles/'];


	%location of field.dat
	pfile = ['/temp/markl_m/MESH3D/',num2str(shot),'/field.dat']; %will be calculated if not present
	fluxdatapath = ['/temp/markl_m/FLUXDATA/',num2str(shot),'/',num2str(time(i)),'/']; %will be calculated if not present
	gpecpath = ['/temp/markl_m/GPEC/', studyname, '/', num2str(shot), '_', num2str(time(i)), '/'];

	%BALANCE CODE
	bal = Balance(runpath, shot, time(i), [studyname, '_', runname], path2preh5);
	bal.setModes(m, n);

	disp('==== ==== ==== ==== ==== ==== ==== ==== ====')
	try
		bal.setCoil(cfile, pfile);
	catch
		fileID = fopen([mpath,'log_', num2str(shot), '.txt'],'a');
		fprintf(fileID, ['Problem in coil for ', num2str(time(i)),'; ', datestr(datetime('now')), ' \n']);
		fclose(fileID);
		disp(['Problem in coil for ', num2str(time(i)),'; ', datestr(datetime('now'))]);
		continue
	end

	disp('==== ==== ==== ==== ==== ==== ==== ==== ====')
	try
		bal.setEqui(gfile, fluxdatapath);
	catch
		fileID = fopen([mpath, 'log_', num2str(shot), '.txt'],'a');
		fprintf(fileID, ['Problem in equi for ', num2str(time(i)), '; ', datestr(datetime('now'))]);
		fclose(fileID);
		disp(['Problem in equi for ', num2str(time(i)), '; ', datestr(datetime('now'))]);
		continue
	end

	disp('==== ==== ==== ==== ==== ==== ==== ==== ====')
	try
		bal.setTMHDCode('GPEC', gpecpath);
	catch
		fileID = fopen([mpath, 'log_',num2str(shot), '.txt'],'a');
		fprintf(fileID, ['Problem in GPEC for ', num2str(time(i)), '; ', datestr(datetime('now')),' \n']);
		fclose(fileID);
		disp(['Problem in GPEC for ', num2str(time(i)), '; ', datestr(datetime('now'))]);
		continue
	end

	disp('==== ==== ==== ==== ==== ==== ==== ==== ====')
	try
		%bal.setProfiles(neprof, Teprof, Tiprof, vtprof, copypath);
		bal.setProfiles(neprof, Teprof, Tiprof, vtprof);
	catch
		fileID = fopen([mpath, 'log_', num2str(shot), '.txt'],'a');
		fprintf(fileID, ['Problem in profiles for ', num2str(time(i)), '; ', datestr(datetime('now')),' \n']);
		fclose(fileID);
		disp(['Problem in profiles for ', num2str(time(i)), '; ', datestr(datetime('now'))]);
		continue
	end

	disp('==== ==== ==== ==== ==== ==== ==== ==== ====')

	bal.setKiLCA(ion_mass);
	bal.setDaEstimation(dapath);
	bal.write();

	system(['mkdir -p ', basepath, 'POST/', studyname, '/', shotname, '/']);
%bal.export2HDF5([basepath, 'POST/', studyname, '/', shotname, '/'], [studyname, '_', runname]);

	cd(mpath);

	system(['ln -sf ', path2preh5, ' ', '/temp/markl_m/BALANCE_4CASES/PRERUNDATA/', num2str(shot),'/',...
	num2str(shot), '_', num2str(time(i)), '_mi_', num2str(ion_mass), '.hdf5']);

	disp(['Finished! :), shot: ', num2str(shot), ' at time: ', num2str(time(i)), ' with ion mass= ', num2str(ion_mass)]);

	tend = toc(tstart);

	fileID = fopen([mpath, 'log_', num2str(shot), '.txt'],'a');
	fprintf(fileID, '------------------------------------------------------------------------------\n');
	fprintf(fileID, ['Finished ', num2str(time(i)), ' successfully at ', datestr(datetime('now')),' \n']);
	fprintf(fileID, ['It took ', num2str(tend), ' seconds\n']);
	fprintf(fileID, '------------------------------------------------------------------------------\n\n');
	fclose(fileID);

end

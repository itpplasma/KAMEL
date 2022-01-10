%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% linear_balance_multiple.m
% -----------------------------------------------------------------------
% description:
% -------------------
% This script runs a linear balance run of available time slices of a
% specific shots. To check what time slices are available, an index is
% created from existing prerun hdf5 files.
% It also creates a log.txt file, which can be useful.
% -----------------------------------------------------------------------
% author: Markus Markl
% created: 22.12.2021

% TODO: change paths
libBalance = '/temp/markl_m/GITHUB/Balance/matlab/balance';
addpath('/temp/markl_m/GITHUB/Balance/utility_scripts/matlab_utility/')
addpath(genpath(libBalance))

studyname = 'linearrun'; % specifies the directory, should describe the intent
% project directory name:
project = 'BALANCE_4CASES';

shot = 33353;

% TODO: change path accordingly
prerundatapath = ['/temp/markl_m/', project,'/PRERUNDATA/', num2str(shot), '/'];

create_index_from_hdf5(shot)

timeh5 = get_shot_times_from_hdf5(shot);
%time = [2000];
ionmass = 2;

% specify mode numbers:
m = [5,6,7];
n = 2 .* ones(size(m));

% define factors for parameter scan
fac_n = [1.0]; % for more factors: e.g. [0.5; 1.0; 1.5];
fac_Ti = [1.0];
fac_Te = [1.0];
fac_vz = [1.0];

fileID =fopen('log.txt','a');
fprintf(fileID, ['Starting run at ', datestr(datetime('now')), ' \n', 'Doing ', num2str(numel(timeh5)), ' different time slices\n\n']);
fclose(fileID);


for i =1:numel(timeh5)

	fileID =fopen('log.txt','a');
	fprintf(fileID, ['Currently running ', num2str(timeh5(i)),'\n\n']);
	fclose(fileID);
	tstart = tic;

	% path to the input hdf5 file, that is generated in the pre run
	path2inp = [prerundatapath, num2str(shot), '_', num2str(timeh5(i)),'_mi_', num2str(ionmass),'.hdf5'];

% TODO: change path accordingly
	runpath = ['/temp/markl_m/',project, '/RUNS/', studyname, '/', num2str(shot), '_', num2str(timeh5(i)), '/'];
	system(['mkdir -p ', runpath]);

   
% create balance object
	bal = Balance(runpath, shot, timeh5(i), studyname, path2inp);
	bal.setModes(m, n);
% load KiLCA configuration from prerun
	try
		bal.load_prerun(ionmass);
		bal.write_kilca();
	catch
		fileID =fopen('log.txt','a');
		fprintf(fileID, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
		fprintf(fileID, ['Problem with reading prerun data for time ', num2str(timeh5(i)), ' at ', datestr(datetime('now')), ' \n']);
		fprintf(fileID, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n');
		fclose(fileID);
		continue
	end
% change ion mass to hydrogen
%bal.kil_vacuum.background.mi = 1;
%bal.kil_flre.background.mi = 1;
%bal.write_kilca();


% check if the configuration namelist file exists in run path, 
% if not, copy it from blueprint path
	if ~exist([runpath, 'balance_conf.nml'], 'file')
		system(['cp /temp/markl_m/GITHUB/Balance/matlab/blueprints/balance_conf.nml '...
			, runpath, 'balance_conf.nml']);
	end

	balancenml = InputFile([runpath,'balance_conf.nml']);
	balancenml.read();
	% set options
	balancenml.BALANCENML.path2inp = path2inp;
	balancenml.BALANCENML.flag_run_time_evolution = false;
	balancenml.BALANCENML.br_stopping = true; % stops code, if slope is large enough
	balancenml.BALANCENML.paramscan = false;
	balancenml.BALANCENML.diagnostics_output = false;
	balancenml.BALANCENML.save_prof_time_step = 10;
	balancenml.BALANCENML.debug_mode = true;
	balancenml.BALANCENML.timing_mode = false; % saves time needed for certain computations
	balancenml.BALANCENML.suppression_mode = false;
	balancenml.BALANCENML.am = ionmass;
	%balancenml.BALANCENML.Nstorage = 1000; % 1000 is default
	balancenml.BALANCENML.flre_path = bal.kil_flre.pathofrun;
	balancenml.BALANCENML.vac_path = bal.kil_vacuum.pathofrun;

	bal.setOptionsNML(balancenml);
	%disp('setting factors')
	bal.setFactors(fac_n, 'fac_n')
	bal.setFactors(fac_Ti, 'fac_Ti')
	bal.setFactors(fac_Te, 'fac_Te')
	bal.setFactors(fac_vz, 'fac_vz')
	%disp('now writing');
	%bal.write();
	disp(['now running shot ', num2str(shot),' at time ', num2str(timeh5(i)),'ms']);
	try
		bal.run();
	catch
		fileID =fopen('log.txt','a');
		fprintf(fileID, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
		fprintf(fileID, ['Problem run for time ', num2str(timeh5(i)), ' at ', datestr(datetime('now')), ' \n']);
		fprintf(fileID, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n');
		fclose(fileID);
		continue
	end
	tend = toc(tstart);

	fileID = fopen('log.txt', 'a');
	fprintf(fileID, '-------------------------------------------------\n');
	fprintf(fileID, ['Finished ',num2str(timeh5(i)), ' successfully at ', datestr(datetime('now')), ' \n']);
	fprintf(fileID, ['It took ', num2str(tend), ' seconds\n']);

	fprintf(fileID, '-------------------------------------------------\n\n');
	fclose(fileID);
	system(['mv ', runpath, 'out/', studyname, '_', num2str(shot), '_', num2str(timeh5(i)), '.hdf5 ', '/temp/markl_m/', project,'/POST_h5/', studyname, '/', num2str(shot),'/']);

end




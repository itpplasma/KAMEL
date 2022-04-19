%
%This script is the balance code part of the whole balance cycle.
% This configuration runs the time evolution of a specific shot/time
% created: Markus Markl, 30.03.2021

%libKiLCA = '~/Dokumente/plasma/code/libneo/matlab/BALANCE/balance/KiLCA_interface/';
libBalance = '/temp/markl_m/GITHUB/Balance/matlab/balance';
addpath('/temp/markl_m/GITHUB/Balance/utility_scripts/matlab_utility/')
%addpath(genpath(libKiLCA))
addpath(genpath(libBalance))

mpath = pwd();

studyname = 'linearrun'; % specifies the directory, should describe the intent
% project directory name:
project = 'BALANCE_4CASES';

shot = 33353;
ionmass = 2;

prerundatapath = ['/temp/markl_m/', project,'/PRERUNDATA/', num2str(shot), '/'];

create_index_from_hdf5(shot)
timeh5 = get_shot_times_from_hdf5(shot);
%time = [2000];%shotinput.data(1,2);
timeh5 = timeh5(timeh5>3250);

% specify mode numbers:
m = [5,6,7];
n = 2 .* ones(size(m));

% define factors for parameter scan
fac_n = [1.0];%[0.5; 1.0; 1.5];%[0.3; 0.5; 0.6;0.7; 0.8; 0.9; 1.0; 1.1; 1.2; ...
	%1.3;1.4;1.5;1.8;2.1;2.4;2.7;3.0];
fac_Ti = [1.0];
fac_Te = [1.0];
fac_vz = [1.0];

fileID =fopen([mpath, '/', num2str(shot),'_log.txt'],'a');
fprintf(fileID, ['Starting run at ', datestr(datetime('now')), ' \n', 'Doing ', num2str(numel(timeh5)), ' different time slices\n\n']);
fclose(fileID);


for i =1:numel(timeh5)

	fileID =fopen([mpath, '/',num2str(shot),'_log.txt'],'a');
	fprintf(fileID, ['Currently running ', num2str(timeh5(i)),'\n\n']);
	fclose(fileID);
	tstart = tic;

	path2inp = [prerundatapath, num2str(shot), '_', num2str(timeh5(i)),'_mi_', num2str(ionmass),'.hdf5'];
	runpath = ['/temp/markl_m/',project, '/RUNS/', studyname, '/', num2str(shot), '_', num2str(timeh5(i)), '/'];
	system(['mkdir -p ', runpath]);

% path to the input hdf5 file, that is generated in the pre run
% same for each shot/time combination, i.e. has to be done only once
   
% create balance object
	bal = Balance(runpath, shot, timeh5(i), studyname, path2inp);
	bal.setModes(m, n);
% load KiLCA configuration from prerun
	try
		bal.load_prerun(ionmass);
		bal.write_kilca();
	catch
		fileID =fopen([mpath, '/',num2str(shot),'_log.txt'],'a');
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
		system(['cp /afs/itp.tugraz.at/user/markl_m/Dokumente/plasma/code/libneo/matlab/BALANCE/blueprints/balance_conf.nml '...
			, runpath, 'balance_conf.nml']);
end

	balancenml = InputFile([runpath,'balance_conf.nml']);
	balancenml.read();
	% set options
	balancenml.BALANCENML.path2inp = path2inp;
	balancenml.BALANCENML.flag_run_time_evolution = false;
	balancenml.BALANCENML.paramscan = false;
	balancenml.BALANCENML.diagnostics_output = false;
	balancenml.BALANCENML.debug_mode = true;
	balancenml.BALANCENML.timing_mode = false; % saves time needed for certain computations
	balancenml.BALANCENML.suppression_mode = false;
	balancenml.BALANCENML.am = 2;
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
		fileID =fopen([mpath, '/',num2str(shot),'_log.txt'],'a');
		fprintf(fileID, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
		fprintf(fileID, ['Problem run for time ', num2str(timeh5(i)), ' at ', datestr(datetime('now')), ' \n']);
		fprintf(fileID, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n');
		fclose(fileID);
		continue
	end
	tend = toc(tstart);

	fileID = fopen([mpath, '/',num2str(shot),'_log.txt'], 'a');
	fprintf(fileID, '-------------------------------------------------\n');
	fprintf(fileID, ['Finished ',num2str(timeh5(i)), ' successfully at ', datestr(datetime('now')), ' \n']);
	fprintf(fileID, ['It took ', num2str(tend), ' seconds\n']);

	fprintf(fileID, '-------------------------------------------------\n\n');
	fclose(fileID);
	if i==1
		system(['mkdir -p ', '/temp/markl_m/', project,'/POST_h5/', studyname, '/', num2str(shot),'/']);
	end
	system(['mv ', runpath, 'out/', studyname, '_', num2str(shot), '_', num2str(timeh5(i)), '.hdf5 ', '/temp/markl_m/', project,'/POST_h5/', studyname, '/', num2str(shot),'/']);

end




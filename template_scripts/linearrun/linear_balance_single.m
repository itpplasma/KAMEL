%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% linear_balance_single.m
% ----------------------------------------------------------------------------
% description:
% ------------
% This script runs a linear balance run of a given time slice for a given shot.
% ----------------------------------------------------------------------------
% author: Markus Markl
% created: 23.12.2021

% TODO: change paths
libBalance = '/temp/markl_m/GITHUB/Balance/matlab/balance';
addpath(genpath(libBalance));

% specify study name, directory in runs is created with this name
studyname = 'linearrun'; 

% overall project directory
project = 'BALANCE_4CASES';

shot = 33353;
time = 2350;

% TODO: change path accordingly
prerundatapath = ['/temp/markl_m/', project, '/PRERUNDATA/', num2str(shot), '/'];

ionmass = 2;

% specify mode numbers
m = [5, 6, 7];
n = 2.*ones(size(m));

% define factors for parameter scan
fac_n = [1.0];
fac_Ti = [1.0];
fac_Te = [1.0];
fac_vz = [1.0];

% path to input (i.e. prerun) hdf5 file
path2inp = [prerundatapath, num2str(shot), '_', num2str(time), '_mi_', ...
	num2str(ionmass), '.hdf5'];

% TODO: change path
runpath = ['/temp/markl_m/', project, '/RUNS/', studyname, '/', ...
	num2str(shot), '_', num2str(time), '/'];

% create run directory
system(['mkdir -p ', runpath]);


% create balance object
bal = Balance(runpath, shot, time, studyname, path2inp);
bal.setModes(m,n);

% load prerun data
bal.load_prerun(ionmass);
bal.kil_vacuum.background.mi = ionmass;
bal.kil_flre.background.mi = ionmass;
bal.write_kilca(1); % change 1 to 0 after first run for a time slice

% check if configuration namelist file exists in run path, if not copy it
if ~exist([runpath, 'balance_conf.nml'], 'file')
	% TODO: change path
	system(['cp /temp/markl_m/GITHUB/Balance/matlab/blueprints/balance_conf.nml ',...
		runpath, 'balance_conf.nml']);
end

% read in configuration file
balancenml = InputFile([runpath, 'balance_conf.nml']);
balancenml.read();

% set options
balancenml.BALANCENML.path2inp = path2inp;
balancenml.BALANCENML.flag_run_time_evolution = false;
balancenml.BALANCENML.br_stopping = false;
balancenml.BALANCENML.paramscan = false;
balancenml.BALANCENML.diagnostics_output = false;
balancenml.BALANCENML.save_prof_time_step = 10;
balancenml.BALANCENML.debug_mode = false;
balancenml.BALANCENML.timing_mode = false;
balancenml.BALANCENML.suppression_mode = true; % if false, profiles are written out
balancenml.BALANCENML.am = ionmass; % ion mass number
balancenml.BALANCENML.flre_path = bal.kil_flre.pathofrun;
balancenml.BALANCENML.vac_path = bal.kil_vaccum.pathofrun;

% write options
bal.setOptionsNML(balancenml);

% set factors
bal.setFactors(fac_n, 'fac_n');
bal.setFactors(fac_Ti, 'fac_Ti');
bal.setFactors(fac_Te, 'fac_Te');
bal.setFactors(fac_vz, 'fac_vz');

% run
bal.run();



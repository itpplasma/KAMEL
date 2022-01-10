%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% balancerun_parscan.m
% ----------------------------------------------------------------------------
% description:
% ------------
% This script is the template for a parameter scan for a given shot/time
% combination done with the balance code.
% ----------------------------------------------------------------------------
% author: Markus Markl
% created: 07.01.2022

% TODO: change path
libBalance = '/temp/markl_m/GITHUB/Balance/matlab/balance';
addpath(genpath(libBalance))

% specify study name, directory in runs is created with this name
studyname = 'nTscan'; 

% overall project directory:
project = 'ELMsuppression_in_hydrogen';

shot = 33353;
time = 2000;
ionmass = 2;

% specify mode numbers:
m = [5, 6, 7];
n = 2 .* ones(size(m));

% TODO: change path
runpath = ['/temp/markl_m/',project, '/RUNS/', studyname, '/', num2str(shot), '_', num2str(time), '/'];
system(['mkdir -p ', runpath]);

% TODO: change path
% path to the input hdf5 file, that is generated in the pre run
% same for each shot/time combination, i.e. has to be done only once
path2inp = ['/temp/markl_m/', project,'/PRERUNDATA/', num2str(shot), '/', num2str(shot), '_', num2str(time),'.hdf5'];

   
% define factors for parameter scan
% in this example, a 2D scan in density and electron temperature is done.
fac_n = linspace(0.3,12,20);
fac_n = transpose(fac_n);
fac_Ti = [1.0]; 
fac_Te = linspace(0.3, 12,20);
fac_Te = transpose(fac_Te);
fac_vz = [1.0];

% create balance object
bal = Balance(runpath, shot, time, studyname, path2inp);
bal.setModes(m, n);

% load KiLCA configuration from prerun, set ion mass
bal.load_prerun(ionmass);

% write KiLCA configuration files
bal.write_kilca(1); % change 1 to 0 after first run for a time slice

% check if the configuration namelist file exists in run path, 
% if not, copy it from blueprint path
if ~exist([runpath, 'balance_conf.nml'], 'file')
	% TODO: change path
	system(['cp /temp/markl_m/GITHUB/Balance/matlab/blueprints/balance_conf.nml '...
		, runpath, 'balance_conf.nml']);
end

balancenml = InputFile([runpath,'balance_conf.nml']);
balancenml.read();
% set options
balancenml.BALANCENML.path2inp = path2inp;
balancenml.BALANCENML.flag_run_time_evolution = false;
balancenml.BALANCENML.br_stopping = false; % stops code, if slope is large enough

% this is the important configuration
balancenml.BALANCENML.paramscan = true;

balancenml.BALANCENML.diagnostics_output = false;
balancenml.BALANCENML.save_prof_time_step = 10;
balancenml.BALANCENML.write_formfactors = false;
balancenml.BALANCENML.debug_mode = true;
balancenml.BALANCENML.timing_mode = false; % saves time needed for certain computations
balancenml.BALANCENML.suppression_mode = false;
balancenml.BALANCENML.am = ionmass;
balancenml.BALANCENML.flre_path = bal.kil_flre.pathofrun;
balancenml.BALANCENML.vac_path = bal.kil_vacuum.pathofrun;

% write options
bal.setOptionsNML(balancenml);

% set factors
bal.setFactors(fac_n, 'fac_n')
bal.setFactors(fac_Ti, 'fac_Ti')
bal.setFactors(fac_Te, 'fac_Te')
bal.setFactors(fac_vz, 'fac_vz')

% run
bal.run();


% This script is the balance code part of the whole balance cycle.
% This configuration runs the time evolution of a specific shot/time
% created: Markus Markl, 30.03.2021

libBalance = '/temp/markl_m/GITHUB/Balance/matlab/balance';

addpath(genpath(libBalance))

studyname = 'timeevol_RMP_paper'; % specifies the directory, should describe the intent
% project directory name:
project = 'BALANCE_4CASES';

shot = 33353;
time = 2900;%shotinput.data(1,2);
ionmass = 2;
%profoffset = 0;
%equioffset = -1;

% specify mode numbers:
m = [5,6,7];
n = 2 .* ones(size(m));


runpath = ['/temp/markl_m/',project, '/RUNS/', studyname, '/', num2str(shot), '_', num2str(time), '/'];
system(['mkdir -p ', runpath]);

% path to the input hdf5 file, that is generated in the pre run
% same for each shot/time combination, i.e. has to be done only once
path2inp = ['/temp/markl_m/', project,'/PRERUNDATA/', num2str(shot), '/', num2str(shot), '_', num2str(time), '_mi_', num2str(ionmass),'_ULBLP_rho_pol.hdf5'];

   
% define factors for parameter scan
fac_n = [1.0];%[0.5; 1.0; 1.5];%[0.3; 0.5; 0.6;0.7; 0.8; 0.9; 1.0; 1.1; 1.2; ...
	%1.3;1.4;1.5;1.8;2.1;2.4;2.7;3.0];
fac_Ti = [1.0];
fac_Te = [1.0];
fac_vz = [1.0];

% create balance object
bal = Balance(runpath, shot, time, studyname, path2inp);
bal.setModes(m, n);
% load KiLCA configuration from prerun
bal.load_prerun();
% change ion mass to hydrogen
%bal.kil_vacuum.background.mi = ionmass;
%bal.kil_flre.background.mi = ionmass;
bal.write_kilca(1);


% check if the configuration namelist file exists in run path, 
% if not, copy it from blueprint path
if ~exist([runpath, 'balance_conf.nml'], 'file')
	system(['cp /temp/markl_m/GITHUB/Balance/matlab/blueprints/balance_conf.nml '...
		, runpath, 'balance_conf.nml']);
end

system(['cp /temp/markl_m/GITHUB/Balance/matlab/blueprints/balance_conf.nml '...
		, runpath, 'balance_conf.nml']);

balancenml = InputFile([runpath,'balance_conf.nml']);
balancenml.read();
% set options
balancenml.BALANCENML.path2inp = path2inp;
balancenml.BALANCENML.flag_run_time_evolution = true;
balancenml.BALANCENML.stop_time_step = 1e-8;
balancenml.BALANCENML.br_stopping = false; % stops code, if slope is large enough
balancenml.BALANCENML.faster_ramp_up = 0; % ramp up faster
%balancenml.BALANCENML.t_max_ramp_up = 0.01; % time when expt. value of ramp up is reached in sec
%balancenml.BALANCENML.timstep_min = 1e-4;

%balancenml.BALANCENML.temperature_limit = 20; % lower temperature limit in eV
balancenml.BALANCENML.antenna_max_stopping = 3;

balancenml.BALANCENML.paramscan = false;
balancenml.BALANCENML.diagnostics_output = false;
balancenml.BALANCENML.save_prof_time_step = 5;
balancenml.BALANCENML.write_formfactors = false;
balancenml.BALANCENML.debug_mode = true;
balancenml.BALANCENML.timing_mode = false; % saves time needed for certain computations
balancenml.BALANCENML.suppression_mode = false;
%balancenml.BALANCENML.am = ionmass;
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
disp('now running');
bal.run();

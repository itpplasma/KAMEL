% parameter scan, balance code
% created: Markus Markl, 30.03.2021

libBalance = '/temp/markl_m/GITHUB/Balance/matlab/balance';

addpath(genpath(libBalance))

studyname = 'gyro_current_parscan'; % specifies the directory, should describe the intent
% project directory name:
project = 'BALANCE_4CASES';

shot = 33353;
time = 2900;
ionmass = 2;

% specify mode numbers:
m = [7];
n = 2 .* ones(size(m));

runpath = ['/temp/markl_m/',project, '/RUNS/', studyname, '/', num2str(shot), '_', num2str(time), '/'];
system(['mkdir -p ', runpath]);

% path to the input hdf5 file, that is generated in the pre run
% same for each shot/time combination, i.e. has to be done only once
path2inp = ['/temp/markl_m/', project,'/PRERUNDATA/', num2str(shot), '/', num2str(shot), '_', num2str(time),'_mi_', num2str(ionmass),'.hdf5'];

   
% define factors for parameter scan
%fac_n = [0.5; 1.0; 1.5];%[0.3; 0.5; 0.6;0.7; 0.8; 0.9; 1.0; 1.1; 1.2; ...
%fac_n = linspace(0.3,12,5);
%fac_n = transpose(fac_n);
fac_n = [1.0];
	%1.3;1.4;1.5;1.8;2.1;2.4;2.7;3.0];
fac_Ti = [1.0]; %[0.5; 1.0; 1.5];
%fac_Te = [0.5; 1.0; 1.5];
%fac_Te = linspace(0.1, 10,20);
%fac_Te = transpose(fac_Te);
fac_Te = [1.0];
fac_vz = [1.0];

% create balance object
bal = Balance(runpath, shot, time, studyname, path2inp);
bal.setModes(m, n);
% load KiLCA configuration from prerun, set ion mass
bal.load_prerun(ionmass);
% write KiLCA configuration files

% Does not have any impact!!!!!
%bal.kil_flre.zones{1}.flre.flag_corr = 1;
%bal.kil_flre.background.vgalsys=-5e7;
%bal.KiLCA_vgalsys=-1e9;
%bal.KiLCA_flab = 1;
bal.write_kilca(0);

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
balancenml.BALANCENML.br_stopping = false; % stops code, if slope is large enough
balancenml.BALANCENML.paramscan = true;
balancenml.BALANCENML.diagnostics_output = false;
balancenml.BALANCENML.save_prof_time_step = 10;
balancenml.BALANCENML.write_formfactors = false;
balancenml.BALANCENML.debug_mode = true;
balancenml.BALANCENML.timing_mode = false; % saves time needed for certain computations
balancenml.BALANCENML.suppression_mode = false;
balancenml.BALANCENML.am = ionmass;
%balancenml.BALANCENML.Nstorage = 1000; % 1000 is default
balancenml.BALANCENML.flre_path = bal.kil_flre.pathofrun;
balancenml.BALANCENML.vac_path = bal.kil_vacuum.pathofrun;
balancenml.BALANCENML.gyro_current_study = 1;

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

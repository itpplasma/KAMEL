%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% template_timeevol_balance_run.m
% ----------------------------------------------------------------------------
% description:
% Template for a time evolution balance run for a single time slice of a shot.
% ----------------------------------------------------------------------------
% Markus Markl, 16.12.2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% add balance class
libBalance = '/temp/markl_m/GITHUB/Balance/matlab/balance/';
addpath(genpath(libBalance))

% define studyname, in this case time evolution abbreviated by timeevol
studyname = 'timeevol'; 

% project directory name:
% TODO: change to your needs
project = 'ELMsuppression_in_hydrogen';

shot = 38186;
time = 2000;

% ion mass can be set specificaly, default is 2, i.e. Deuterium 
%ionmass = 2;

% specify mode numbers:
m = [5,6,7];
n = 2 .* ones(size(m));

% define the path where the run is executed
% TODO: change path to your needs (i.e. /temp folder), write permissions are required
runpath = ['/temp/markl_m/', project, '/RUNS/', studyname, '/', num2str(shot), '_', num2str(time), '/'];
system(['mkdir -p ', runpath]);

% path to the input hdf5 file, that is generated in the pre run
% TODO: change accordingly, only reading, i.e. file path from /temp w/o write permission possible
path2inp = ['/temp/markl_m/', project, '/PRERUNDATA/', num2str(shot), '/', num2str(shot), '_', num2str(time),'.hdf5'];

% define factors for parameter scan, need to be set to at least one value
fac_n = [1.0];
fac_Ti = [1.0];
fac_Te = [1.0];
fac_vz = [1.0];

% for multiple factors: (i.e. column vector)
%fac_n = [0.5; 1.0; 1.5];

% create balance object
bal = Balance(runpath, shot, time, studyname, path2inp);
bal.setModes(m, n);
% load KiLCA configuration from prerun hdf5 file
bal.load_prerun();
% write KiLCA configuration files
bal.write_kilca();
% for different ion mass use:
%bal.write_kilca(ionmass);


% check if the configuration namelist file exists in run path, 
% if not, copy it from blueprint path
% TODO: change path accordingly (could work without changing it)
if ~exist([runpath, 'balance_conf.nml'], 'file')
	system(['cp /temp/markl_m/GITHUB/Balance/matlab/blueprints/balance_conf.nml '...
		, runpath, 'balance_conf.nml']);
end

% read configuration blueprint
balancenml = InputFile([runpath,'balance_conf.nml']);
balancenml.read();
%-----------------------------------------------------------------
% Define the configuration of the balance run
% ----
% set the input hdf5 file
balancenml.BALANCENML.path2inp = path2inp;

% flag for the time evolution
balancenml.BALANCENML.flag_run_time_evolution = true;

% br_stopping: if true, stops the code if the slope of the radial magnetic field
% is large enough. This is a simple detection of bifurcation. 
% Saves computation time.
balancenml.BALANCENML.br_stopping = false; % stops code, if slope is large enough

% Parameter scan flag
balancenml.BALANCENML.paramscan = false;

% Diagnostics output flag. Always set to false. Diagnostics output is not needed.
balancenml.BALANCENML.diagnostics_output = false;

% Suppression flag. This suppresses additional output, e.g. the profiles. 
% It only writes out the radial magnetic field at the resonant surface.
balancenml.BALANCENML.suppression_mode = false;

% Step intervals at which data is written out.
balancenml.BALANCENML.save_prof_time_step = 10;

% Debug flag
balancenml.BALANCENML.debug_mode = true;

% set ion mass in balance code
%balancenml.BALANCENML.am = ionmass;

% set KiLCA flre and vacuum path of run
balancenml.BALANCENML.flre_path = bal.kil_flre.pathofrun;
balancenml.BALANCENML.vac_path = bal.kil_vacuum.pathofrun;

%-----------------------------------------------------------------
% write the configuration to the file
bal.setOptionsNML(balancenml);

% set the factors
bal.setFactors(fac_n, 'fac_n')
bal.setFactors(fac_Ti, 'fac_Ti')
bal.setFactors(fac_Te, 'fac_Te')
bal.setFactors(fac_vz, 'fac_vz')

% execute the balance run
bal.run();

% done :)

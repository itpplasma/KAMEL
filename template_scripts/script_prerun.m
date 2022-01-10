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

mpath = pwd();

%Runs to make
shot = 39711; 
time = 2000;  
% profoffset is used, if names of profiles differ from time.
profoffset = 2;
equioffset = 0;
shotname = [num2str(shot), '_', num2str(time)];

%proftype = '_ULBLP_rho_pol'; %'_MICDU' or '_MMARKL_rho_pol'
proftype = '_MMARKL_rho_pol';

studyname = ['DEMO', proftype];
% TODO:change path
basepath = '/temp/markl_m/DEMONSTRATION_PROJ/';
datapath = [basepath, 'DATA/'];
path2preh5 = [basepath,'PRERUNDATA/',num2str(shot), '/', num2str(shot), '_', num2str(time),'_mi_', num2str(ion_mass) '.hdf5'];

% use copyinitial to use existing profiles to a shot/time.
%copyinitial = ['/temp/markl_m/BALANCE_4CASES/RUNS/',studyname, '/33353_2670/NRef_timeevol/profiles/'];

m = [5,6,7];
n = 2 .* ones(size(m));

%##########################################################################
% 1) PRE-RUN
%##########################################################################

runname = 'prerun';
runpath = [basepath, 'RUNS_PRE/', studyname,'/', num2str(shot), '_', num2str(time),'/',runname,'/'];

%input
gfile  = [datapath, 'EQUI/', num2str(shot),'/g',num2str(shot),'.',num2str(time+equioffset),'_EQH'];
%gfile  = [datapath, 'DATA/EQUI/', num2str(shot),'/g',num2str(shot),'.',num2str(time),'_COCOS3'];
cfile  = [datapath, 'COIL/', num2str(shot),'/',num2str(shot),'.',num2str(time),'_coil.dat'];
%TODO: change path
dapath = '/temp/markl_m/DA/ASTRA/';

filehead = [datapath, 'PROF/', num2str(shot),'/',num2str(shot),'.',num2str(time+profoffset)];
neprof = [filehead,'_ne_PED',proftype,'.dat'];
Teprof = [filehead,'_Te_PED',proftype,'.dat'];
Tiprof = [filehead,'_Ti_PED',proftype,'.dat'];
vtprof = [filehead,'_vt_PED',proftype,'.dat'];

%location of field.dat
pfile = ['/temp/markl_m/MESH3D/',num2str(shot),'/field.dat']; %will be calculated if not present
fluxdatapath = ['/temp/markl_m/FLUXDATA/',num2str(shot),'/',num2str(time),'/']; %will be calculated if not present
gpecpath = ['/temp/markl_m/GPEC/', studyname, '/', num2str(shot), '_', num2str(time), '/'];

%BALANCE CODE
bal = Balance(runpath, shot, time, [studyname, '_', runname], path2preh5);
bal.setModes(m, n);
bal.setCoil(cfile, pfile);
bal.setEqui(gfile, fluxdatapath);
bal.setTMHDCode('GPEC', gpecpath);
bal.setProfiles(neprof, Teprof, Tiprof, vtprof);
bal.setKiLCA(ion_mass);
bal.setDaEstimation(dapath);
bal.write();

system(['mkdir -p ', basepath, 'POST/', studyname, '/', shotname, '/']);
%bal.export2HDF5([basepath, 'POST/', studyname, '/', shotname, '/'], [studyname, '_', runname]);

cd(mpath);

system(['ln -sf ', path2preh5, ' ', '/temp/markl_m/inputdata_h5/', num2str(shot),'/',...
	num2str(shot), '_', num2str(time), '.hdf5']);

disp(['Finished! :), shot: ', num2str(shot), ' at time: ', num2str(time)]);

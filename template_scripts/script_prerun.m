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

ion_mass = 1; % mass is ion_mass * proton mass

% path to KiLCA interface and balance class. Needs to be changed individually.
libKiLCA = '~/Dokumente/plasma/code/libneo/matlab/KiLCA_interface/';
libBalance = '~/Dokumente/plasma/code/libneo/matlab/BALANCE/balance';

addpath(genpath(libKiLCA))
addpath(genpath(libBalance))

mpath = pwd();

%Runs to make
shot = 39711; %or 33353 2900ms
time = 2000;  %or 33120 5500ms
profoffset = 2;
equioffset = 0;
shotname = [num2str(shot), '_', num2str(time)];

%proftype = '_ULBLP_rho_pol'; %'_MICDU' or '_MMARKL_rho_pol'
proftype = '_MMARKL_rho_pol';

studyname = ['ELMSuppr_H', proftype];
basepath = '/temp/markl_m/ELMsuppression_in_hydrogen/';
datapath = [basepath, 'DATA/'];
path2preh5 = [basepath,'PRERUNDATA/',num2str(shot), '/', num2str(shot), '_', num2str(time), '.hdf5'];

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

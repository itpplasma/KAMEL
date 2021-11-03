%#################################################################
% script_balanceinit.m
%#################################################################
% Description:
%-----------------------------------------------------------------
% This script is an experimental script to test out what is
% necessary to do multiple balance configuration runs.
%#################################################################

% author: Markus Markl
% created: 10.02.2021

% path to the KiLCA interface and the Balance class
libKiLCA = '~/Dokumente/plasma/code/libneo/matlab/KiLCA_interface/';
libBalance = '~/Dokumente/plasma/code/libneo/matlab//BALANCE/balance';

addpath(genpath(libKiLCA))
addpath(genpath(libBalance))

% safe current directory
mpath = pwd();

studyname = 'Multipleshotrun_test';
system(['mkdir -p ~/BalanceProposal/', studyname, '/']); % create output directory

% specify shot and time:
% import from file with delimiter ' ' and with 1 header line
% could contain more than 1 set of shot and time number
shotinput = importdata('/temp/markl_m/AUG/SHOTS/shotnmbr.dat', ' ', 1);

shot = shotinput.data(:,1);
time = shotinput.data(:,2);

% specify mode numbers:
m = 5;
n = 2;

dapath = '/temp/ulbl_p/DA/ASTRA';

for i = 1:numel(shot)

   runpath = ['/temp/markl_m/BalanceProposal/', studyname, '/', num2str(shot(i)), '_', num2str(time(i)), '/'];
   % input files
   gfile = ['/temp/markl_m/AUG/SHOTS/', num2str(shot(i)),'/g', num2str(shot(i)), '.', num2str(time(i)), '_EQH'];
   filehead = ['/temp/markl_m/AUG/SHOTS/', num2str(shot(i)), '/', num2str(shot(i)), '.', num2str(time(i))];
   cfile = [filehead, '_coil.dat'];
   % profile data
   neprof = [filehead, '_ne_PED_ULBLP_rho_pol.dat'];
   Teprof = [filehead, '_Te_PED_ULBLP_rho_pol.dat'];
   Tiprof = [filehead, '_Ti_PED_ULBLP_rho_pol.dat'];
   vtprof = [filehead, '_vt_PED_ULBLP_rho_pol.dat'];
   %location of field.dat, will be calculated if not present
   fluxdatapath = ['/temp/markl_m/FLUXDATA/', num2str(shot(i)), '/', num2str(time(i)), '/'];
   % path to gpec
   gpecpath = ['/temp/markl_m/GPEC/balanceprop/', num2str(shot(i)), '_', num2str(time(i)), '/'];

%copy = '/temp/ulbl_p/BALANCE_2020/TimeEvolution/33133_3000/profiles/';



   % ---------- Run the Balance Code with the balance class ------------

   % create balance object
   bal = Balance(runpath, shot(i), time(i), studyname);
   % run the balance configuration
   bal.run_balance_init(m, n, cfile, gfile, fluxdatapath, gpecpath, neprof, Teprof, Tiprof, vtprof, dapath);
end

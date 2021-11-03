%#################################################################
% readshots.m
%#################################################################
% Description:
%-----------------------------------------------------------------
% Script to read the shots that are to be run through the balance
% code. They are read through the existing file names.
%#################################################################

% author: Markus Markl
% created: 08.02.2021

% path to the KiLCA interface and the Balance class
libKiLCA = '~/Dokumente/plasma/code/libneo/matlab/KiLCA_interface/';
libBalance = '~/Dokumente/plasma/code/libneo/matlab//BALANCE/balance';

addpath(genpath(libKiLCA))
addpath(genpath(libBalance))

% safe current directory
mpath = pwd();

studyname = 'readshotsviadir';
system(['mkdir -p ~/BalanceProposal/', studyname, '/']); % create output directory

% specify shot and time:
%shot = 33133;
%time = 3000;

shotspath = '/temp/markl_m/AUG/SHOTS'

checkdata = 0;
shotnumbers = [];
shots = dir(shotspath);
for i = 1:numel(shots)
  %disp(shots(i).name);
  % check if the directory name is a number
  % to do: must check at the same time if necessary files are present
  if (~isnan(str2double(shots(i).name)))
     % check if there are anyfiles
     nextdir = dir([shotspath,'/', shots(i).name]);
     %disp(shots(~cellfun('isempty', {shots.date})))
     for j = 1:numel(nextdir)
       if ~(strcmp(nextdir(j).name, '.') || strcmp(nextdir(j).name, '..'))
           disp(['yep, there is data in ', shotspath,'/', shots(i).name])
           checkdata = 1;
           break
       end
     end

     if checkdata
       disp([shots(i).name, ' not empty'])
       shotnumbers = [shotnumbers str2double(shots(i).name)];
       disp(['Added ', shots(i).name, ' to the shots']);
       checkdata = 0;
     else
       disp(['No data found in ', shotspath,'/', shots(i).name])
     end

   %  if ~isempty(nextdir.date)%(nextdir(~cellfun('isempty', {nextdir.name})))
%        disp([shots(i).name, ' not empty'])
%        shotnumbers = [shotnumbers str2double(shots(i).name)];
%        disp(['Added ', shots(i).name, ' to the shots']);
%     end

  end
  % ? check if all necessary files are present, or if any

  % run balance init

end

%disp('The list of shots is ');
%disp(shotnumbers);

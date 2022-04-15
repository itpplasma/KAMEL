
%##########################################################################
% script_resonances_from_h5index.m
%##########################################################################
% description:
%--------------------------------------------------------------------------
% Get the resonant values (zero crossings) of electric, fluid resonance
% without dTe and with dTe velocities.
%##########################################################################
%author:   Markus Markl
%created:  14.03.2022


addpath('/temp/markl_m/GITHUB/Balance/utility_scripts/matlab_utility/');

shot = 33353;

prerundatapath = ['/temp/markl_m/BALANCE_4CASES/PRERUNDATA/',num2str(shot),'/'];
outpath = ['./result/', num2str(shot), '/'];

system(['rm ', outpath, num2str(shot), '_fluid_resonance.dat']);
system(['rm ', outpath, num2str(shot), '_fluid_resonance_wdTe.dat']);
system(['rm ', outpath, num2str(shot), '_electric_resonance.dat']);
system(['rm ', outpath, num2str(shot), '_resonant_surfaces.dat']);



%create_index(shot);
create_index_from_hdf5(shot, prerundatapath);
%[time, studyname] = get_shot_times(shot);
timeh5 = get_shot_times_from_hdf5(shot);

%[fluid_resonance_h5, res_surf_h5, timeh5] = get_fluid_resonance_from_hdf5(shot, timeh5);
[fluid_resonance_h5, fluid_resonance_wdTe_h5, electric_resonance_h5, res_surf_h5, timeh5] = get_all_resonances_from_hdf5(shot, timeh5, outpath);


disp('!- done -!')

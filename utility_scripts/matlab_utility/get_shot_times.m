% get_shot_times.m
%
% This function is intended to retrieve the times and corresponding fit types
% of available shot runs from a (shot).index text file.
%

% author: Markus Markl
% created: 30.06.2021

function [time, type] = get_shot_times(shot, path);
	% time ... vector of available time slices
	% type ... study and type of fit, e.g. 4Cases_extended_MMARKL_rho_pol 

	if(nargin <2 || isempty(path))
		path = '/temp/markl_m/BALANCE_4CASES/DATA/INDEX/';
	end

	[time, type] = textread([path, num2str(shot), '.index'], '%u %s');

end

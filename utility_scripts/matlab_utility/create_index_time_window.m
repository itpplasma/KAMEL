%-------------------------------------------------------------------------
% create_index.m
%-------------------------------------------------------------------------
% This function is intended to create an index file from available run files
% of a certain shot and time that where created with plots fitted with different time
% windows. The index file contains the time window %
%-------------------------------------------------------------------------

% author: Markus Markl
% created: 30.06.2021

%-------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function time_windows = create_index_time_window(shot, time, studyname, path, spath);
	% shot ... shot number for which the index file should be created
	% studyname ... study to be considered
	% path ... path where the directory system is for which the index is 
	% should be created. Can be multiple paths. Default is the 
	% bif_time_window directory, in particular the MMARKL study.
	% spath ... path where the index file should be saved to

	runname = [num2str(shot), '_', num2str(time)];
	
	if (nargin<4 || isempty(path))
		path = '/temp/markl_m/bif_time_window/POST/';
	end

	if (nargin<4 || isempty(studynames))
		studynames = {'bif_time_window_MMARKL_rho_pol'};
	end

	if (nargin<4 || isempty(spath))
		spath = '/temp/markl_m/bif_time_window/DATA/INDEX/';
	end

	% remove existing index file
	if exist([spath, runname, '.index'])
		disp('Removing existing index file')
		delete([spath, runname, '.index']);
	end
		

	for i=1:numel(studynames)
		workingpath = [path, studynames{i}, '/', runname, '/'];
		files = dir(workingpath);
		files = files(arrayfun(@(f) f.isdir == true, files)); % remove files
		files = files(arrayfun(@(f) ~contains(f.name, '.'), files)); % remove . and .. dirs
		fnames = {files.name}; % save only directory names
		%disp(fnames)
		
		if isempty(fnames)
			continue
		end

		% handle case where only 1 file exists differently (name cell array
		% is different)
		if numel(fnames) == 1
			time_windows = fnames{1};
			%disp(times_windows)
			fileID = fopen([spath, num2str(shot), '.index'], 'a');
			fprintf(fileID, '%s %s\n', times);
			fclose(fileID);
		else
			time_windows = fnames;
			%disp(time_windows)
			fileID = fopen([spath, runname, '.index'], 'a'); 
			for j = 1:numel(time_windows)
				fprintf(fileID, '%s\n', time_windows{j});
			end
			fclose(fileID);
		end

	end

	disp(['Index file for shot ', num2str(shot), ' created']);

end

% create_index.m
%
% This function is intended to create an index file from available run files
% of a certain shot. The index file contains the time and the studyname 
% (which again contains the type of profile fit).
%

% author: Markus Markl
% created: 30.06.2021


function create_index(shot, studynames, path, spath)
	% shot ... shot number for which the index file should be created
	% studynames ... studies to be included
	% path ... path where the directory system is for which the index
	% should be created. Can be multiple paths. Default is the 
	% BALANCE_4CASES directory, in particular the MMARKL and ULBLP studies.
	% spath ... path where the index file should be saved to
	
	if (nargin<4 || isempty(path))
		path = '/temp/markl_m/BALANCE_4CASES/POST/';
	end

	if (nargin<4 || isempty(studynames))
		studynames = {'4Cases_extended_MMARKL_rho_pol', '4Cases_extended_ULBLP_rho_pol'};
	end

	if (nargin<4 || isempty(spath))
		spath = '/temp/markl_m/BALANCE_4CASES/DATA/INDEX/';
	end

	% remove existing index file
	if exist([spath, num2str(shot), '.index'])
		disp('Removing existing index file')
		delete([spath, num2str(shot), '.index']);
	end
		

	for i=1:numel(studynames)
		workingpath = [path, studynames{i}, '/'];
		files = dir(workingpath);
		files = files(arrayfun(@(f) f.isdir == true, files)); % remove files
		files = files(arrayfun(@(f) ~contains(f.name, '.'), files)); % remove . and .. dirs
		files = files(arrayfun(@(f) contains(f.name, num2str(shot)), files)); % get only directories for according shot number
		fnames = {files.name}; % save only directory names
		%disp(fnames)
		
		if isempty(fnames)
			continue
		end
		if numel(fnames) == 1
			stpair = split(fnames, '_');
			times = stpair(2);
			%disp(times{1})
			fileID = fopen([spath, num2str(shot), '.index'], 'a');
			fprintf(fileID, '%s %s\n', times{1}, studynames{i});
			fclose(fileID);
		else
			stpair = split(fnames, '_'); % split name into shot and time numbers
			%disp(stpair)
			times = stpair(:,:,2);
			%disp(times)
			fileID = fopen([spath, num2str(shot), '.index'], 'a'); 
			for j = 1:numel(times)
				fprintf(fileID, '%s %s\n', times{j}, studynames{i});
			end
			fclose(fileID);
		end

	end

	disp(['Index file for shot ', num2str(shot), ' created']);

end

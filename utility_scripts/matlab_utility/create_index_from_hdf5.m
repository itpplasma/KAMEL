% create_index.m
%
% This function is intended to create an index file from available run files
% of a certain shot. The index file contains the time and the studyname 
% (which again contains the type of profile fit).
%

% author: Markus Markl
% created: 30.06.2021


function create_index_from_hdf5(shot, path, spath)
	% shot ... shot number for which the index file should be created
	% path ... path where the directory system is for which the index
	% should be created. Can be multiple paths. Default is the 
	% BALANCE_4CASES directory, in particular the MMARKL and ULBLP studies.
	% spath ... path where the index file should be saved to
	
	if (nargin<4 || isempty(path))
		path = '/temp/markl_m/BALANCE_4CASES/PRERUNDATA/';
	end

	if (nargin<4 || isempty(spath))
		spath = '/temp/markl_m/BALANCE_4CASES/DATA/INDEX/';
	end

	% remove existing index file
	if exist([spath, num2str(shot), '.indexh5'])
		disp('Removing existing index file')
		delete([spath, num2str(shot), '.indexh5']);
	end
		

	workingpath = [path, num2str(shot),'/'];
	files = dir(workingpath);
	files = files(arrayfun(@(f) contains(f.name, num2str(shot)), files)); % get only files for according shot number
	files = files(arrayfun(@(f) contains(f.name, 'hdf5'), files)); % get only hdf5 files 
	fnames = {files.name}; % save only file names
		
	if isempty(fnames)
		error([path, ' is empty'])
	end
	if numel(fnames) == 1
		stpair = split(fnames, '_');
		times = stpair(2);
		%disp(times{1})
		fileID = fopen([spath, num2str(shot), '.indexh5'], 'a');
		fprintf(fileID, '%s\n', times{1});
		fclose(fileID);
	else
		stpair = split(fnames, '_'); % split name into shot and time numbers
		%disp(stpair)
		times = stpair(:,:,2);
		%disp(times)
		fileID = fopen([spath, num2str(shot), '.indexh5'], 'a'); 
		for j = 1:numel(times)
			fprintf(fileID, '%s\n', times{j});
		end
		fclose(fileID);
	end

	disp(['Index file for shot ', num2str(shot), ' created']);

end




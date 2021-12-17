%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% write_resonance_to_file.m
%#######################################################################
% description:
% ----------------------------------------------------------------------
% This function is used to write resonance values to a file.
% For multiple resonance values for one time slice, it writes multiple
% lines.
%#######################################################################
% author: Markus Markl
% created: 30.11.2021

function write_resonance_to_file(time, value, fileID)

	formatstring = '%u %u\n';
	for i=1:numel(value)
		fprintf(fileID, formatstring, time, value(i));
	end

	fclose(fileID);

end

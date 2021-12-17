
function write_to_file(time, value, fileID)

	formatstring = '%u';
	for i=1:numel(value)
		formatstring = strcat(formatstring, ' %u');
	end
	formatstring = strcat(formatstring, '\n');
	%disp(formatstring)
	fprintf(fileID, formatstring, time, value);
	fclose(fileID);

end

function generate_parabolic_profiles(path, q0, n0, Te0, Ti0, Vz0, Er0, Vth0, rmin, rmax, num, a, const)

	if (nargin < 14 || isempty(const))
		const = '';	
	end

	r = linspace(rmin, rmax, num);
	
	n  = n0  .* (1- r.^2/a^2);
	Te = Te0 .* (1- r.^2/a^2);
	Ti = Ti0 .* (1- r.^2/a^2);
	Vz = Vz0 .* (1- r.^2/a^2);
	Er = Er0 .* (1- r.^2/a^2);
	Vth = Vth0 .* (1- r.^2/a^2);
	q  = -(1.05 + q0 .* r.^2/a^2);

	if ~strcmp(const,'')
		eval([const '=' const '0 .* ones(size(r));'])
	end

	fileID = fopen([path, '/q.dat'], 'w');
	fprintf(fileID, '%6.4f %6.4f\n', [r; q]);
	fclose(fileID);

	fileID = fopen([path, '/n.dat'], 'w');
	fprintf(fileID, '%6.4f %6.4f\n', [r; n]);
	fclose(fileID);

	fileID = fopen([path, '/Te.dat'], 'w');
	fprintf(fileID, '%6.4f %6.4f\n', [r; Te]);
	fclose(fileID);

	fileID = fopen([path, '/Ti.dat'], 'w');
	fprintf(fileID, '%6.4f %6.4f\n', [r; Ti]);
	fclose(fileID);

	fileID = fopen([path, '/Vz.dat'], 'w');
	fprintf(fileID, '%6.4f %6.4f\n', [r; Vz]);
	fclose(fileID);

	fileID = fopen([path, '/Er.dat'], 'w');
	fprintf(fileID, '%6.4f %6.4f\n', [r; Er]);
	fclose(fileID);

	fileID = fopen([path, '/Vth.dat'], 'w');
	fprintf(fileID, '%6.4f %6.4f\n', [r; Vth]);
	fclose(fileID);
end

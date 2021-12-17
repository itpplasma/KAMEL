%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find_resonance.m
%##############################################################################
% Description:
%------------------------------------------------------------------------------
% Finds the resonance (zero crossing) of a given profile and returns the 
% radial value.
%##############################################################################
% author: Markus Markl
% created: 30.11.2021


function [r_resonance, more] = find_resonance(prof, r, r_res)

	zero_crossing_index = @(v) find(diff(sign(v)));

	pos = zero_crossing_index(prof);
	%r_resonance = r(zero_crossing_index(prof));
	if numel(pos) > 1 
		more = numel(pos);
		r_resonance = zeros(numel(pos), 1);
		for i=1:numel(pos)
			disp(['pos = ', num2str(pos(i))]);
			r_resonance(i) = interp1(prof(pos(i):pos(i)+1),r(pos(i):pos(i)+1),0);
		end
	else
			more = 1;
			r_resonance = interp1(prof(pos:pos+1),r(pos:pos+1),0);
	end


%	if numel(r_resonance) > 1
%		more = numel(r_resonance);
%		[a, ind1] = min(abs(r_res - r_resonance), [], 2);
%		[b, ind2] = min(a, [], 1);
%		r_resonance = r_resonance(ind2);
%	else
%		more = 1;
%	end




end

function applyPropertiesToHandle(h, pr)
	fields = fieldnames(pr);
	for ii = 1:length(fields)
		if isprop(h, fields{ii})
			set(h, fields{ii}, pr.(fields{ii}));
		else
			warning('matlabHelpers:applyPropertiesToHandle:noSuchProperty', 'Property ''%s'' does not exist. Ignoring ...', fields{ii});
		end
	end
end

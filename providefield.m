function S = providefield(S, varargin)
% function S = providefield(S, name, value)
%
% If 'S.name' is nonexistent, it is set to 'value'
%
% Works with a variable number of names, i.e.,
%
%     providefield(S, name1, name2, name3, value)
%
% is also possible.

names = varargin(1:nargin-2);
value = varargin{end};

try
	v = getfield(S, names{:});
catch ME
	if (strcmp(ME.identifier,'MATLAB:nonExistentField'))
		S = setfield(S, names{:}, value);
	else
		rethrow(ME);
	end
end

end

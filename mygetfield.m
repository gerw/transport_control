function val = mygetfield( S, name, default)
% function val = mygetfield( S, name, default)
%
% Like
%
%   val = getfield( S, name )
%
% but you can provide a default value.

if isfield( S, name )
	val = getfield( S, name );
else
	val = default;
end

end

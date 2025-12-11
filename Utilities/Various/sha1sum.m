function sum = sha1sum(varargin)

    if mrstPlatform('matlab')

        % Use persistent Java MessageDigest object to avoid repeated initialization
        % persistent md

        % if isempty(md)
        %     md = java.security.MessageDigest.getInstance('SHA-1');
        % else
        %     md.reset(); % Reset for reuse
        % end

        md = java.security.MessageDigest.getInstance('SHA-1');

        % Serialize the input arguments to byte arrays
        bytes = cellfun(@(x) getByteStreamFromArray(x), varargin, 'UniformOutput', false);

        % Sort to get consistent hash, i.e. sha1sum('a', 'b') ==
        % sha1sum('b', 'a')
        flat = sort([bytes{:}]);

        % Hash
        md.update(flat);
        hash = typecast(md.digest(), 'uint8');

        % Convert to lowercase hex string
        sum = lower(reshape(dec2hex(hash), 1, []));

    else

        error('Only supported on MATLAB platform');

    end

end

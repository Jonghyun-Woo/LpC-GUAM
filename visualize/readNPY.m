function data = readNPY(filename)
% Minimal NumPy .npy reader (no toolbox dependency). Returns a MATLAB array
% indexed like the original numpy array. Supports v1.0/2.0 headers, numeric
% dtypes, C- and Fortran-order.
    fid = fopen(filename, 'r');
    if fid == -1
        error('readNPY:open', 'Could not open %s', filename);
    end
    cleaner = onCleanup(@() fclose(fid));

    magic = fread(fid, 6, 'uint8=>char')';
    if ~strcmp(magic, sprintf('%cNUMPY', 147))
        error('readNPY:magic', '%s is not a valid .npy file', filename);
    end
    verMajor = fread(fid, 1, 'uint8');
    fread(fid, 1, 'uint8');
    if verMajor >= 2
        headerLen = fread(fid, 1, 'uint32=>double');
    else
        headerLen = fread(fid, 1, 'uint16=>double');
    end
    header = fread(fid, headerLen, 'uint8=>char')';

    descr        = regexp(header, '''descr''\s*:\s*''([^'']+)''', 'tokens', 'once');
    fortranTok   = regexp(header, '''fortran_order''\s*:\s*(\w+)', 'tokens', 'once');
    shapeTok     = regexp(header, '''shape''\s*:\s*\(([^)]*)\)', 'tokens', 'once');
    descr        = descr{1};
    fortranOrder = strcmpi(fortranTok{1}, 'True');
    shape        = str2num(['[' shapeTok{1} ']']); %#ok<ST2NM>
    if isempty(shape), shape = 1; end

    endian = 'l';
    if descr(1) == '>', endian = 'b'; end
    typeMap = struct('f4','single','f8','double', ...
                     'i1','int8','i2','int16','i4','int32','i8','int64', ...
                     'u1','uint8','u2','uint16','u4','uint32','u8','uint64');
    key = descr(regexp(descr, '[a-zA-Z]\d') + (0:1));
    if ~isfield(typeMap, key)
        error('readNPY:dtype', 'Unsupported dtype %s', descr);
    end
    mtype = typeMap.(key);

    raw = fread(fid, prod(shape), ['*' mtype], 0, endian);

    if numel(shape) > 1
        if fortranOrder
            data = reshape(raw, shape);
        else
            data = permute(reshape(raw, fliplr(shape)), numel(shape):-1:1);
        end
    else
        data = raw;
    end
end

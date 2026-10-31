function F = load_velocity_npy(filename)
% Read the float64 C-order cubic NPY files exported by this experiment.
fid=fopen(filename,'rb','ieee-le');
if fid<0, error('Cannot open velocity field'); end
cleanup=onCleanup(@() fclose(fid));
magic=fread(fid,6,'*uint8')';
if ~isequal(magic,uint8([147,78,85,77,80,89]))
    error('Invalid NPY magic');
end
version=fread(fid,2,'*uint8');
if version(1)==1
    header_length=fread(fid,1,'uint16');
elseif version(1)==2 || version(1)==3
    header_length=fread(fid,1,'uint32');
else
    error('Unsupported NPY version');
end
header=char(fread(fid,header_length,'*uint8')');
if isempty(strfind(header,'''<f8''')) || isempty(strfind(header,'False'))
    error('Expected little-endian float64 C-order data');
end
shape=regexp(header,'\(\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*,?\s*\)','tokens','once');
if isempty(shape), error('Expected a three-dimensional array'); end
dims=cellfun(@str2double,shape);
if ~all(dims==dims(1)), error('Expected a cubic field'); end
values=fread(fid,prod(dims),'*double');
if numel(values)~=prod(dims), error('Truncated velocity field'); end
F=permute(reshape(values,dims),[3,2,1]);
end

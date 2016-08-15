function atf = LoadAtf(atfFile)
% atf = LoadAtf(atfFile)
% Open an Axon Text File and read its data into a useful structure
% Inputs:
%   -atfFile:  string specifying file name
% Outputs:
%   -atf: structure with fields
%      version: ATF version number
%      header: cell array with header info
%      data: structure with named fields corresponding to data columns in ATF
%      units: structure with strings specifying the units of each data column

% open the file
fid = fopen(atfFile);
if fid == -1
  error('File %s does not exist', atfFile)
end

try
  % get the data from the Axon Text File
  [version, header, labels, units, dataMat] = getAtfData(fid);
catch err
  fclose(fid);
  rethrow(err);
end
fclose(fid);

% organize the data into a useful structure
atf = organizeAtfData(version, header, labels, units, dataMat);
end

function [version, header, labels, units, dataMat] = getAtfData(fid)
% ensure that it's a valid ATF file
line = fgetl(fid);
version = sscanf(line, 'ATF\t%f', 1);
if ~any(version)
  error('Corrupted file or not an Axon Text File')
end

% get the header info
line = fgetl(fid);
headerSize = sscanf(line, '%d %d');
numHeader = headerSize(1); numDataCols = headerSize(2);
header = cell(1, numHeader);
for n = 1:numHeader
  header{n} = fgetl(fid);
end

% get information about what the traces are
labelsLine = fgetl(fid);
quoteInds = strfind(labelsLine, '"');
labels = cell(1, numDataCols);
units  = cell(1, numDataCols);
for n = 1:numDataCols
  i1 = quoteInds(2 * n - 1) + 1;
  i2 = quoteInds(2 * n) - 1;
  label = labelsLine(i1:i2);
  i1 = strfind(label, '(');
  if any(i1)
    units{n} = strtrim( label(i1+1:end-1) );
    label = strtrim( label(1:i1-2) );
  else
    units{n} = '';
  end
  splitLabel = strsplit(label, {' ', '-'});
  labels{n} = splitLabel{end};
end

% read in the data
dataMat = fscanf(fid, '%f', [numDataCols, inf]);
end


function atf = organizeAtfData(version, header, labels, units, dataMat)
% organize the data into a useful structure

unitsStruct = struct();
data = struct();
for n=1:size(dataMat, 1)
  unitsStruct.(labels{n}) = units{n};
  data.(labels{n}) = dataMat(n,:);
end

atf = struct(...
  'version', version, ...
  'units', unitsStruct, ...
  'data', data ...
);
atf.('header') = header;
end
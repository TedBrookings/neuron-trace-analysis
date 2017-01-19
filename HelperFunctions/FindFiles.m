%fileNames = FindFiles( searchDirs, fileType )
% Search directory (or cell array of directories) for files of given type
function fileNames = FindFiles( searchDirs, fileType, varargin )
  parser = inputParser();
  parser.addParameter( 'recursive', true )

  parser.parse( varargin{:} )
  options = parser.Results;
  
  % ensure that searchDirs is a cell array of strings
  if ischar( searchDirs )
    searchDirs = { searchDirs };
  end
  % if fileType is not empty, ensure it starts with a '.'
  if ~isempty( fileType ) && fileType(1) ~= '.'
    fileType = ['.', fileType];
  end
  % search
  fileNames = {};
  while ~isempty( searchDirs )
    dirName = searchDirs{1};
    searchDirs(1) = [];
    dirList = dir( dirName );
    for n = 1:numel( dirList )
      listing = dirList(n);
      if listing.isdir
        if ~options.recursive || ismember( listing.name, {'.', '..'} )
          continue
        end
        searchDirs = [searchDirs, fullfile( dirName, listing.name )]; %#ok<AGROW>
      else
        [~, ~, fType] = fileparts( listing.name );
        if strcmp( fType, fileType )
          fileNames = [fileNames, fullfile( dirName, listing.name )]; %#ok<AGROW>
        end
      end
    end
  end
end

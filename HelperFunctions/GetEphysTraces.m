% [voltageTrace, currentTrace, time] = GetEphysTraces( ephysData )
% return traces from loaded .abf file by units
function varargout = GetEphysTraces( ephysData, varargin )
  parser = inputParser();
  parser.addParameter( 'voltageUnits', {'mV'} )
  parser.addParameter( 'currentUnits', {'pA', 'nA'} )
  
  parser.parse( varargin{:} )
  options = parser.Results;
  
  if ischar( ephysData )
    fileName = ephysData;
    ephysData = LoadAbf( fileName );
  end
  
  if nargout >= 3
    varargout{3} = ephysData.time;
  end
  if nargout >= 2
    varargout{2} = getTrace( ephysData, options.currentUnits );
  end
  varargout{1} = getTrace( ephysData, options.voltageUnits );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trace = getTrace( ephysData, wantedUnits )
  fieldName = getFieldName( ephysData, wantedUnits );
  trace = ephysData.data.(fieldName);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fieldName = getFieldName( ephysData, wantedUnits )
  checkUnits = lower( wantedUnits );
  fNames = fieldnames( ephysData.units );
  for n = 1:numel( fNames )
    fieldName = fNames{n};
    if ismember( lower( ephysData.units.(fieldName) ), checkUnits )
      return
    end
  end
  error( 'Unable to find field with units: %s', strjoin( wantedUnits ) )
end
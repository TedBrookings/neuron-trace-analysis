function minis = GetMinis( dT, v, varargin )
  parser = inputParser();
  parser.KeepUnmatched = true;
  parser.addParameter( 'bracketWidth', 25.0 ) % was 50, change with slowTimeFactor
  parser.addParameter( 'slowTimeFactor', 25.0 ) % was 200, change with bracketWidth
  parser.addParameter( 'minCutoffDiff', 0.001 )
  parser.addParameter( 'minSpikeHeight', 0.0 )
  parser.addParameter( 'minSpikeAspect', 0.0 )
  parser.addParameter( 'pFalseSpike', 0.05 )
  parser.addParameter( 'recursive', true )
  parser.addParameter( 'outlierFraction', 0.0 )
  parser.addParameter( 'noiseCheckQuantile', 0.5 ) % was 0.25
  parser.addParameter( 'minSpikeWidth', 0.1 ) % <- note this stupid conservative due to terrible width measurement in small minis
  parser.addParameter( 'useDerivatives', false )

  parser.parse( varargin{:} )
  options = parser.Results;
  unmatchedFields = fieldnames( parser.Unmatched );
  for n = 1:numel( unmatchedFields )
    options.(unmatchedFields{n}) = parser.Unmatched.(unmatchedFields{n});
  end
  
 
  if options.useDerivatives
    % this is here as a place-holder so that derivatives-based detection
    % can have different options. At the moment it's unnecessary since I
    % haven't tuned derivatives-based detection
    miniOptions = struct( ...
      'useDerivatives', true, ...
      'discountNegativeDerivatives', true ...
    );
  else
    miniOptions = struct( ...
      'useDerivatives', false ...
    );
  end
  fNames = fieldnames( miniOptions );
  for n = 1:numel( fNames )
    fName = fNames{n};
    if isempty( options.(fName) )
      options.(fName) = miniOptions.(fName);
    end
  end


  
  %if isfield( options, 'findMinis' ) && options.findMinis
  %  options = rmfield( options, 'findMinis' );
  %end
  options.findMinis = true;
  minis = GetSpikes( dT, v, options );
end
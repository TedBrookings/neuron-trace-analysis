function minis = GetMinis( dT, v, varargin )
  parser = inputParser();
  parser.KeepUnmatched = true;
  parser.addParameter( 'bracketWidth', [] )
  parser.addParameter( 'minCutoffDiff', 0.001 )
  parser.addParameter( 'minSpikeHeight', 0.0 )
  parser.addParameter( 'minSpikeAspect', 0.0 )
  parser.addParameter( 'pFalseSpike', [] )
  parser.addParameter( 'discountNegativeDeriv', true )
  parser.addParameter( 'recursive', true )
  parser.addParameter( 'outlierFraction', 0.0 )
  parser.addParameter( 'noiseCheckQuantile', [] )
  parser.addParameter( 'slowTimeFactor', 20.0 )
  parser.addParameter( 'useDerivatives', false )

  parser.parse( varargin{:} )
  options = parser.Results;
  unmatchedFields = fieldnames( parser.Unmatched );
  for n = 1:numel( unmatchedFields )
    options.(unmatchedFields{n}) = parser.Unmatched.(unmatchedFields{n});
  end
  
 
  if options.useDerivatives
    miniOptions = struct( ...
      'bracketWidth', 100.0, ...
      'pFalseSpike', 0.25, ...
      'noiseCheckQuantile', 0.55 ...
    );
  else
    miniOptions = struct( ...
      'bracketWidth', 25.0, ...
      'pFalseSpike', 0.05, ...
      'noiseCheckQuantile', 0.1 ...
    );
  end
  fNames = fieldnames( miniOptions );
  for n = 1:numel( fNames )
    fName = fNames{n};
    if isempty( options.(fName) )
      options.(fName) = miniOptions.(fName);
    end
  end


  
  if isfield( options, 'findMinis' ) && options.findMinis
    options = rmfield( options, 'findMinis' );
  end
  minis = GetSpikes( dT, v, options );
end
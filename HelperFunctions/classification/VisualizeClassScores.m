% varargout = VisualizeClassScores( classScores, labels, varargin )
% Make dot-density plot visualizing "class scores" for multidimensional
% classifier (class scores are nominally 0-1 numbers specifying how close
% a point is to the centroid of each class). As a measure of fidelity, Z'
% is superposed over the data for each class.
% INPUTS
%  classScores: numData x numClasses matrix of class scores
%  labels: numData x 1 list of labels of class identity (cell of strings or
%          array of numbers)
% OUTPUTS:
%  fig: handle to figure
%  zPrime: 1 x numClasses list of Z'
function varargout = VisualizeClassScores( classScores, labels, varargin )
  parser = inputParser();
  parser.addParameter( 'title', 'Class scores' )
  parser.addParameter( 'figBackgroundColor', [0 0 0] )
  parser.addParameter( 'colorMap', 'hsv' )
  parser.addParameter( 'xSpread', 0.15 )
  parser.addParameter( 'windowStyle', 'docked' )
  okayAnnotation = @(a) ismember( a, {'% correct', 'zPrime'} );
  parser.addParameter( 'annotation', '% correct', okayAnnotation )
  parser.addParameter( 'annotationFontSize', 12 )
  
  parser.parse( varargin{:} )
  options = parser.Results;


  [uniqueLabels, ~, classInd] = unique( labels );
  numClasses = numel( uniqueLabels );
  colors = num2cell( feval( options.colorMap, numClasses ), 2 );
  
  x = classInd + options.xSpread .* randn( size( classInd ) );
  for k = numel( labels ):-1:1
    y(k) = classScores(k,classInd(k));
  end
  
  fig = NamedFigure( options.title, 'WindowStyle', options.windowStyle, ...
                     'Color', options.figBackgroundColor );
  clf( fig )
  ax = axes( 'Parent', fig, 'LooseInset', [0 0 0 0], ...
             'OuterPosition', [0 0 1 1], 'Color', fig.Color, ...
             'XColor', 1.0 - fig.Color, 'YColor', 1.0 - fig.Color );
  hold( ax, 'on' )
  for n = 1:numClasses
    plot( ax, x(classInd==n), y(classInd==n), 'o', 'Color', colors{n} )
  end
  axis( ax, 'tight' )
  yRange = ylim( ax );
  for n = numClasses:-1:1
    switch options.annotation
      case '% correct'
        pCorrect = 100 * sum( classScores(classInd==n,n) > 0.5 ) ...
                 / sum( classInd == n );
        aStr = sprintf( '%.1f%%', pCorrect );
      case 'zPrime'
        zPrime(1,n) = getZPrime( classScores(classInd==n,n), ...
                                 classScores(classInd~=n,n) );
        aStr = sprintf( 'Z'' = %.2f', zPrime(n) );
    end
    text( ax, n, yRange(2), aStr, ...
          'Color', 1.0 - options.figBackgroundColor, ...
          'FontName', 'Arial', 'FontSize', options.annotationFontSize, ...
          'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom' )
  end
  plot( ax, [0.5, numClasses + 0.5], [0.5 0.5], 'w--' )

  title( ax, options.title )
  ax.XTick = 1:numel( uniqueLabels );
  ax.XTickLabel = uniqueLabels;
  ylabel( ax, 'Score of correct class' )
  
  if nargout > 0
    varargout = { fig, zPrime };
  else
    varargout = {};
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function zPrime = getZPrime( vals1, vals2 )
  zPrime = 1 - 3 * (std( vals1 ) + std( vals2 )) ...
                 / abs( mean( vals1 ) - mean( vals2 ) );
end
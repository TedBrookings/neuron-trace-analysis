% VisualizeClassification( projectedProperties, classLabels, varargin )
% Make scatter plot visualizing projectedProperties.
% -If projected properties is 2D or more, plot 1st and 2nd dimension as
%  scatter plot
% -If projected properties is 1D, plot 1st dimension as x, and make y
% cumulative by point number from 0 to 1
% INPUTS
%  projectedProperties: numData x numProperties matrix of
%                    dimensionality-reduced properties output by classifier
%  classLabels: numData x 1 list of labels of class membership (cell array
%               of strings or vector of numbers)
% OUTPUTS
%  fig: handle to figure
function VisualizeClassification( projectedProperties, classLabels, varargin )
  parser = inputParser();
  parser.addParameter( 'title', 'Visualize Classification' )
  parser.addParameter( 'colorMap', 'hsv' )
  parser.addParameter( 'figBackgroundColor', [0 0 0] )
  parser.addParameter( 'windowStyle', 'docked' )
  parser.parse( varargin{:} )
  options = parser.Results;
  
  allMarkers = {'o' '+' '*' '.' 'x' 's' 'd' '^' 'v' '>' '<' 'p' 'h'};
  badMarkers = {'.'};
  allMarkers = setdiff( allMarkers, badMarkers );
  numMarkers = numel( allMarkers );
  [allLabels, ~, labelInd] = unique( classLabels );
  numClasses = numel( allLabels );
  numProperties = size( projectedProperties, 2 );
  colors = num2cell( feval( options.colorMap, numClasses ), 2 );

  fig = NamedFigure( options.title, 'Color', options.figBackgroundColor, ...
                     'WindowStyle', options.windowStyle );
  clf( fig )
  ax = axes( 'Parent', fig, 'LooseInset', [0 0 0 0], ...
             'OuterPosition', [0 0 1 1], 'Color', fig.Color, ...
             'XColor', 1.0 - fig.Color, 'YColor', 1.0 - fig.Color );
  hold( ax, 'on' );
  
  for n = 1:numClasses
    marker = allMarkers{mod( n-1, numMarkers ) + 1};
    color = colors{n};
    classMembers = labelInd == n;
    x = projectedProperties(classMembers,1);
    if numProperties >= 2
      y = projectedProperties(classMembers,2);
    else
      x = sort( x );
      yStop = n * 0.5; yStart = yStop - 0.5;
      y = reshape( linspace( yStart, yStop, numel( x ) ), size( x ) );
    end
    plot( ax, x, y, marker, 'Color', color, ...
          'DisplayName', allLabels{n} )
  end
  legend( ax, {}, 'Location', 'Best', 'Color', ax.Color, ...
          'TextColor', 1.0 - ax.Color )
  axis( ax, 'tight' )
  
  if nargout > 0
    varargout = { fig };
  else
    varargout = {};
  end
end
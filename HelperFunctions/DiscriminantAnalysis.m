%[classifyFunc, projectedProperties, classLabels] ...
%    = DiscriminantAnalysis( propertiesMat, labels, varargin )
% INPUTS:
%   propertiesMat: numData x numProperties matrix of doubles,
%                  all properties must be finite
%   labels: numData x 1 array of labels, either strings or numbers,
%           labeling the "class" of each data point
% OPTIONS:
%   kernel:  function that allows for nonlinear discriminants by redefining
%            the dot-product. defaults to Euclidean (linear discriminant)
%   numFinalDimensions: Maximum dimensions of final data. Actual final
%   dimension will be the smallest of: this option,
%                                      initial dimension of data
%                                      number of classes
% OUTPUTS:
%   classifyFunc( properties ) returns as first argument,
%                                a lower-dimensional projection of the
%                                properties that yields the best separation
%                              returns as second argument,
%                                the best estimate of the class label
%   projectedProperties   a lower-dimensional projection of the training
%                         data
%   classLabels           estimated class labels of the training data
function varargout = DiscriminantAnalysis( propertiesMat, labels, varargin )
  parser = inputParser();
  parser.addParameter( 'kernel', @(a,b) a * b' )
  parser.addParameter( 'numFinalDimensions', Inf )
  parser.addParameter( 'floatTol', 1e-12 )
  parser.addParameter( 'regularization', 1e-9 )
  parser.addParameter( 'minConditionNumber', 0.25 )
  parser.addParameter( 'preprocess', true )
  
  parser.parse( varargin{:} )
  options = parser.Results;
  
  assert( nargout <= 3, 'Too many output arguments' )
  varargout = cell( 1, nargout );
  
  % preprocess data if requested
  [propertiesMat, options] = preprocessProperties( propertiesMat, options );
  
  % find the classes (unique labels),
  %      the rows that are members of each class
  % and update the numFinalDimensions to be consistent with # of classes
  [classes, classMembers, options] = getClassInfo( labels, propertiesMat, ...
                                                   options );
  
  % get between-class and within-class kernel dot-products
  [allPairwiseDotMat, betweenClassDotMat, withinClassDotMat] = ...
    getBetweenAndWithinDotProducts( propertiesMat, classMembers, options );

  % solve for eigenvectors of discriminant analysis problem
  %  betweenDots*betweenDots' * vec = lambda * withinDots*withinDots' * vec
  eigenVecs = solveEigenProblem( betweenClassDotMat, withinClassDotMat, ...
                                 options );
  
  % return 1) function handle that does projection and classification
  %        2) projected training data and classification of training data
  [varargout{:}] = getClassifyFunc( eigenVecs, propertiesMat, ...
                                    classes, classMembers, ...
                                    allPairwiseDotMat, options );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preprocess data if requested
function [propertiesMat, options] = preprocessProperties( propertiesMat, ...
                                                          options )
  if ~options.preprocess
    options.shift = []; options.scale = [];
    return
  end
  options.shift = mean( propertiesMat, 1 );
  propertiesMat = bsxfun( @minus, propertiesMat, options.shift );
  options.scale = std( propertiesMat, 1 );
  propertiesMat = bsxfun( @rdivide, propertiesMat, options.scale );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find the classes (unique labels),
%      the rows that are members of each class
% and update the numFinalDimensions to be consistent with # of classes
function [classes, classMembers, options] = getClassInfo( labels, ...
                                                   propertiesMat, options )
  if ~iscolumn( labels ), labels = labels'; end
  [classes, ~, classIndex] = unique( labels );
  numClasses = numel( classes );
  for n = numClasses:-1:1
    classMembers{n} = find( classIndex == n );
  end
  numInitialDimensions = size( propertiesMat, 2 );
  options.numFinalDimensions = min( [options.numFinalDimensions, ...
                                     numClasses - 1, ...
                                     numInitialDimensions] );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get all pairwise, between-class, and within-class kernel dot-products
% allPairwiseDot: numData x numData matrix of dot products
% betweenClassDotMat: numData x numClasses matrix of dot products

function [allPairwiseDotMat, betweenClassDotMat, withinClassDotMat] = ...
    getBetweenAndWithinDotProducts( propertiesMat, classMembers, options )
  allPairwiseDotMat = getAllPairwiseDotMat( propertiesMat, options );
  % get betweenClassDotMat: numData x numClasses matrix of dot products
  % and withinClassDotMat: numData x numData matrix of dot products
  [betweenClassDotMat, withinClassDotMat] ...
    = getClassDotMatrices( allPairwiseDotMat, classMembers );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get allPairwiseDot: numData x numData matrix of dot products
function allPairwiseDotMat = getAllPairwiseDotMat( propertiesMat, options )
  numData = size( propertiesMat, 1 );
  kernel = options.kernel;
  try
    % try vectorized solution first
    allPairwiseDotMat = kernel( propertiesMat, propertiesMat );
    assert( isequal( size( allPairwiseDotMat ), [numData numData] ) )
  catch
    % have to use nested for loops
    for ind1 = numData:-1:1
      vec1 = propertiesMat(ind1,:);
      allPairwiseDotMat(ind1,ind1) = kernel( vec1, vec1 );
      for ind2 = numData:-1:ind1+1
        dot = kernel( vec1, propertiesMat(ind2,:) );
        allPairwiseDotMat(ind1,ind2) = dot;
        allPairwiseDotMat(ind2,ind1) = dot;
      end
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get betweenClassDotMat: numData x numClasses matrix of dot products
% and withinClassDotMat: numData x numData matrix of dot products
function [betweenClassDotMat, withinClassDotMat] ...
    = getClassDotMatrices( allPairwiseDotMat, classMembers )
  % for each data point, get average dot to all data points
  allMeanDots = mean( allPairwiseDotMat, 2 );
  % for each data point, get average dot to members of each class  
  numClasses = numel( classMembers ); numData = numel( allMeanDots );
  % allocate space for withinClassDotMat
  withinClassDotMat(numData, numData) = 0;
  % calculate within and between class dot matrices
  for col = numClasses:-1:1
    dataIndices = classMembers{col};
    coef = sqrt( numel( dataIndices ) );
    % for each data point, get average dot to all members of this class
    classDots = mean( allPairwiseDotMat(:,dataIndices), 2 );
    % between-class dot product is proportional to difference between dot
    % to class members and dot to all points
    betweenClassDotMat(:,col) = coef .* (classDots - allMeanDots);
    % for the members of each class, withinClassDotMat is difference
    % between all pairwise dots and mean of dots to members of their class
    withinClassDotMat(:,dataIndices) = ...
      bsxfun( @minus, allPairwiseDotMat(:,dataIndices), classDots );
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve for eigenvectors of discriminant analysis problem
%  betweenDots*betweenDots' * vec = lambda * withinDots*withinDots' * vec
function varargout = solveEigenProblem( betweenClassDotMat, ...
                                        withinClassDotMat, options )
  withinClassDotMat = makeWellConditioned( withinClassDotMat, options );
  
  % [U,V,X,C,S] = gsvd( A, B );
  % A = U*C*X'
  % B = V*S*X'
  % C and S non-negative diagonal so that C'*C + S'*S = I
  % U and V are unitary
  [~, ~, X, C, S] = gsvd( betweenClassDotMat', withinClassDotMat' );
  % X * C' * C * X' * vec = lambda * X * S' * S * X' * vec
  % vec = X / ( X' * X );
  %vec = real( X * pinv( X' * X ) );
  vec = X;
  lC = diag( C' * C ); lS = diag( S' * S );
  lambda = lC ./ lS;
  [lambda, sortInd] = sort( lambda, 'descend' );
  lambda = lambda';
  vec = vec(:,sortInd);
  
  numDims = min( size( vec, 2 ), options.numFinalDimensions );
  vec(:,numDims+1:end) = []; lambda(numDims+1:end) = [];
  
  for col = 1:size( vec, 2 )
    vecNorm = norm( vec(:,col) );
    vec(:,col) = vec(:,col) ./ vecNorm;
  end
  
  varargout = {vec, lambda};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mat = makeWellConditioned( mat, options )
  initialRCond = rcond( mat )
  if rcond( mat ) < options.minConditionNumber
    reg = options.regularization .* eye( size( mat ) );
    mat = mat + reg;
    while rcond( mat ) < options.minConditionNumber
      reg = reg .* 10;
      mat = mat + reg;
    end
  end
  finalRCond = rcond( mat )
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% return 1) function handle that does projection and classification
%        2) projected training data and classification of training data
function varargout = ...
    getClassifyFunc( eigenVecs, propertiesMat, classes, classMembers, ...
                     allPairwiseDotMat, options )
  projectedProperties = allPairwiseDotMat * eigenVecs;
  
  numClasses = numel( classes );
  for n = numClasses:-1:1
    classInds = classMembers{n};
    classMeanProps(n,:) = mean( projectedProperties(classInds,:), 1 );
  end
  
  kernel = options.kernel; shift = options.shift; scale = options.scale;
  classifyFunc = ...
    @(properties) projectProperties( properties, propertiesMat, eigenVecs, ...
                                     kernel, shift, scale, ...
                                     classes, classMeanProps );
  if nargout < 3
    classLabels = [];
  else
    meanClassDistance = pdist2( projectedProperties, classMeanProps );
    [~, minInd] = min( meanClassDistance, [], 2 );
    classLabels = classes(minInd);    
  end
  varargout = { classifyFunc, projectedProperties, classLabels };
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = projectProperties( newProperties, trainedProperties, ...
                                        eigenVecs, kernel, shift, scale,...
                                        classes, classMeanProps )
  numNew = size( newProperties, 1 );
  numTrained = size( trainedProperties, 1 );
  if ~isempty( shift )
    newProperties = bsxfun( @minus, newProperties, shift );
  end
  if ~isempty( scale )
    newProperties = bsxfun( @rdivide, newProperties, scale );
  end
  
  try
    % try to do vectorized first
    dotMat = kernel( newProperties, trainedProperties );
    assert( isequal( size( dotMat ), [numNew, numTrained] ) )
  catch
    % do nested for loops
    for row = numNew:-1:1
      for col = numTrained:-1:1
        dotMat(row,col) = kernel( newPropertiesRow, ...
                                  trainedProperties(col,:) );
      end
    end
  end
  projectedProps = dotMat * eigenVecs;
  if nargout > 1
    meanClassDistance = pdist2( projectedProps, classMeanProps );
    [~, minInd] = min( meanClassDistance, [], 2 );
    classLabels = classes(minInd);
  else
    classLabels = [];
  end
  varargout = { projectedProps, classLabels };
end
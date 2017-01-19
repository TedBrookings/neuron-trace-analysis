% varargout = GetClassifierScore( trainingMat, labels, varargin )
% Get a multi-class classifier score for each data point. Hold out one or
% more points and train a classifier on remaining points, then score the
% held out points on a nominally 0-1 scale for how close they are to each
% class' centroid.
% INPUTS
%  trainingMat: numData x numProperties matrix of data points to train
%  classifier on.
%  labels: numData x 1 list of labels of class identity (cell of strings or
%          array of numbers)
% OUTPUTS
%  classScores: numData x numClasses matrix of class scores
%  scoreFunc: function handle classScores = scoreFunc( testData )
%          input
%            testData: numTestData x numProperties matrix of untrained data
%          output
%            classScores: numTestData x numClasses matrix of class scores
function varargout = GetClassifierScore( trainingMat, labels, varargin )
  parser = inputParser();
  parser.addParameter( 'classifierFunc', [] )
  parser.addParameter( 'kernel', 'gaussian' )
  parser.addParameter( 'maxNumJackknife', 100 )

  parser.parse( varargin{:} )
  options = parser.Results;
  
  if isempty( options.classifierFunc )
    options.classifierFunc = ...
      @(mat, labels) DiscriminantAnalysis( mat, labels, ...
                                           'kernel', options.kernel );
  end
  
  numData = size( trainingMat, 1 );
  numChunks = min( numData, options.maxNumJackknife );
  chunkEdges = 1 + round( (0:numChunks) .* (numData / numChunks) );
  ProgressBar( 'Jackknifed classifier score', numChunks )
  for chunk = numChunks:-1:1
    low = chunkEdges(chunk); high = chunkEdges(chunk+1) - 1;
    classScores(low:high,:) = ...
      getJackknifeScore( low, high, numData, trainingMat, labels, ...
                         options );
    ProgressBar( 'Jackknifed classifier score' )
  end
  
  classifierFunc = options.classifierFunc( trainingMat, labels );
  scoreFunc = @(mat) getScores( mat, classifierFunc );
  
  varargout = { classScores, scoreFunc };
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function classifierScores = ...
      getJackknifeScore( low, high, numData, trainingMat, ...
                         labels, options )
  testInds = low:high;
  trainInds = [1:(low-1), (high+1):numData];
  scoreFunc = options.classifierFunc( trainingMat(trainInds,:), ...
                                     labels(trainInds) );

  [~, ~, classifierScores] = scoreFunc( trainingMat(testInds,:) );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function scores = getScores( testMat, classifierFunc )
  [~, ~, scores] = classifierFunc( testMat );
end
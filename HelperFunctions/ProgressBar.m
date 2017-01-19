function ProgressBar(name, totalTicks, varargin)
% ProgressBar(name, totalTicks, updatePeriod, keepUnderscores)
% Create or update a progress indicator.
% Note: in a parallel environment, uses a temporary file
%  -To create the progress indicator, pass the name, the total number of
%   ticks (updates to be completed) and if desired, specify the period to
%   update estimated completion time.
%  -To increment progress, pass only the name of an existing indicator.
%   The indicator automatically closes if progress is complete.
%  INPUTS:
%    name:  string, describing the progress indicator
%   OPTIONAL:
%    totalTicks:  the total number of progress updates until completion
%    updatePeriod:  interval (in seconds) between updated estimates of
%                   completion time.  Defaults to 5s.
%    keepUnderscores:  interpret underscores as characters, not subscripts.
%                      Defaults to true
%    useFirstJobName:  name the progress bar after the first job, otherwise
%                      name it 'Progress'. Default to false

if nargin == 1
  % increment progress bar
  incrementJobTicks(name);
else
  % define passable options and default values
  parser = inputParser();
  parser.addParameter('updatePeriod', 5)
  parser.addParameter('keepUnderscores', true)
  parser.addParameter('errorOnPreexisting', false)
  parser.addParameter('useFirstJobName', false)
  % get the options
  parser.parse(varargin{:})
  options = parser.Results;
  
  % monitor a new job
  monitorNewJob(name, totalTicks, options)
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function monitorNewJob(name, totalTicks, options)
% Monitor a new job
% 1. Information about the job is stored in a temporary file:
%   a) message (ID)
%   b) ticks (number of ticks so far)
% 2. A figure is created with
%   a) a waitbar-like figure
%   b) a timer
%   c) a callback function actived by the timer, which monitors the
%      temporary file and updates the figure appropriately

invalidChars = unique(name(regexp(name, '[^ .a-zA-Z_0-9\-]')));
if any(invalidChars)
  error('ProgressBar name ''%s'' contains invalid character(s): %s', ...
        name, invalidChars)
end

% create an empty timer file to store ticks, overwriting any existing file
timerFileName = [tempdir, name, '_ProgressBarTimer.txt'];
fid = fopen(timerFileName, 'w');
if fid < 0
  error('Couldn''t create temporary file: %s', timerFileName)
end
% close the file, leaving it empty
fclose(fid);

% create a structure to store job information
if options.keepUnderscores
  name = realUnderscores(name);
end
job = struct(...
  'message', name, ...
  'ticks', 0, ...
  'maxTicks', totalTicks, ...
  't0', tic, ...
  'tSinceUpdate', 0.0, ...
  'tRemaining', 'N/A', ...
  'timerFileName', timerFileName ...
);

% check if there is an existing progress timer
progressTimer = timerfind('Tag', 'ProgressBarTimer');
if ~isempty(progressTimer)
  % there is already a job being monitored, update the list of new jobs
  progressFigure = findobj('Tag', 'ProgressBar');
  newJobs = [get(progressFigure, 'UserData'), job];
  set(progressFigure, 'UserData', newJobs);
else
  % this is the first job, create a new figure/timer system

  existingFigs = findobj('Tag', 'ProgressBar');
  if ~isempty(existingFigs)
    % if a progress bar already exists, close it
    close(existingFigs)
  end
  
  % store the data necessary for updates
  if options.useFirstJobName
    if options.keepUnderscores
      updateName = realUnderscores(name);
    else
      updateName = name;
    end
  else
    updateName = 'Progress';
  end
  update = struct(...
    'name', updateName, ...
    'period', options.updatePeriod, ...
    'jobs', job...
  );
  % draw the figure
  update.figureHandle = drawBar(job, updateName);
  
  % create the timer and set the callback functions
  %   the callback function checks the temporary file and updates the
  %   figure appropriately. Once the last tick is called, each job closes;
  %   once the last job closes, the whole figure/timer system closes
  % update structure is stored in the timer's UserData
  t = timer('Tag', 'ProgressBarTimer', ...
            'Name', updateName, ...
            'timerfcn', @popupCallback, ...
            'stopfcn', @popupStopCallback, ...
            'period', update.period, ...
            'StartDelay', update.period, ...
            'ExecutionMode', 'FixedDelay', ...
            'UserData', update);
  start(t);
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function incrementJobTicks(name)
% Increment the ticks in the progress bar

% Increment ticks using a temporary file
% This method works even in a parallel environment (e.g. with parfor)
fileName = [tempdir, name, '_ProgressBarTimer.txt'];
fid = fopen(fileName, 'a');
if fid > 0
  fwrite(fid, ' ', 'char*1');
  fclose(fid);  
else
  error('couldn''t append to file %s', fileName)
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function popupCallback(obj, event, string_arg) %#ok<INUSD>
% monitor the progress via the temporary file, and update the progress bar
% appropriately

% get the update structure from the timer's (obj) UserData
update = get(obj, 'UserData');
jobs = update.jobs;
% get the new data
figureHandle = update.figureHandle;
if ~ishandle(figureHandle)
  % ProgressBar closed, just shut everything down
  stop(obj)
  return
end
% Add unique new jobs to be monitored
newJobs = get(figureHandle, 'UserData');
set(figureHandle, 'UserData', []);
anyNewJobs = false;
for job = newJobs
  % A very fast inner loop could have repeated, so must check for
  % uniqueness
  compareBools = strcmp(job.timerFileName, {jobs.timerFileName});
  if any(compareBools)
    % A new job overwrites / replaces an old one with the same timerFileName
    jobs(find(compareBools, 1)) = job;
  else
    % A completely new job has been added
    jobs = [jobs, job];
    anyNewJobs = true;
  end
end


for n = length(jobs):-1:1
  % get the number of ticks, as the length of the file
  fid = fopen(jobs(n).timerFileName, 'r');
  if fid >= 0
    fseek(fid, 0, 'eof');
    ticks = ftell(fid);
    if jobs(n).ticks == ticks
      jobs(n).tSinceUpdate = jobs(n).tSinceUpdate + update.period;
    else
      jobs(n).tSinceUpdate = 0.0;
      jobs(n).ticks = ticks;
    end
    fclose(fid);
  else
    % no timer file exists, so just delete this job
    jobs(n) = [];
    continue;
  end
  
  if jobs(n).ticks == jobs(n).maxTicks
    % This job is complete
    % delete the timer file
    if exist(jobs(n).timerFileName, 'file')
      delete(jobs(n).timerFileName);
    end
    if anyNewJobs
      % if there are new jobs waiting to take its place, delete it
      jobs(n) = [];
    end
  end
end

% save the list of jobs
update.jobs = jobs;
set(obj, 'UserData', update);

if isempty(jobs)
  % progress is complete, stop the timer
  stop(obj)
else
  % otherwise draw the progress bar figure
  drawBar(jobs, update.name, update.figureHandle);
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handle = drawBar(jobs, name, handle)

fontSize = 14;
pointsPerPixel = 72/get(0,'ScreenPixelsPerInch');

maxMessageWidth = 0;
for job = jobs
  maxMessageWidth = max(maxMessageWidth, length(job.message));
end

width = pointsPerPixel * fontSize * max(36, 0.6 * maxMessageWidth);
height = pointsPerPixel * fontSize * (2.0 + 4.0 * length(jobs));

if nargin < 3 || isempty(handle)
  handle = figure('Resize', 'off', ...
                  'MenuBar', 'none', ...
                  'Name', name, ...
                  'Tag', 'ProgressBar', ...
                  'NumberTitle', 'off', ...
                  'IntegerHandle', 'off', ...
                  'HandleVisibility', 'on', ...
                  'Interruptible', 'off', ...
                  'DockControls', 'off', ...
                  'BusyAction', 'queue');
  position = get(handle, 'Position');
  position(3) = width;
  position(4) = height;
  set(handle, 'Position', position)
  xlim([0 width])
  ylim([0 height])
  set(gca, 'Visible', 'off', ...
           'FontSize', fontSize,...
           'Position', [0.05, 0.05, 0.9, 0.9])
else
  set(0, 'CurrentFigure', handle);
  position = get(handle, 'Position');
  for c = get(gca, 'Children')
    delete(c)
  end
  if width > position(3)
    position(3) = width;
    xlim([0, width]);
  end
  if height > position(4)
    position(4) = height;
    ylim([0, height])
  end
  set(handle, 'Position', position)
end

dY = height / (2 * length(jobs) + 1);
y = height - dY * 0.75;
ticksX = pointsPerPixel * fontSize;
timeX = width - ticksX;
for job = jobs
  text(0, y, job.message)
  y = y - dY;
  barY = y - dY / 2;
  barX = width * job.ticks / job.maxTicks;
  rectangle('Position', [0, barY, width, dY], 'LineWidth', 1)
  if barX > 0
    rectangle('Position', [0, barY, barX, dY], 'FaceColor', 'r', ...
              'EdgeColor', 'none')
  end
  text(ticksX, y, sprintf('%d / %d', job.ticks, job.maxTicks))

  % estimate remaining time
  tNow = toc(job.t0);
  ticks = max(0.5, job.ticks / (1 - job.tSinceUpdate / tNow));
  tFinish = (job.maxTicks - ticks ) * tNow / ticks;
  tRemaining = [getTimeString(tFinish), ' remaining'];
  text(timeX, y, tRemaining, 'HorizontalAlignment', 'right')
  y = y - dY;
end
drawnow
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function popupStopCallback(obj, event, string_arg) %#ok<INUSD>
% shutdown the progress bar

% get the update structure
update = get(obj, 'UserData');
jobs = update.jobs;

% delete the timer
delete(obj)

% close the waitbar
if ishandle(update.figureHandle)
  jobs = [jobs, get(update.figureHandle, 'UserData')];
  close(update.figureHandle);
end

% delete any temporary files
for job = jobs
  if exist(job.timerFileName, 'file')
    delete(job.timerFileName);
  end
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function timeString = getTimeString(tFinish)
% convert a time in seconds into a human-friendly string
s = round(tFinish);
if s < 60
  timeString = sprintf('%ds', s);
  return
end
m = floor(s/60);
s = s - m * 60;
if m < 60
  timeString = sprintf('%dm%ds', m, s);
  return
end
h = floor(m/60);
m = m - h * 60;
if h < 24
  timeString = sprintf('%dh%dm%ds', h, m, s);
  return
end

d = floor(h/24);
h = h - d * 24;
timeString = sprintf('%gd%dh%dm%ds', d, h, m, s);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function texSafeString = realUnderscores(titleStr)
%  Protects the underscores in titleStr from being interpreted as
%  subscript commands by inserting a '\' before them.

ind = strfind(titleStr, '_');
texSafeString = '';
previous = 1;
for n = 1:length(ind)
  if ind(n) > 1
    texSafeString = [texSafeString, titleStr(previous:(ind(n)-1))]; %#ok<*AGROW>
  end
  texSafeString = [texSafeString, '\'];
  previous = ind(n);
end
texSafeString = [texSafeString, titleStr(previous:end)];
end

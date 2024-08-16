function pyplot(varargin)
    % PYPlot A drop-in replacement for MATLAB's plot function using Python's matplotlib.
    %
    % Usage:
    %   pyplot(x1, y1, 'Style1', x2, y2, 'Style2', ..., 'Title', 'My Plot', 'XLabel', 'X-axis', 'YLabel', 'Y-axis', 'Legend', {'Label1', 'Label2', ...})
    %
    % Parameters:
    %   x1, y1, x2, y2, ... - X and Y data for multiple lines
    %   'Style1', 'Style2', ... - Line styles for each line (required)
    %   'Title' - (Optional) Title of the plot
    %   'XLabel' - (Optional) Label for the X-axis
    %   'YLabel' - (Optional) Label for the Y-axis
    %   'Legend' - (Optional) Cell array of legend labels

    % Parse optional parameters
    p = inputParser;
    addParameter(p, 'Title', '', @(x) ischar(x) || isstring(x));
    addParameter(p, 'XLabel', '', @(x) ischar(x) || isstring(x));
    addParameter(p, 'YLabel', '', @(x) ischar(x) || isstring(x));
    addParameter(p, 'Legend', {}, @(x) iscell(x) && all(cellfun(@(y) ischar(y) || isstring(y), x)));

    % Find the position of the first optional parameter
    optParamIdx = find(cellfun(@(x) ischar(x) && any(strcmp(x, {'Title', 'XLabel', 'YLabel', 'Legend'})), varargin), 1);
    if isempty(optParamIdx)
        optParamIdx = nargin + 1;
    end

    parse(p, varargin{optParamIdx:end});

    titleStr = char(p.Results.Title); % Convert to character array
    xlabelStr = char(p.Results.XLabel); % Convert to character array
    ylabelStr = char(p.Results.YLabel); % Convert to character array
    legendLabels = p.Results.Legend;

    % Extract data and styles
    numArgs = optParamIdx - 1;
    numLines = numArgs / 3;
    x = cell(1, numLines);
    y = cell(1, numLines);
    styles = cell(1, numLines);

    for i = 1:numLines
        x{i} = varargin{(i - 1) * 3 + 1};
        y{i} = varargin{(i - 1) * 3 + 2};
        styles{i} = varargin{(i - 1) * 3 + 3};
    end

    % Save data and parameters to a .mat file
    dataFile = 'temp_plot_data.mat';
    save(dataFile, 'x', 'y', 'styles', 'titleStr', 'xlabelStr', 'ylabelStr', 'legendLabels');

    % Call the Python plotting script
    system(['python3 plot_script.py ', dataFile, '&']);
end

currpath = pwd();

% This line sets the sabl/src/ code path
addpath(genpath([currpath, '/src/']));

% This line adds sabl/src and adigator paths to the Matlab default paths
savepath;

fprintf([ ...
         '     The automatic_differentiation module in Sabl requires ADiGator. \n', ...
         '     Please visit https://github.com/matt-weinstein/adigator to download it. \n', ...
         '     This can be omitted, but the automatic_differentiation tests should be disabled in that case. \n']);

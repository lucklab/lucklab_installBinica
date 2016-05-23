% Unzip 'Binica_mac_intel.zip' to eeglab directory
binica_zip = 'Binica_mac_intel.zip';
[path_eeglab, ~] = fileparts(which('eeglab'));
system(['unzip ' binica_zip , ' -d ' , path_eeglab ])
system(['rm -rf ', path_eeglab '/__MACOSX'])

% Copy `startup.m` to Matlab directory
path_matlab      = userpath();           % Get Matlab's startup directory
path_matlab      = path_matlab(1:end-1); % remove colon (:) at end of string
filename_startup = 'startup.m';
if (exist(fullfile(path_matlab, filename_startup), 'file'))
    % Append startup info
else
    % Copy startup.m file to Matlab startup directory
    system(['cp ', filename_startup, ' ', path_matlab]);
end
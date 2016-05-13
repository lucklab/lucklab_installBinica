% 1. Download BINICA appropriate for your computer Mac OS X PPC, Mac OS X Intel,
% 
% Source URL:http://sccn.ucsd.edu/wiki/Binica
% 
%% In MATLAB:
%      Locate Matlab?s Bash path

[status, result] = system('echo $PATH'); display(result);
% eval('!echo $PATH')
% getenv('PATH');


fullfile(pwd, 'binica_osx_fat')
 
% Copy the ica_osx file to one of those paths
% Terminal: sudo cp ica_osx /bin/

%      Check to make sure ica_osx function is now available in Matlab
%      eval('!ica_osx')

%      if not
%          /bin/bash: ica_osx: command not found
%      if yes

% To use one of these programs from within Matlab (and EEGLAB)
% download the file and place them in the eeglab directory (you may create a subfolder for them or uncompress them in the function subfolder).
% edit the icadefs.m file to specify the file name of the executable you intend to use (variable ICABINARY).
% add the path to the binary function both to your Matlab path and to your Unix path, otherwize the system will not be able to locate the executable file.
% From the command line, you may use the binica() function that will call the binary executable. From the EEGLAB graphical interface, run ICA using the 'binica' option of the Tools > Run ICA graphic interface (see the tutorial for how to compute ICA components).
%%
% 
% $$e^{\pi i} + 1 = 0$$
% 
%   for x = 1:10
% 
%  PREFORMATTED
% 
% * ITEM1
% * ITEM2
% 
%  TEXT
% 
%       disp(x)
%   end
% 
% 


% To modify the system path across shell and MATLAB sessions, add the following commands to the MATLAB startup file as described in Startup Options in MATLAB Startup File.

path1 = getenv('PATH')
path1 = [path1 ':/Users/jtkarita-admin/Desktop/erplab_installBinica/binica_osx_fat']
setenv('PATH', path1)
system('echo $PATH')
system('ica_osx')
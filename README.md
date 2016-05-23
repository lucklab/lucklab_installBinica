# Install BINICA

To install:

- In Matlab, un LuckLab BINICA Installer script: `install_binica`
- Edit the `icadefs.m` file to point to the `ica_osx` binary file.
- `.../eeglab13_x_x/functions/sigprocfunc/icadefs.m` (line 119)
  - ```matlab
    ...
    % INSERT location of ica executable (LINUX ONLY) for binica.m below
    eeglab_p = fileparts(which('eeglab'));
    ICABINARY = 'ica_osx';  % <<< EDIT HERE
    
    try
       set(0,'defaultaxesfontsize',AXES_FONTSIZE);
    ...

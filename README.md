# Install BINICA

To install:
- Download the [`lucklab_installBinica` folder](https://github.com/lucklab/lucklab_installBinica) to your computer
- In Matlab, change directory to the `lucklab_installBinica` folder
- Run the `install_binica` script
  - *This will unzip the `binica_osx_fat` folder into your EEGLAB folder, and add it to your Matlab path and Unix path*
- Edit `icadefs.m`:  Set the `ICABINARY` variable to `ICABINARY = fullfile(eeglab_p, 'binica_osx_fat', 'ica_osx');` binary file.
  - To open `icadefs.m`: `eeglab; open('icadefs.m');` 
  - Example: `.../eeglab13_x_x/functions/sigprocfunc/icadefs.m` (line 119)
    ```matlab
      ...
      % INSERT location of ica executable (LINUX ONLY) for binica.m below
      eeglab_p = fileparts(which('eeglab'));
     ICABINARY = fullfile(eeglab_p, 'binica_osx_fat', 'ica_osx');;  % <<< EDIT HERE
      
      try
         set(0,'defaultaxesfontsize',AXES_FONTSIZE);
      ...
      ```
- Quit and restart Matlab

# Install BINICA

To install:
- download the file and place them in the eeglab directory (you may create a subfolder for them or uncompress them in the function subfolder).
- edit the icadefs.m file to specify the file name of the executable you intend to use (variable ICABINARY).
- add the path to the binary function both to your Matlab path and to your Unix path, otherwize the system will not be able to locate the executable file.
- From the command line, you may use the binica() function that will call the binary executable. From the EEGLAB graphical interface, run ICA using the 'binica' option of the Tools > Run ICA graphic interface (see the tutorial for how to compute ICA components).


- copy `startup.m` and `binica_osx_fat` to your Matlab startup directory
  - `userpath()`
- edit the icadefs.m file to specify the file name of the executable you intend to use 
 - `eeglab; open('icadefs.m');`
 - Change the variable ICABINARY to point to `ica_osx` (and not `ica_linux`)
 - `ICABINARY = 'ica_osx'`

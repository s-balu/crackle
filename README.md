Crackle is a fully C-based version of the Grackle cooling and chemistry code.  Information on Grackle is available here:
https://grackle.readthedocs.io/en/latest/index.html

Currently crackle ONLY works with Ewan Jones's kiara branch of Grackle, which includes the subgrid ISM dust+H2 model described here: 
https://ui.adsabs.harvard.edu/abs/2024MNRAS.535.1293J/abstract
Ewan's Grackle branch is here:
https://github.com/EwanBJones98/grackle
Don't forget to checkout the kiara branch.

To use Crackle, replace/add the codes in this repo to those in the above Grackle repo under src/clib (see below).  Then when putting it in your C code, set use_grackle=2.  The API is otherwise the same.

INSTALLATION INSTRUCTIONS

Clone crackle into your /src/clib/ directory of grackle.  It should then appear in the crackle subdirectory.

There are 6 files in grackle that crackle overwrites, in addition to adding many new ones.  From the /src/lib  directory, run the following commands to archive them in a separate directory, and then make symbolic links to (or copy over) all the crackle files:

mkdir grackle_orig

mv ./calculate_cooling_time.c grackle_orig

mv ./calculate_temperature.c grackle_orig

mv ./calculate_pressure.c grackle_orig

mv ./grackle_chemistry_data.h grackle_orig

mv ./grackle_types.h grackle_orig

mv ./solve_chemistry.c grackle_orig

ln -s crackle/* .

Then proceed with making grackle as usual.  
To update crackle, just git pull the changes into the crackle subdirectory and all should proceed as normal (if you've made symlinks).

NCI-Gadi specific install instruction:

Edit the `src/clib/Make.mach.gadi` file by changing the `MACH_INSTALL_PREFIX` to the location of grackle_kiara.
Then install grackle_kiara as:
```
cd grackle
./configure
cd src/clib/
## edit Make.mach.gadi
## edit path to install library MACH_INSTALL_PREFIX = $(HOME)/<grackle location>#
make machine-gadi
make show-config
make clean
make
make install
```

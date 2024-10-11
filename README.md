Crackle is a fully C-based version of the Grackle cooling and chemistry code.  Information on Grackle is available here:
https://grackle.readthedocs.io/en/latest/index.html

Currently crackle ONLY works with Ewan Jones's kiara branch of grackle, which is unfortunately private until Ewan gets his thesis papers out.  However, it should be straightforward to adapt it to the main grackle branch.  If I have time at some point I will try to do this, but I will also happily accept pull requests:)

To use Crackle, replace/add the codes in this repo to those in the Grackle repo under src/clib (see below).  Then when putting it in your C code, set use_grackle=2.  The API is otherwise the same.

INSTALLATION INSTRUCTIONS

Clone crackle into your /src/clib/ directory of grackle.  It should then appear in the crackle subdirectory.

There are 5 files in grackle that crackle overwrites, in addition to adding many new ones.  From the /src/lib  directory, run the following commands to archive them in a separate directory, and then make symbolic links to (or copy over) all the crackle files:

mkdir grackle_orig
mv ./calculate_cooling_time.c grackle_orig
mv ./calculate_temperature.c grackle_orig
mv ./grackle_chemistry_data.h grackle_orig
mv ./grackle_types.h grackle_orig
mv ./solve_chemistry.c grackle_orig
ln -s crackle/* .

Then proceed with making grackle as usual.  
To update crackle, just git pull the changes into the crackle subdirectory and all should proceed as normal (if you've made symlinks).

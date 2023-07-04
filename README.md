# Apollonian
Methods for computing various things related to (primarily integral) Apollonian circle packings.

## Installation
System requirements:
* Linux/Mac: none that I am aware of;
* Windows: you need to be running PARI/GP through WSL. See this [tutorial](https://pari.math.u-bordeaux.fr/PDF/PARIwithWindows.pdf) for how to set this up.

After cloning into the repository, call "make". If your PARI/GP is not installed in the default location of "/usr/local", then find the folder X, where:
* X/lib contains the file "libpari.so";
* X/lib/pari contains the file "pari.cfg";
* X/include/pari contains around 18 .h header files for PARI/GP.
Searching for the file pari.cfg should be sufficient to find this folder.

Afterwards, make a new file called "pari_loc.txt", containing the path to X (which should start with "/"). Now you can call ``make'' and it will compile properly.

## Using the package
Start the program with "gp apol", and call "?apol" to access the help.

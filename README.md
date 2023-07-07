# Apollonian
Methods for computing various things related to (primarily integral) Apollonian circle packings. This includes:
* Methods to generate integral packings, and perform the basic operations (Apollonian moves, reduction, etc.);
* Automatic generation of pictures of circle packings in LaTeX;
* Finding (all) occurences of a given curvature in a packing;
* Finding all curvatures in a given range in a packing.

Methods relating to searching for curvatures have been optimized for speed and memory as much as possible (aside from parallelization).
In addition, there are a few methods related to quadratic forms, Möbius maps, and data collection in the package.

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

## Papers

This package has been used to generate the data for several papers:
* [The Apollonian staircase](https://academic.oup.com/imrn/advance-article/doi/10.1093/imrn/rnad007/7072814?utm_source=authortollfreelink&utm_campaign=imrn&utm_medium=email&guestAccessKey=72a0a7b9-45f7-47f0-900c-b13d85dac729)
* [The Local-Global Conjecture for Apollonian circle packings is false](https://arxiv.org/abs/2307.02749): data on missing curvatures is stored in the GitHub repository [Apollonian Missing Curvatures](https://github.com/JamesRickards-Canada/Apollonian-Missing-Curvatures).

## Bibtex

If you use this code in a project, please let me know! A suggested Bibtex entry would be:
```
@misc{Apollonian,
	AUTHOR = {Rickards, James},
	TITLE = {Apollonian},
	YEAR = {2023},
	PUBLISHER = {GitHub},
	JOURNAL = {GitHub repository},
	HOWPUBLISHED = {\url{https://github.com/JamesRickards-Canada/Apollonian}},
}
```

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

You need to know where the version of PARI/GP you want to use is installed, in particular, the file ```pari.cfg```. The default location is ```/usr/local/lib/pari/pari.cfg```.
* If this is the correct location, call ```make``` to build the project.
* Otherwise, call ```make setup``` to search for the location of the file. By default the program searches in ```/usr```, but there is a chance it is not installed there (this sometimes happens on a server). If this is the case, you can supply an alternate location.
* If the program finds potential matches, it will ask you to confirm which files are correct, and saves them to ```paricfg_loc.txt```. Once this step completes, a call to ```make``` will compile the project! Modifying the program (e.g. via ```git pull```) won't require redoing this setup, unless the version of PARI/GP you use changes.

## Using the package
Start the program with ```gp apol```, and call ```?apol``` to access the help.

## Papers

This package has been used to generate the data for several papers:
* [The Apollonian staircase](https://academic.oup.com/imrn/advance-article/doi/10.1093/imrn/rnad007/7072814?utm_source=authortollfreelink&utm_campaign=imrn&utm_medium=email&guestAccessKey=72a0a7b9-45f7-47f0-900c-b13d85dac729): I have since overhauled this repository, so the code relevant to this project is temporarily absent. Rolling back to commit 0eab06e54ce02f367c6f24b0fff1862b789f1633 (on June 1st 2023) will give access to the methods used to generate the data for this paper. I will also eventually add them back in, whereupon this warning will be deleted.
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

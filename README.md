
# LIMA

Hello there!
My name is Daniel Johansen, and I am the sole author of this repository.


**If you just want to see the cool stuff, go to Quantom->Quantom->engine.cu/bodies.cuh. Herein lies the main components which makes the simulations go fast**


## LIMA backstory

LIMA is a collection of programs, centered around LIMA_MD which is a GPU-based Molecular Dynamics simulator, close to capable of keeping up with GROMACS - in terms of speed only!

LIMA was developed in 1 year, as part of the Masters Thesis of Computer  Science at University of Southern Denmark (SDU).
The long-term goal of LIMA is to become a MD cloud-provider for MD users.

If you want to read about LIMA_MD, or how to make an efficient MD simulator, I will upload my Thesis here August 2022.
LIMA_MD takes 3.34 ms/step compared to GROMACS 0.27 with the same setup. Disabling 2 temporary and poorly implemented features (Explicit Solvation and Bonded Lennard-Jones Verification), LIMA_MD reaches an impressive 0.17 ms/step!




## How to install/run

I have only had luck using this on Manjaro linux. It should be as easy as cloning the repo, and running install.sh from the root directory

LIMA_MD is the main program. The old development name was Quantom, and you will still see that name around the repo for VisualStudio reason..









If you want to get in touch with me, feel free to write me at daniel@lima-dynamics.com

// Daniel

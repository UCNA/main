The purpose of this file is to assist in downloading, compiling, and using the geant4 simulations which are on the github.

Thanks to Michael Mendenhall for all your help in getting this up and running. I hope this helps others who are interested. Also feel free to make comments/suggestions for things I may have missed and also problems others run into. The more detailed this document becomes the better.

Michael Brown
August 2013

/////////////////////////////////
///////// Downloading ///////////
/////////////////////////////////

If you haven't opened this file from the github initially, go to https://github.com/hickerson/UCNA. Here you will find a repository of all pertinent code used for analysis and simulation. 

I suggest downloading the ZIP file of this entire directory using the "Download ZIP" button on the right side of the screen.

Once you have downloaded and unzipped this file into a directory of your choosing, type 'make' in the terminal while in this directory. Note: I did this as I was unsure as to what the simulations were going to need. This may not be necessary, as it seems all the code required for running the simulation is within ucnG4_dev.

/////////////////////////////////
////////// Compiling ////////////
/////////////////////////////////

Once you have the files downloaded, cd into the directory ucnG4_dev. Here you find the simulation code. ucnG4_prod.cc houses the main() function. The rest of the code is in the /src directory and headers in include/. 

As suggested in the geant4 manual, create a build directory wherever you would like. Using CMAKE, you can then generate the appropriate Makefile to compile the code.

In the build directory you choose, type the following:

cmake -DGeant4_DIR=/path/to/Geant4_install/lib64/Geant4-9.5.1 /path/to/ucnG4_dev

The first path is to your installation of Geant4. The second path is that to the simulation code. Running this will automatically create the Makefile in your build directory. Once this is done, type 'make' in the build directory. This should compile the code and create executables in the /bin directory.

////////////////////////////////
////////// Running /////////////
////////////////////////////////

Note:

You need to set some environment variables prior to running the simulation. (note these are for my machine. Replace the paths with those for yours)

for csh or tcsh shells:

setenv G4WORKDIR /tmp_mnt/home62/grads/mabrow05/UCNA/geant4work
setenv G4BINDIR /tmp_mnt/home62/grads/mabrow05/UCNA/mb_sims/ucnG4_dev-build/bin

for bash shells:

export G4WORKDIR=/tmp_mnt/home62/grads/mabrow05/UCNA/geant4work
export G4BINDIR=/tmp_mnt/home62/grads/mabrow05/UCNA/mb_sims/ucnG4_dev-build/bin

The first directory is to where you want the geant output to be saved. The second is where the binaries are held (the bin directory within your build directory).....


First things first, you may want to check if your code is all compiled correctly and runs. Using the generic geant4 macro test.mac (given to me by Michael Mendenhall), you can simply change the local environment variables within this file and run:

./ucnG4_prod <macro filename> [physics list]

Here physics list is not mandatory as it will default to Livermore. This was more of a check to see that everything is working. Now all you have produced from this are raw root files of the simulated data. To actually work with the simulated data in terms of pertinent plots, you need to run the UCNA_MC_Analyzer program:

./UCNA_MC_Analyzer <filename containing list of raw root files> <output root file name> [saveall|undead]

You literally need a file with a list of the raw root files you want to analyze, ie. list.txt. In this file you may only have one raw root file listed in this case, since you only generated one when you ran the test macro.

Now you can either manipulate this macro to do what you want, or you can utilize the more complex (but robust) methods which were developed in previous analyses... (this is below)

****** GeantSimManager.py ******

Using the script GeantSimManager.py will allow you to create a list of macros, run jobs, and then analyze the raw root files, all with one keystroke. It is very convenient once you have it working.I suggest that if you aren't familiar with python, you read up a little on class and function definitions within python. Also, take the time to read through this python script and understand the process.

You need to go to the downloaded directory /Scripts which came with the directories from the Github. Within these scripts you will find:

GeantGenMacroTemplate.mac
GeantJobTemplate.sub (not required as this was never implemented)
GeantSimManager.py

The above environment variables take care of the paths needed by these scripts. By running GeantSimManager.py, the script will create a geant4 macro which runs in geant4 and saves all the pertinent information, while also running the UCNA_MC_Analyzer program on the raw data. All the information from the Analyzer will be saved where you set the geant4 work directory.

Make sure you have GNU Parallel installed on your machine prior to using GeantSimManager.py method. This is utilized within the code to use all the cores on your machine simultaneously. I found that going to the launch_sims function and setting the nClusters equal to the number of cores I had (or wanted to delegate to simulation) was the easiest way to manage my cores. 

Another important note.... You need to note what shell you are using. In line 208, within the launch_sims function, some output is piped to the "joblog". Depending on the shell, this syntax changes:

jobsout.write(onejob+" > %s 2>&1\n"%self.settings["joblog"]) (bash)
			or
jobsout.write(onejob+" >&  %s \n"%self.settings["joblog"])   (csh)

I had to add the second line into the code once I realized this problem.

Now, all the options for what type of simulation you would like to run is in the main function at the bottom. Here there are a series of IF statements. Each one represents a different subset of simulation processes. You choose which one you want to run, and set this if statement to 1 rather than 0. This is just evaluating 0 and 1 as a simple boolean expression (0 being false and 1 being true).

From here you just go to your terminal and run ./GeantSimManager.py. 

You can then go to your geant4 work directory and go to output. Here you will find all the analyzed files as well as the raw root files. Have Fun!



# ACTI-Model
The ACTI Model (Agent-based Ctl-Tumor Interaction Model) 

The ACTI Model is an agent-based model that simulates interaction between Cytotoxic T-Lymphocytes and tumor cells. It is callibrated on melanoma cells and conventionnal target cells in collaboration with Pr. Salvatore Valitutti (CRCT - INSERM) and the results of the study it has been used is published in [Nature Scientific Reports](https://www.nature.com/articles/s41598-019-48711-2). If you use this model, please cite the paper :)

# Installation
## MecaCell
The ACTI Model is based on the open-source simulation platform MecaCell developed in our research group. It is available on [github](https://github.com/jdisset/MecaCell). Please follow the installation instruction of the master branch before proceed.
## ACTI Model
You will need cmake 2.8+ and a recent C++ compiler (C++11 support is required). QT5.2+ is necessary if you want a graphical interface.
  * Clone this repository (git clone git@github.com:scussatb/ACTI-Model.git)
  * In the ACTI-Model directory, create a directory "build"
  * In this new directory, run 
    cmake ..
    make

# Run the model
To run the model, in the build directory:
./src/mecaBase ctlTumorRatio killDur inhRad inhProb simDur tumTyp
where
  * ctlTumorRatio is the ratio ctl/tumor in the experiment
  * killDur is the duration of 1 kill for a CTL
  * inhRad is the inhibition radius when a tumor cell is killed
  * inhProb is the probability for a CTL to be inhibited if a tumor dies in the inhibition radius
  * simDur is the duration of the simulation in hours
  * tumType is the tumor cell line (D10 or JY)


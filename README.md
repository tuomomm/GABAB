# GABAB

Biochemically detailed model of GABAB activation for NEURON simulator.  
Tuomo Maki-Marttunen, 2022-2023  

This entry contains the python scripts needed for running a biochemically
detailed model of GABAB activation. The reaction rate parameters have been
fit to data on GABABR-mediated effects (both pre- and post-synaptic), see
Isaacson et al. 1993 ("Local and diffuse synaptic actions of GABA in the
hippocampus") and Olpe et al. 1994 ("Contribution of presynaptic GABA-B
receptors to paired-pulse depression of GABA-responses in the hippocampus").

Files included:  
- IC_singlecompartment.xml       #Initial concentrations (non-fitted) in NeuroRD format
- Reactions.xml                  #Reactions reaction rates (non-fitted) in NeuroRD format
- drawfig_fit.py                 #Python file to run the simulation
- makeNeuronModelExtfilename.py  #Python file to convert the NeuroRD model into NEURON model
- mesh_general.out               #A mesh-file (NeuroRD format) needed for determining volume
- model_nrn_extfilename.py       #The model converted into a NEURON model
- mytools.py                     #Python file with some generic tools

To run the model, first convert the NeuroRD format into a NEURON model (this
has been done, but should be redone if some initial concentrations or reaction
rates are changed):  

``python3 makeNeuronModelExtfilename.py``

Then, run the model and draw the results with the following script:  

``python3 drawfig_fit.py``

The script runs three simulations, one post-synaptic simulation (here,
of interest are the amounts of Gibg-bound GIRK channels) and two pre-
synaptic simulations (here, of interest are the amount of Gibg-bound
vs non-bound VGCCs). The script outputs four important lines:  

paramdict = {'k[0]': 111.531, 'k[1,9]': 708.687, 'k[3,4,5,6,7,8]': 165.616, 'k[15,17,19,21,23]': 10.797, 'k[16,18,20,22,24]': 93.218, 'OnlyExp0_RGS': 2.903, 'OnlyExp1_RGS': 1.728, 'OnlyExp2_RGS': 2.235, 'gaba_flux': 2797.567}  
- This is the parameter dictionary, showing the alterations of each
  parameter as obtained from a multi-objective optimization (not
  described here)  

python3 model_nrn_extfilename.py 1010000.0 1e-07 1000000.0 1 1 50.0 2797.567 1 1 None VGCC,RGS 0.0,2.903 0,1,9,2,10,3,4,5,6,7,8,15,17,19,21,23,16,18,20,22,24 111.531,708.687,708.687,86.60863827,86.60863827,165.616,165.616,165.616,165.616,165.616,165.616,10.797,10.797,10.797,10.797,10.797,93.218,93.218,93.218,93.218,93.218 fitXXXXXXXXXX_0.mat 5  
- This is the script for running the post-synaptic simulation (the command is not only printed
  but also executed inline in the script drawfig_fit.py). The command describes the simulation
  parameters (length of simulation, start of stimulation, etc.), including the parameter
  alterations (processed values of paramdict). The results will be saved to fitXXXXXXXXXX_0.mat.
  In the post-synaptic simulation, no VGCCs are present (VGCC concentration is multiplied by 0).  

python3 model_nrn_extfilename.py 1010000.0 1e-07 1000000.0 1 1 50.0 2797.567 1 1 None GIRK,RGS 0.0,1.728 0,1,9,2,10,3,4,5,6,7,8,15,17,19,21,23,16,18,20,22,24 111.531,708.687,708.687,86.60863827,86.60863827,165.616,165.616,165.616,165.616,165.616,165.616,10.797,10.797,10.797,10.797,10.797,93.218,93.218,93.218,93.218,93.218 fitXXXXXXXXXX_0.mat 5  
- The same as above, but with slightly different parameters (no GIRK channels but a default
  amount of VGCCs)  

python3 model_nrn_extfilename.py 1010000.0 1e-07 1000000.0 1 1 50.0 2797.567 1 1 None GIRK,RGS 0.0,2.235 0,1,9,2,10,3,4,5,6,7,8,15,17,19,21,23,16,18,20,22,24 111.531,708.687,708.687,86.60863827,86.60863827,165.616,165.616,165.616,165.616,165.616,165.616,10.797,10.797,10.797,10.797,10.797,93.218,93.218,93.218,93.218,93.218 fitXXXXXXXXXX_0.mat 5  
- The same as above, but with slightly different parameters (different amount of RGS proteins)  

The simulation results will be saved in fig_fit.eps (summary of the three simulations, plotted
 against the experimental data), and fig_fit_allspecies_syn?.eps (three simulations plotted
 separately, species by species).






#A script for simulating the GABAB model in three scenarios:
#  1) A post-synaptic terminal-like compartment
#  2) A pre-synaptic terminal-like (I->E-type) compartment
#  3) A pre-synaptic terminal-like (E->E-type) compartment
#  The script runs the model with an optimized parameter set and plots the results (blue), comparing to experimental data (red)
#    
#Tuomo Maki-Marttunen, 2019-2023, CCBY 4.0

import matplotlib
matplotlib.use('Agg')
from pylab import *
import scipy.io
import os
import time
from os.path import exists
import re
import random
from scipy.optimize import curve_fit
import mytools
import pickle

fixedBlock = 'Leak,AC5'
fixedBlockCoeff = '1.33,0.0'
fixedAltered = '1'
fixedAlteredCoeff = '1'

mesh_input_file = open('mesh_general.out','r')
mesh_firstline = mesh_input_file.readline()
mesh_secondline = mesh_input_file.readline()
mesh_values = mesh_secondline.split()
my_volume = float(mesh_values[-2])*1e-15 #litres
mesh_input_file.close()

#These are fixed across experiments:
Duration         = 1e6+10000
tolerance        = 1e-7
gaba_input_onset   = 1e6
gaba_input_N       = 1
gaba_input_freq    = 1
gaba_input_dur     = 50.0 #Long duration due to assumed volume transmission
gaba_input_flux    = 5000.0 # Onset time constant (post-syn) should depend on concentration. This maybe corresponds to a saturating concentration
Ntrains          = 1
trainT           = 1
record_T         = 5

#Data for PPR
data_exppr_scaled = [[147.51475085724883, 0.7605217633066633], #Local and diffuse synaptic actions of GABA in the hippocampus. J.S.Isaacson, J.M.Solis, R.A.Nicoll, 1993
                         [198.375088814062, 0.7229294739118345],
                         [299.67563559976526, 0.664576318309598],
                         [497.7047357202435, 0.7310354947329398],
                         [693.0153532482777, 0.8064038800160637],
                         [892.7743968366749, 0.903557196255908],
                         [1295.6040900806277, 0.9651926106700441],
                         [1591.7704117883293, 1.0000293472552593]]
data_exppr_ydata_scaled = [data_exppr_scaled[i][1] for i in range(0,len(data_exppr_scaled))]
data_exppr_xdata = [data_exppr_scaled[i][0] for i in range(0,len(data_exppr_scaled))]
data_antiexppr_ydata_scaled = [1-data_exppr_ydata_scaled[i] for i in range(0,len(data_exppr_ydata_scaled))]                                                              #1-PPR(control)/PPR(GABAB_blocked) should closely reflect the temporal dynamics of GABAB activation

data_inppr_xdata = [350, 550, 750, 1000, 2000, 3000, 4000]
data_inppr_ydata_nonscaled = [0.4978165938864628, 0.6462882096069869, 0.6965065502183405, 0.7467248908296943, 0.8275109170305677, 0.8668122270742359, 0.8973799126637554,  ] #https://link.springer.com/content/pdf/10.1007/BF00169135.pdf (Olpe et al. 1994)
data_inppr_ydata_GABABblocked = [0.8078602620087335, 0.7991266375545851, 0.7991266375545851, 0.8296943231441048, 0.8362445414847161, 0.868995633187773, 0.9301310043668121 ] #
data_inppr_ydata_scaled = [data_inppr_ydata_nonscaled[i]/data_inppr_ydata_GABABblocked[i] for i in range(0,len(data_inppr_ydata_nonscaled))]                                 #Rationale: Certain fraction of PPR is caused by GABAB and the rest by vesicle depletion etc.
data_antiinppr_ydata_scaled = [1-data_inppr_ydata_scaled[i] for i in range(0,len(data_inppr_ydata_scaled))]                                                                  #1-PPR(control)/PPR(GABAB_blocked) should closely reflect the temporal dynamics of GABAB activation

data_postsyn_scaled = [[77.25316584201036, 0.12989552390326564],
                       [97.62754955662777, 0.27964900975021945],
                       [120.74912053813799, 0.5282774023543535],
                       [138.750612423254, 0.7288741564007608],
                       [156.63910945934614, 0.8842729708376185],
                       [196.90764878022938, 0.9916886991909393],
                       [234.23125984842804, 0.9211331264703673],
                       [271.3853786430907, 0.7827806443354708],
                       [311.0253841162788, 0.6387828336106707],
                       [355.6724738368372, 0.4976187218340474],
                       [392.88309005601195, 0.38186520950392605],
                       [432.6007794879039, 0.26894098226069135],
                       [469.88907966578245, 0.18426105341213533],
                       [577.1282535675034, 0.07993061410052119],
                       [679.5227733173257, 0.03773851402945802],
                       [731.8817614837636, -0.018666219395389358]]
data_postsyn_xdata = [data_postsyn_scaled[i][0] for i in range(0,len(data_postsyn_scaled))]
data_postsyn_ydata_scaled = [data_postsyn_scaled[i][1] for i in range(0,len(data_postsyn_scaled))]

nrnfactor = 6.022e23*my_volume*1e-9 #how many molecules is 1nM in the given volume (0.5 fl)
experiments_blockeds = [ [['VGCC', 0.0]], [['GIRK', 0.0]], [['GIRK', 0.0]] ] #Three experiment types, where in the first one the only concentration change is setting VGCC to zero, and the second and third one where the only change is setting GIRK to zero
experiments_altereds = [ ['', ''], ['', ''], ['', ''] ]                      #All rates same in the two experiments


#A function that runs all three simulations. This same function, combined with an error function that calculated the difference between data and predicted dynamics of the GABAB-mediated entity, was used when fitting the parameters.
def run_model(parameters,deleteFiles=True,rankID=0):
  data = []
  gaba_flux = gaba_input_flux
  thisAltered = ''
  thisAlteredCoeff = ''
  thisBlock = ''
  thisBlockCoeff = ''
  paramkeys = list(parameters.keys())
  experimentwiseBlockeds = []

  #Search the parameters from the parameter dictionary and determine how to change the reaction rate parameters
  for iparam in range(0,len(paramkeys)):
    if paramkeys[iparam].find('k[') > -1:
      iks = paramkeys[iparam][2:-1].split(',')
      thisAltered = thisAltered + ','.join(iks)+','
      thisAlteredCoeff = thisAlteredCoeff + ','.join([str(parameters[paramkeys[iparam]]) for i in iks])+','
      if paramkeys[iparam].find('k[1,9]') > -1: #since we know the affinity of gaba + GABABR <-> gabaGABABR, use it to calculate the backward rate
        thisAltered = thisAltered + '2,10,'
        thisAlteredCoeff = thisAlteredCoeff + str(5.555*parameters[paramkeys[iparam]]*0.00011/0.005)+','+str(5.555*parameters[paramkeys[iparam]]*0.00011/0.005)+','
    elif paramkeys[iparam].find('gaba_flux') > -1:
      gaba_flux = parameters[paramkeys[iparam]]
    elif 'OnlyExp' in paramkeys[iparam]:
      blockediexp = int(paramkeys[iparam][7])
      blockeds = paramkeys[iparam][9:].split(',')
      experimentwiseBlockeds.append([blockediexp,paramkeys[iparam][9:],str(parameters[paramkeys[iparam]])])
    else:
      blockeds = paramkeys[iparam].split(',')
      thisBlock = thisBlock +','.join(blockeds)+','
      thisBlockCoeff = thisBlockCoeff + ','.join([str(parameters[paramkeys[iparam]]) for i in blockeds])+','

  if len(thisAltered) > 0 and thisAltered[-1] == ',':
    thisAltered = thisAltered[0:-1]
    thisAlteredCoeff = thisAlteredCoeff[0:-1]
  
  DataAll = []
  filenames = []
  timenow = time.time()
  timetmp = time.time()

  #Run the three experiments, only change the RGS concentration according to the given model parameters (now saved in blockediexp and blockeds)
  nplotted_experiment = 0
  for iexperiment_type in range(0,3):
    Blocked = thisBlock+','.join([experiments_blockeds[iexperiment_type][i][0] for i in range(0,len(experiments_blockeds[iexperiment_type]))])
    BlockedCoeff = thisBlockCoeff+','.join([str(experiments_blockeds[iexperiment_type][i][1]) for i in range(0,len(experiments_blockeds[iexperiment_type]))])
    Altered = thisAltered + experiments_altereds[iexperiment_type][0]
    AlteredCoeff = thisAlteredCoeff + experiments_altereds[iexperiment_type][1]
    for iexperimentwiseBlock in range(0,len(experimentwiseBlockeds)):
      if experimentwiseBlockeds[iexperimentwiseBlock][0] == iexperiment_type:
        blockeds = experimentwiseBlockeds[iexperimentwiseBlock][1].split(',')
        Blocked = Blocked + ','+experimentwiseBlockeds[iexperimentwiseBlock][1]
        BlockedCoeff = BlockedCoeff + ','+','.join([str(experimentwiseBlockeds[iexperimentwiseBlock][2]) for x in blockeds])
    if len(Altered) > 0 and Altered[0] == ',':
      Altered = Altered[1:]
      AlteredCoeff = AlteredCoeff[1:]

    #Run the simulation. Use a short file name described in randomString to save the file, and delete the file afterwards
    randomString = ''.join(random.choice('ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789') for _ in range(0,10))+'_'+str(rankID)
    print('python3 model_nrn_extfilename.py '+str(Duration)+' '+str(tolerance)+' '+str(gaba_input_onset)+' '+str(gaba_input_N)+' '+str(gaba_input_freq)+' '+str(gaba_input_dur)+' '+
       str(gaba_flux)+' '+str(Ntrains)+' '+str(trainT)+' None '+Blocked+' '+BlockedCoeff+' '+Altered+' '+AlteredCoeff+' fit'+randomString+'.mat '+str(record_T))
    os.system('python3 model_nrn_extfilename.py '+str(Duration)+' '+str(tolerance)+' '+str(gaba_input_onset)+' '+str(gaba_input_N)+' '+str(gaba_input_freq)+' '+str(gaba_input_dur)+' '+
         str(gaba_flux)+' '+str(Ntrains)+' '+str(trainT)+' None '+Blocked+' '+BlockedCoeff+' '+Altered+' '+AlteredCoeff+' fit'+randomString+'.mat '+str(record_T))
    print('Exp.  '+str(iexperiment_type)+', ID='+str(rankID)+' done in '+str(time.time()-timenow)+' sec')

    timetmp = time.time()
    filename = 'fit'+randomString+'.mat'
    filenames.append(filename)
    if not exists(filename):
      print('Error: filename = '+filename+' does not exist, Exp. ='+str(iexperiment))
      DataAll.append([])
      continue

    A = scipy.io.loadmat(filename)
    DataAll.append(A['DATA'])
    times = A['DATA'][0]
    headers = A['headers']    
  if deleteFiles:
    for filename in filenames:
      os.system('rm '+filename)
  return [DataAll,headers]

#Plot the results
f,axarr = subplots(3,1)
for iax in range(0,3):
  axarr[iax].set_position([0.14,0.5-0.21*iax,0.84,0.16])
for iax in range(0,len(axarr)):
  for tick in axarr[iax].xaxis.get_major_ticks() + axarr[iax].yaxis.get_major_ticks():
    tick.label.set_fontsize(4)

#These parameters were obtained from a multi-objective optimisation algorithm (emoo.py, Bahl et al. 2012)
paramTxt = '99.210,706.796,145.360,5.214,20.302,4.551,2.016,3.098,642.686'
paramTxt = '111.531,708.687,165.616,10.797,93.218,2.903,1.728,2.235,2797.567'
paramdict = {'k[0]': float(paramTxt.split(',')[0]), 'k[1,9]': float(paramTxt.split(',')[1]), 'k[3,4,5,6,7,8]': float(paramTxt.split(',')[2]), 'k[15,17,19,21,23]': float(paramTxt.split(',')[3]), 'k[16,18,20,22,24]': float(paramTxt.split(',')[4]), 'OnlyExp0_RGS': float(paramTxt.split(',')[5]), 'OnlyExp1_RGS': float(paramTxt.split(',')[6]), 'OnlyExp2_RGS': float(paramTxt.split(',')[7]), 'gaba_flux': float(paramTxt.split(',')[8])}

#Run the model:
print('paramdict = '+str(paramdict))
A = run_model(paramdict,True,0)

DataAll = A[0]
headers = A[1]
timesAll = [DataAll[i][0] for i in range(0,len(DataAll))]
VGCCboundsAll = [0*x for x in timesAll]
for i in range(0,len(headers)):
  if headers[i].find(' ') != -1:
    headers[i] = headers[i][0:headers[i].find(' ')]
  if 'VGCCGibg' in headers[i]:
    VGCCboundsAll = [DataAll[iexp][i] for iexp in range(0,len(timesAll))]
iGIRK = [i for i in range(0,len(headers)) if 'GIRKGibg' in headers[i]]
nGibgs = [1 if headers[i][-1] == 'g' else int(headers[i][-1]) for i in iGIRK]
coeffs = [0.0, 0.01, 0.06, 0.26, 1.0]
GIRKcondsAll = [0*x for x in timesAll]
GIRK1boundsAll = [0*x for x in timesAll]
for iiGIRK in range(0,len(iGIRK)):
  for iexp in range(0,len(timesAll)):
    GIRKcondsAll[iexp] = GIRKcondsAll[iexp] + coeffs[nGibgs[iiGIRK]]*DataAll[iexp][iGIRK[iiGIRK]]
    if iiGIRK == 0:
      GIRK1boundsAll[iexp] = GIRK1boundsAll[iexp] + DataAll[iexp][iGIRK[0]]

GIRKcond_interpolated = mytools.interpolate(timesAll[0],GIRKcondsAll[0]/max(GIRKcondsAll[0]),[1e6+x for x in data_postsyn_xdata])
VGCCbound_ex_interpolated = mytools.interpolate(timesAll[1],VGCCboundsAll[1]/max(VGCCboundsAll[1]),[1e6+x for x in data_exppr_xdata])
VGCCbound_in_interpolated = mytools.interpolate(timesAll[2],VGCCboundsAll[2]/max(VGCCboundsAll[2]),[1e6+x for x in data_inppr_xdata])

axarr[0].plot(timesAll[0]-1e6,GIRKcondsAll[0]/max(GIRKcondsAll[0]),'b-')
axarr[0].plot(data_postsyn_xdata,GIRKcond_interpolated,'b.')
axarr[0].plot(data_postsyn_xdata,data_postsyn_ydata_scaled,'r.')
axarr[1].plot(timesAll[1]-1e6,VGCCboundsAll[1]/max(VGCCboundsAll[1]),'b-')
axarr[1].plot(data_exppr_xdata,VGCCbound_ex_interpolated,'b.')
axarr[1].plot(data_exppr_xdata,[x*max(VGCCboundsAll[1]/max(VGCCboundsAll[1]))/max(data_antiexppr_ydata_scaled) for x in data_antiexppr_ydata_scaled],'r.')
axarr[2].plot(timesAll[2]-1e6,VGCCboundsAll[2]/max(VGCCboundsAll[2]),'b-')
axarr[2].plot(data_inppr_xdata,VGCCbound_in_interpolated,'b.')
axarr[2].plot(data_inppr_xdata,[x*max(VGCCboundsAll[2]/max(VGCCboundsAll[2]))/max(data_antiinppr_ydata_scaled) for x in data_antiinppr_ydata_scaled],'r.')

for iax in range(0,len(axarr)):
  pos = axarr[iax].get_position()
  f.text(pos.x0 - 0.03 - 0.02, pos.y1 - 0.01+0.02*(iax<2), chr(ord('A')+iax), fontsize=12)

axarr[0].set_ylabel('Normalized\n GIRK cond.', fontsize=7)
axarr[1].set_ylabel('Normalized\n VGCC binding\n(excitatory)', fontsize=7)
axarr[2].set_ylabel('Normalized\n VGCC binding\n(inhibitory)', fontsize=7)
axarr[2].set_xlabel('time (ms)',fontsize=7)

axarr[0].set_xlim([0,4500])
axarr[1].set_xlim([0,4500])
axarr[1].set_xlim([0,4500])

f.savefig('fig_fit.eps')

#Plot each species in three extra figures, one for each simulation
ylims_dict = {'GABABR': 0.0004, 'Gi': 0.0026, 'RGS': 0.001*max(paramdict['OnlyExp0_RGS'],paramdict['OnlyExp1_RGS'],paramdict['OnlyExp2_RGS']), 'GIRK': 0.001, 'VGCC': 0.0001}
ylims_dict_keys = list(ylims_dict.keys())
for ifig in range(0,3):
  f,axs = subplots(5,4)
  for iay in range(0,5):
    for iax in range(0,4):
      axs[iay,iax].set_position([0.08+0.24*iax,0.8-0.1825*iay,0.18,0.115])
      for tick in axs[iay,iax].xaxis.get_major_ticks() + axs[iay,iax].yaxis.get_major_ticks():
        tick.label.set_fontsize(4)
  axarr = sum([axs[i].tolist() for i in range(0,len(axs))]+[[]])
  for iax in range(0,len(axarr)):
    axarr[iax].tick_params(axis='y',direction='out',length=0.8,width=0.3)
    axarr[iax].tick_params(axis='x',direction='out',length=0.8,width=0.3)
    axarr[iax].tick_params(axis='both', which='major', pad=0.12)
    for axis in ['top','bottom','left','right']:
      axarr[iax].spines[axis].set_linewidth(0.3)
  
  for iax in range(0,20):
    axarr[iax].plot((DataAll[ifig][0]-1000000)/1000, DataAll[ifig][1+iax],lw=0.3)
    axarr[iax].set_title(headers[1+iax],fontsize=5)
    axarr[iax].set_xlabel('time (s)',fontsize=5)
    axarr[iax].set_ylabel('conc. (mM)',fontsize=5)
    thisYlim = -1
    for ikey in range(0,len(ylims_dict_keys)):
      if ylims_dict_keys[ikey] in headers[1+iax]:
        if thisYlim == -1:
          thisYlim = ylims_dict[ylims_dict_keys[ikey]]
        else:
          thisYlim = min(thisYlim,ylims_dict[ylims_dict_keys[ikey]])
    if thisYlim > 0:
      axarr[iax].set_ylim([0,thisYlim*1.05])
    axarr[iax].set_xlim([0,10])
  if ifig == 0:
    f.suptitle('Postsynaptic, no VGCCs',fontsize=8)
  if ifig == 1:
    f.suptitle('Presynaptic, EE synapse, no GIRK channels',fontsize=8)
  if ifig == 2:
    f.suptitle('Presynaptic, IE synapse, no GIRK channels',fontsize=8)

  f.savefig('fig_fit_allspecies_syn'+str(ifig)+'.eps')

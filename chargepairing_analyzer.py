import numpy as np
import math
import random
from math import log10, floor
import os
import string
import subprocess
from scipy.optimize import curve_fit
from functools import partial
from lheinfo import get_xsec_witherror
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MaxNLocator
import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from scipy import stats
from scipy.interpolate import interp1d
from scipy.optimize import fsolve, brentq
import matplotlib.lines as mlines
import threading
from threading import Thread
import time
import logging
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing
from joblib import Parallel, delayed

# xgboost stuff
from xgboost_root_varfiles_module import *


############################
# LOCATIONS AND PARAMETERS #
############################
print('HEFT Higgs XSEC Fitting and Analysis -- Charge Pairing version')

#################
# RUN FLAGS:
#################

# read the fit or perform it and write it?
DoFit = True

# Rerun?
ReRunHerwig = False

# Run herwig/analysis in the first place?
RunHerwig = False

# Rerun Analysis?
ReRunAnalysis = False

# Rerun XGBOOST Analysis?
ReRunAnalysisXGBOOST = False

# do the training and write it or not?
DoTraining = False

###############
# PARAMS
###############

# choose the model - HEFT2/HEFT3 or C3-D4 model only 
#MODEL = 'HEFT2' # HEFT2 or C3D4ONLY or HEFT3
#MODEL = 'HEFT3'
MODEL = 'C3D4ONLY'
#MODEL = 'HEFT4'
#MODEL = 'HEFT6'
#MODEL = 'HEFT4C3D4'
print('MODEL=', MODEL)

# choose the type of smearing: NONE, ATLAS, CMS
#SMEARING = 'NONE'
#SMEARING = 'ATLAS'
SMEARING = 'CMS'

# Systematics
Systematics = 0.0 # the alpha value for the systematics
# b-tagging rate
btagging = 0.85

# energy and luminosity
Energy = 13.6 # energy
Luminosity = 3000 # integrated luminosity in /fb to calculate signif
#Energy = 100
#Luminosity = 20000

# ENERGY RESCALING HERE:
DoRescaling = False
EnergyToRescale = 10000
ERESCALE = 1 # not a switch
RESCALETAG = ''
if DoRescaling is True:
    ERESCALE = EnergyToRescale**2/100**2
    RESCALETAG = '_RescaleE' + str(EnergyToRescale)


# K-factors for signal and backgrounds
KFAC_SIGNAL = 2.0
KFAC_BACKGROUNDS = 2.0

# change the KFACTOR ON THE BACKGROUND
CHANGE_KFAC = False
KFACTAG = ''
KFAC_BACKGROUNDS_NEW = 3.0
if CHANGE_KFAC is True:
    KFAC_BACKGROUNDS = KFAC_BACKGROUNDS_NEW
    KFACTAG = '_KFACBKG' + str(KFAC_BACKGROUNDS)

# change the KFACTOR ON THE SIGNL
CHANGE_KFAC_SIGNAL = True
KFACTAG = ''
KFAC_SIGNAL_NEW = 2.5
if CHANGE_KFAC_SIGNAL is True:
    KFAC_SIGNAL = KFAC_SIGNAL_NEW
    KFACTAG = '_KFACSIG' + str(KFAC_SIGNAL)


###############
# END OF PARAMS
###############

    
# array of variables:
variables = {}
variables[0] = 'c3'
variables[1] = 'ct2'
variables[2] = 'ct3'
variables[3] = 'd4'
variables[4] = 'ct1'
variables_latex = {}
variables_latex[0] = 'c_3'
variables_latex[1] = 'c_{t2}'
variables_latex[2] = 'c_{t3}'
variables_latex[3] = 'd_4'
variables_latex[4] = 'c_{t1}'

# constraints on these (fractional):
constraints = {}
constraints[100] = {}
constraints[100][0] = 5/100
constraints[100][1] = 0.1
constraints[100][2] = -1
constraints[100][3] = -1
constraints[13.6] = {}
constraints[13.6][0] = 50/100
constraints[13.6][1] = -1
constraints[13.6][2] = -1
constraints[13.6][3] = -1

# Input file templates for LO, MC@NLO and FxFx:
# the real files have a .in.template extension
HW_template = ['','', '']
HW_template[0] = 'Templates/HW-LO.in' # 0th element is LO

# The reduction factor of the number of events between the LHE file and the actual HW run for each process:
Reduction_Fac = [ '', '', '' ]
Reduction_Fac[0] = 0.999

# Branching ratios:
BR_z_ellell = 3.3632E-2 #  Z -> lepton lepton (one flavour)
BR_w_ellnu = 10.86E-2 # W -> lepton+neutrino (one flavour)
BR_z_vv = 0.2 # Z -> neutrino neutrino (all flavours)
BR_z_qq = 0.116 + 0.156 + 0.1203 + 0.1512 # Z -> qq
BR_z_bb = 0.150998
BR_h_bb = 0.5824
BR_h_gamgam = 0.00229

# chi-sq values in 2D for one and two sigma:
onesigma = 2.278868566376729
twosigma = 5.99

# debug flag
debug = True

# define the process under investigation:
Process = 'gg_hhh'

# the number of runs and tests for fitting
Nruns = 205

# The number of free coefficients to fit in the ME for each process
NCoeffs = {}
if MODEL == 'HEFT2':
    NCoeffs['gg_hhh'] = 18
elif MODEL == 'HEFT3':
    NCoeffs['gg_hhh'] = 8
elif MODEL == 'C3D4ONLY' or MODEL == 'HEFT4C3D4':
    NCoeffs['gg_hhh'] = 9
elif MODEL == 'HEFT4':
    NCoeffs['gg_hhh'] = 25
elif MODEL == 'HEFT6':
    NCoeffs['gg_hhh'] = 80
    
# directory for plots:
plot_dir = 'plots/'

# directory for fits:
fit_dir = 'fits/'

# Dictionaries to hold the fit coefficients and their covariance:
popt = {}
pcov = {}

# Directory for the pickle results
ResultsDir = '/mnt/hdd/Projects/ChargePairing/PickleResults/'

# Constraints directory
ConstraintsDir = 'Constraints/'

# MG5_aMC sub-dir:
if MODEL == 'HEFT2':
    MGLocation = '/home/apapaefs/Projects/GlobalHHH100/MG5_aMC_v2_9_22/' # hhh with 2 insertions in the HEFT
elif MODEL == 'C3D4ONLY' or MODEL == 'HEFT3' or MODEL == "HEFT4":
    MGLocation = '/home/apapaefs/Projects/GlobalHHH100/MG5_aMC_v2_9_24/' # hhh with 2 insertions in the HEFT
elif MODEL == 'HEFT4C3D4' or MODEL == "HEFT6":
    MGLocation = '/home/apapaefs/Projects/GlobalHHH100/MG5_aMC_v2_9_26/' # hhh with 2 insertions in the HEFT

# Analysis executable:
ExecutableSmear = {}
#ExecutableSmear[100] = 'Code/HwSimPostAnalysis_smear_100_example' # to be replaced with the full analysis including smearing
smearing_tag = ''
if SMEARING == 'NONE':
    ExecutableSmear[100] = 'Code/HwSimPostAnalysis_smear_100_variables'
    smearing_tag = ''
elif SMEARING == 'ATLAS':
    ExecutableSmear[100] = 'Code/HwSimPostAnalysis_smear_100_variables_ATLAS'
    smearing_tag = 'ATLAS'
elif SMEARING == 'CMS':
    ExecutableSmear[100] = 'Code/HwSimPostAnalysis_smear_100_variables_CMS'
    ExecutableSmear[13.6] = 'Code/HwSimPostAnalysis_smear_100_variables_CMS'
    smearing_tag = 'CMS' 



# the MG5 subdirectory for each process
ProcLocations = {}
if MODEL == 'HEFT2':
    ProcLocations['gg_hhh'] = 'gg_hhh_mheft2l2_restricted/' # hhh with squared truncation
elif MODEL == 'HEFT3':
    ProcLocations['gg_hhh'] = 'gg_hhh_mheft2l3_morerestricted/' # hhh with cubic truncation
elif MODEL == 'C3D4ONLY': 
    ProcLocations['gg_hhh'] = 'gg_hhh_c3d4/' # hhh with 2 insertions in the HEFT (no truncation)
elif MODEL == 'HEFT4C3D4':
    ProcLocations['gg_hhh'] = 'gg_hhh_restricted5new_heft4/' # hhh with 2 insertions in the HEFT (no truncation)
    #ProcLocations['gg_hhh'] = 'gg_hhh_full_mheft4/'
elif MODEL == 'HEFT4': 
    ProcLocations['gg_hhh'] = 'gg_hhh_restricted_mheft4/' # hhh with 2 insertions in the HEFT (no truncation)
elif MODEL == 'HEFT6': 
    ProcLocations['gg_hhh'] = 'gg_hhh_restricted5_mheft6_new/' # hhh with 2 insertions in the HEFT (no truncation)



# The numbering tag for the run:
if MODEL == 'HEFT2':
    RunNum = '11' # 100 TeV event generation # NEW FOR GLOBAL HHH - HEFT
elif (MODEL == 'C3D4ONLY' and Energy==100) or (MODEL == 'HEFT4C3D4' and Energy==13.6): # C3D4ONLY was 100 TeV, HEFT4C3D4 is 13.6 TeV
    RunNum = '10' # 100 TeV event generation # NEW FOR GLOBAL HHH - C3-D4 MODEL ONLY
elif MODEL == 'HEFT3':
    RunNum = '12'
elif MODEL == 'HEFT4': # 13.6 event generation - HEFT4 restricted (c3,d4,ct2,ct3)
    RunNum = '13'
elif MODEL == 'HEFT6': # 13.6 event generation - HEFT6 restricted (c3,d4,ct2,ct3,ct1)
    RunNum = '14'
elif MODEL == 'C3D4ONLY' and Energy==13.6:
    RunNum = '15' # 13.6 TeV event generation # NEW FOR GLOBAL HHH - C3-D4 MODEL ONLY


# SELECT FINAL STATE HERE:
FinalState = '6b'
if FinalState == '6b':
    FinalState6b = ''
    FinalStatebtau = '#'
    FinalStatebgamma = '#'

# Background Location:
BackgroundLocation = 'Backgrounds/events/'
Backgrounds = []
Backgrounds.append('all_events_6b')
Backgrounds.append('pp_zbbbb')
Backgrounds.append('pp_zzbb')
Backgrounds.append('pp_hzbb')
Backgrounds.append('pp_hhbb')
Backgrounds.append('pp_hbbbb')
Backgrounds.append('pp_hzz')
Backgrounds.append('pp_hhz')
Backgrounds.append('pp_zzz')
Backgrounds.append('gg_hzz')
Backgrounds.append('gg_zzz')
Backgrounds.append('gg_hhz')
Backgrounds_xsec = {}

Backgrounds_xsec[(100, 'all_events_6b')] = 28.328254252903694E3 # cross section for 6b background in fb (100 TeV)
Backgrounds_xsec[(100, 'pp_zbbbb')] = 958.3291282 # cross section for zbbbb background in fb (100 TeV) # 
Backgrounds_xsec[(100, 'pp_zzbb')] = 30.18859257 # cross section for pp_zzbb background in fb (100 TeV) # 
Backgrounds_xsec[(100, 'pp_hzbb')] = 5.417507336 # cross section for pp_hzbb background in fb (100 TeV) #
Backgrounds_xsec[(100, 'pp_zzz')] = 0.4773830264  # cross section for gg_zzz background in fb (100 TeV) # 
Backgrounds_xsec[(100, 'pp_hzz')] = 0.392990544 # cross section for pp_hzz background in fb (100 TeV)
Backgrounds_xsec[(100, 'pp_hhz')] = 0.2149781325 # cross section for pp_hhbb background in fb (100 TeV) # 
Backgrounds_xsec[(100, 'pp_hhbb')] = 0.04761220149 # cross section for pp_hhbb background in fb (100 TeV) #
Backgrounds_xsec[(100, 'pp_hbbbb')] = 1.92239859 # cross section for pp_hbbbb background in fb (100 TeV) # 
Backgrounds_xsec[(100, 'gg_hzz')] = 0.09506002389 # cross section for gg_hzz background in fb (100 TeV) #
Backgrounds_xsec[(100, 'gg_zzz')] = 0.01372856589  # cross section for gg_zzz background in fb (100 TeV) # 
Backgrounds_xsec[(100, 'gg_hhz')] = 0.1700475286  # cross section for gg_hhz background in fb (100 TeV) #
# initial total weight of events (before the analysis that created the _var.root files):
initial_S_SM = 100000
initial_S = 9990

if Energy == 100:
    #xsS=0.0028783E3 # signal cross section at 100 TeV in fb (SM)
    xsS=0.0028783
elif Energy == 13.6:
    xsS = 5.7563e-05 # signal cross section at 13.6 TeV in PB (SM)

signal_SM_file = './Herwig/events/HW-8_SM_6b_var.smear' + smearing_tag + '.root'

# location of the _var root files for the backgrounds:
Background_files = {}
Background_files[(100, 'all_events_6b')] = './Backgrounds/events/HW-all_events_6b_100_var.smear' + smearing_tag + '.root'
Background_files[(100, 'pp_zbbbb')] = './Backgrounds/events/HW-pp_zbbbb_100_var.smear' + smearing_tag + '.root'
Background_files[(100, 'pp_zzbb')] = './Backgrounds/events/HW-pp_zzbb_100_var.smear' + smearing_tag + '.root'
Background_files[(100, 'pp_hzbb')] = './Backgrounds/events/HW-pp_hzbb_100_var.smear' + smearing_tag + '.root'
Background_files[(100, 'pp_hhbb')] = './Backgrounds/events/HW-pp_hhbb_100_var.smear' + smearing_tag + '.root'
Background_files[(100, 'pp_hbbbb')] = './Backgrounds/events/HW-pp_hbbbb_100_var.smear' + smearing_tag + '.root'
Background_files[(100, 'pp_hzz')] = './Backgrounds/events/HW-pp_hzz_100_var.smear' + smearing_tag + '.root'
Background_files[(100, 'pp_zzz')] = './Backgrounds/events/HW-pp_zzz_100_var.smear' + smearing_tag + '.root'
Background_files[(100, 'pp_hhz')] = './Backgrounds/events/HW-pp_hhz_100_var.smear' + smearing_tag + '.root'
Background_files[(100, 'gg_hzz')] = './Backgrounds/events/HW-gg_hzz_100_var.smear' + smearing_tag + '.root'
Background_files[(100, 'gg_zzz')] = './Backgrounds/events/HW-gg_zzz_100_var.smear' + smearing_tag + '.root'
Background_files[(100, 'gg_hhz')] = './Backgrounds/events/HW-gg_hhz_100_var.smear' + smearing_tag + '.root'

# initial weight of Monte Carlo events (at the start of the analysis that generated the var root files):
initial_B = {}
initial_B['all_events_6b'] = 864960
initial_B['pp_zbbbb'] = 200000
initial_B['pp_zzbb'] = 200000
initial_B['pp_hzbb'] = 200000
initial_B['pp_hzz'] = 200000
initial_B['pp_zzz'] = 200000
initial_B['pp_hhz'] = 200000
initial_B['pp_hhbb'] = 200000
initial_B['pp_hbbbb'] = 200000
initial_B['gg_zzz'] = 100000
initial_B['gg_hzz'] = 200000
initial_B['gg_hhz'] = 200000

# initial actual (i.e. at luminosity) number of events for backgrounds
initial_NB = {}

# background ids:
idB = {}
idB['all_events_6b'] = 1
idB['pp_zbbbb'] = 2
idB['pp_zzbb'] = 3
idB['pp_hzbb'] = 4
idB['pp_hhz'] = 5
idB['pp_hzz'] = 6
idB['pp_zzz'] = 7
idB['pp_hhbb'] = 8
idB['pp_hbbbb'] = 9
idB['gg_hzz'] = 10
idB['gg_zzz'] = 11
idB['gg_hhz'] = 12
    


# factors to apply to signal and background (K-factors and BRs)
sig_factors = KFAC_SIGNAL * BR_h_bb**3 * btagging**6 * ERESCALE
bkg_factors = KFAC_BACKGROUNDS * btagging**6 * ERESCALE # BRs already applied. The k-factor is uniform



# Herwig input file sub-dir and output for the events
HerwigLocation = 'Herwig/'
HerwigOutputLocation = HerwigLocation + 'events/'
HerwigOutputDirectory = HerwigOutputLocation



#########################################################
# FUNCTIONS                                             # 
#########################################################

# function to get template
def getTemplate(basename):
    with open('%s.template' % basename, 'r') as f:
        templateText = f.read()
    return string.Template( templateText )

# write a filename
def writeFile(filename, text):
    with open(filename,'w') as f:
        f.write(text)

# round to a certain number of significant figures
def round_sig(x, sig=4):
    if x == 0.:
        return 0.
    if math.isnan(x) is True:
        print('Warning, NaN!')
        return 0.
    return round(x, sig-int(floor(log10(abs(x))))-1)

# gaussian function
def gaussian(x, mu, delta):
    return 1./(np.sqrt(2.*np.pi)*delta)*np.exp(-np.power((x - mu)/delta, 2.)/2)

# function for Higgs boson triple production in the HEFT:
# only c3, d4, ct2, ct3 are assumed to be relevant
def func_CX(couplings=[], *coeffs, procname):
    #print('couplings=', couplings)
    Msq = 0
    if procname == 'gg_hhh':
        if MODEL == 'HEFT2':
            S1, S2, A1, A2, B1, B2, C1, C2, D1, D2, E1, E2, F1, F2, L1, L2, N1, N2 = [float(coef) for coef in coeffs]
            c3, d4, cg1, cg2, ct1, cb1, ct2, cb2, ct3, cb3 = couplings
            Msq = A1**2*c3**2 + 2*A1*B1*c3*d4 + 2*A1*L1*ct2*c3 + 2*A1*N1*ct3*c3 + 2*A1*S1*c3 + A2**2*c3**2 + 2*A2*B2*c3*d4 + 2*A2*L2*ct2*c3 + 2*A2*N2*ct3*c3 + 2*A2*S2*c3 + B1**2*d4**2 + 2*B1*L1*ct2*d4 + 2*B1*N1*ct3*d4 + 2*B1*S1*d4 + B2**2*d4**2 + 2*B2*L2*ct2*d4 + 2*B2*N2*ct3*d4 + 2*B2*S2*d4 + 2*C1*S1*c3**2 + 2*C2*S2*c3**2 + 2*D1*S1*d4**2 + 2*D2*S2*d4**2 + 2*E1*S1*ct2**2 + 2*E2*S2*ct2**2 + 2*F1*S1*ct3**2 + 2*F2*S2*ct3**2 + L1**2*ct2**2 + 2*L1*N1*ct2*ct3 + 2*L1*S1*ct2 + L2**2*ct2**2 + 2*L2*N2*ct2*ct3 + 2*L2*S2*ct2 + N1**2*ct3**2 + 2*N1*S1*ct3 + N2**2*ct3**2 + 2*N2*S2*ct3 + S1**2 + S2**2
        elif MODEL == 'HEFT3':
            S1, B1, C1, D1, E1, F1, L1, N1 = [float(coef) for coef in coeffs]
            c3, d4, cg1, cg2, ct1, cb1, ct2, cb2, ct3, cb3 = couplings
            Msq = S1 + B1 * c3**3 + C1 * c3**2 * d4 + D1 * c3**2 + E1 * d4**2 + F1 * c3 * d4 + L1 * d4 + N1 * c3
        elif MODEL == 'C3D4ONLY' or MODEL == 'HEFT4C3D4': 
            S1, A1, B1, C1, D1, E1, F1, L1, N1 = [float(coef) for coef in coeffs]
            c3, d4, cg1, cg2, ct1, cb1, ct2, cb2, ct3, cb3 = couplings
            Msq = S1 + A1 * c3**4 + B1 * c3**3 + C1 * c3**2 * d4 + D1 * c3**2 + E1 * d4**2 + F1 * c3 * d4 + L1 * d4 + N1 * c3
        elif MODEL == 'HEFT4':
            A, B, C, D, E, F, G, H, J, K, L, M, N, O, P, Q, R, S, T, W, X, Y, Z, ZZ, WW= [float(coef) for coef in coeffs]
            c3, d4, cg1, cg2, ct1, cb1, ct2, cb2, ct3, cb3 = couplings
            Msq = A*c3**2*d4 + B*c3**2*ct2**2 + C*c3**2*ct2 + D * c3**2*ct3 + E*c3**2 + F*c3**4 + G*c3**3*ct2 + H*c3**3 + J*c3*ct2*d4 + K*c3*d4 + L*c3*ct2*ct3 + M*c3*ct2 + N*c3*ct2**2 + O*c3*ct3 + P*c3 + Q*d4**2 + R*ct2*d4 + S*ct3*d4 + T*d4 + W*ct2**2 + X*ct2*ct3 + Y*ct2 + Z*ct3**2 + ZZ*ct3 + WW
        elif MODEL == 'HEFT6':
            A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U, V, W, X, Y, Z, AA, AB, AC, AD, AE, AF, AG, AH, AI, AJ, AK, AL, AM, AN, AO, AP, AQ, AR, AS, AT, AU, AV, AW, AX, AY, AZ, BA, BB, BC, BD, BE, BF, BG, BH, BI, BJ, BK, BL, BM, BN, BO, BP, BQ, BR, BS, BT, BU, BV, BW, BX, BY, BZ, CA, CB = [float(coef) for coef in coeffs]
            c3, d4, cg1, cg2, cb3, cb1, ct2, cb2, ct3, ct1 = couplings # notice change of order here
            Msq = A*c3**2*ct1*d4 + B*c3**2*ct1**2*d4 + C*c3**2*d4 + D*c3**2*ct1**4 + E*c3**2*ct1**2*ct2 + F*c3**2*ct1**2 + G*c3**2*ct1*ct2 + H*c3**2*ct1*ct3 + I*c3**2*ct1 + J*c3**2*ct1**3 + K*c3**2*ct2**2 + L*c3**2*ct2 + M*c3**2*ct3 + N*c3**2 + O*c3**4*ct1**2 + P*c3**4*ct1 + Q*c3**4 + R*c3**3*ct1*ct2 + S*c3**3*ct1 + T*c3**3*ct1**2 + U*c3**3*ct1**3 + V*c3**3*ct2 + W*c3**3 + X*c3*ct1*ct2*d4 + Y*c3*ct1*d4 + Z*c3*ct1**2*d4 + AA*c3*ct1**3*d4 + AB*c3*ct2*d4 + AC*c3*d4 + AD*c3*ct1*ct2 + AE*c3*ct1*ct2**2 + AF*c3*ct1*ct3 + AG*c3*ct1 + AH*c3*ct1**2*ct2 + AI*c3*ct1**2*ct3 + AJ*c3*ct1**2 + AK*c3*ct1**3*ct2 + AL*c3*ct1**3 + AM*c3*ct1**4 + AN*c3*ct1**5 + AO*c3*ct2*ct3 + AP*c3*ct2 + AQ*c3*ct2**2 + AR*c3*ct3 + AS*c3 + AT*ct1**2*d4**2 + AU*ct1*d4**2 + AV*d4**2 + AW*ct1*ct2*d4 + AX*ct1*ct3*d4 + AY*ct1*d4 + AZ*ct1**2*ct2*d4 + BA*ct1**2*d4 + BB*ct1**3*d4 + BC*ct1**4*d4 + BD*ct2*d4 + BE*ct3*d4 + BF*d4 + BG*ct1**2*ct2**2 + BH*ct1**2*ct2 + BI*ct1**2*ct3 + BJ*ct1**2 + BK*ct1**4*ct2 + BL*ct1**4 + BM*ct1**6 + BN*ct1**3*ct2 + BO*ct1**3*ct3 + BP*ct1**3 + BQ*ct1*ct2*ct3 + BR*ct1*ct2 + BS*ct1*ct2**2 + BT*ct1*ct3 + BU*ct1 + BV*ct1**5 + BW*ct2**2 + BX*ct2*ct3 + BY*ct2 + BZ*ct3**2 + CA*ct3 + CB
    return Msq


# function for Higgs boson triple production in the HEFT (PLOT VERSION)
def func_t_CX(c3, d4, ct2, ct3, coeffs, procname):
    if procname == 'gg_hhh':
        if MODEL == 'HEFT2':
            S1, S2, A1, A2, B1, B2, C1, C2, D1, D2, E1, E2, F1, F2, L1, L2, N1, N2 = [float(coef) for coef in coeffs]
            Msq = A1**2*c3**2 + 2*A1*B1*c3*d4 + 2*A1*L1*ct2*c3 + 2*A1*N1*ct3*c3 + 2*A1*S1*c3 + A2**2*c3**2 + 2*A2*B2*c3*d4 + 2*A2*L2*ct2*c3 + 2*A2*N2*ct3*c3 + 2*A2*S2*c3 + B1**2*d4**2 + 2*B1*L1*ct2*d4 + 2*B1*N1*ct3*d4 + 2*B1*S1*d4 + B2**2*d4**2 + 2*B2*L2*ct2*d4 + 2*B2*N2*ct3*d4 + 2*B2*S2*d4 + 2*C1*S1*c3**2 + 2*C2*S2*c3**2 + 2*D1*S1*d4**2 + 2*D2*S2*d4**2 + 2*E1*S1*ct2**2 + 2*E2*S2*ct2**2 + 2*F1*S1*ct3**2 + 2*F2*S2*ct3**2 + L1**2*ct2**2 + 2*L1*N1*ct2*ct3 + 2*L1*S1*ct2 + L2**2*ct2**2 + 2*L2*N2*ct2*ct3 + 2*L2*S2*ct2 + N1**2*ct3**2 + 2*N1*S1*ct3 + N2**2*ct3**2 + 2*N2*S2*ct3 + S1**2 + S2**2
        elif  MODEL == 'HEFT3':
            S1, B1, C1, D1, E1, F1, L1, N1 = [float(coef) for coef in coeffs]
            Msq = S1 + B1 * c3**3 + C1 * c3**2 * d4 + D1 * c3**2 + E1 * d4**2 + F1 * c3 * d4 + L1 * d4 + N1 * c3
        elif MODEL == 'C3D4ONLY' or MODEL == 'HEFT4C3D4': 
            S1, A1, B1, C1, D1, E1, F1, L1, N1 = [float(coef) for coef in coeffs]
            Msq = S1 + A1 * c3**4 + B1 * c3**3 + C1 * c3**2 * d4 + D1 * c3**2 + E1 * d4**2 + F1 * c3 * d4 + L1 * d4 + N1 * c3
        elif MODEL == 'HEFT4':
            A, B, C, D, E, F, G, H, J, K, L, M, N, O, P, Q, R, S, T, W, X, Y, Z, ZZ, WW = [float(coef) for coef in coeffs]
            Msq = A*c3**2*d4 + B*c3**2*ct2**2 + C*c3**2*ct2 + D * c3**2*ct3 + E*c3**2 + F*c3**4 + G*c3**3*ct2 + H*c3**3 + J*c3*ct2*d4 + K*c3*d4 + L*c3*ct2*ct3 + M*c3*ct2 + N*c3*ct2**2 + O*c3*ct3 + P*c3 + Q*d4**2 + R*ct2*d4 + S*ct3*d4 + T*d4 + W*ct2**2 + X*ct2*ct3 + Y*ct2 + Z*ct3**2 + ZZ*ct3 + WW
        elif MODEL == 'HEFT6':
            A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U, V, W, X, Y, Z, AA, AB, AC, AD, AE, AF, AG, AH, AI, AJ, AK, AL, AM, AN, AO, AP, AQ, AR, AS, AT, AU, AV, AW, AX, AY, AZ, BA, BB, BC, BD, BE, BF, BG, BH, BI, BJ, BK, BL, BM, BN, BO, BP, BQ, BR, BS, BT, BU, BV, BW, BX, BY, BZ, CA, CB = [float(coef) for coef in coeffs]
            Msq = A*c3**2*ct1*d4 + B*c3**2*ct1**2*d4 + C*c3**2*d4 + D*c3**2*ct1**4 + E*c3**2*ct1**2*ct2 + F*c3**2*ct1**2 + G*c3**2*ct1*ct2 + H*c3**2*ct1*ct3 + I*c3**2*ct1 + J*c3**2*ct1**3 + K*c3**2*ct2**2 + L*c3**2*ct2 + M*c3**2*ct3 + N*c3**2 + O*c3**4*ct1**2 + P*c3**4*ct1 + Q*c3**4 + R*c3**3*ct1*ct2 + S*c3**3*ct1 + T*c3**3*ct1**2 + U*c3**3*ct1**3 + V*c3**3*ct2 + W*c3**3 + X*c3*ct1*ct2*d4 + Y*c3*ct1*d4 + Z*c3*ct1**2*d4 + AA*c3*ct1**3*d4 + AB*c3*ct2*d4 + AC*c3*d4 + AD*c3*ct1*ct2 + AE*c3*ct1*ct2**2 + AF*c3*ct1*ct3 + AG*c3*ct1 + AH*c3*ct1**2*ct2 + AI*c3*ct1**2*ct3 + AJ*c3*ct1**2 + AK*c3*ct1**3*ct2 + AL*c3*ct1**3 + AM*c3*ct1**4 + AN*c3*ct1**5 + AO*c3*ct2*ct3 + AP*c3*ct2 + AQ*c3*ct2**2 + AR*c3*ct3 + AS*c3 + AT*ct1**2*d4**2 + AU*ct1*d4**2 + AV*d4**2 + AW*ct1*ct2*d4 + AX*ct1*ct3*d4 + AY*ct1*d4 + AZ*ct1**2*ct2*d4 + BA*ct1**2*d4 + BB*ct1**3*d4 + BC*ct1**4*d4 + BD*ct2*d4 + BE*ct3*d4 + BF*d4 + BG*ct1**2*ct2**2 + BH*ct1**2*ct2 + BI*ct1**2*ct3 + BJ*ct1**2 + BK*ct1**4*ct2 + BL*ct1**4 + BM*ct1**6 + BN*ct1**3*ct2 + BO*ct1**3*ct3 + BP*ct1**3 + BQ*ct1*ct2*ct3 + BR*ct1*ct2 + BS*ct1*ct2**2 + BT*ct1*ct3 + BU*ct1 + BV*ct1**5 + BW*ct2**2 + BX*ct2*ct3 + BY*ct2 + BZ*ct3**2 + CA*ct3 + CB
    return Msq


# function to read the mg5 event cross sections                
def read_files(runnum, mgloc, procloc, procname, CouplingsArray, nruns):
    X = []
    Z = []
    ZERR = []
    XSEC = {}
    for coups in CouplingsArray:
        #print(coups)
        lhe = 'run_' + procname + '_' + str(runnum) + '_' + '_'.join((coups)) + '/unweighted_events.lhe.gz'
        lhefile = mgloc + '/' + procloc + 'Events/' + lhe
        print('lhefile read=', lhefile)
        #TestBool = True
        #if TestBool is False:
        if os.path.exists(lhefile) is False:
            print('Error, lhe file or summary file:', lhefile, 'does not exist!')
            exit()
        else:
            #zgrepcommand = 'zgrep "Integrated weight" ' + lhefile
            #p = subprocess.Popen(zgrepcommand, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, cwd='.')
            #for line in iter(p.stdout.readline, b''):
            #    xsec = float(line.split()[5])
            #print(coups, xsec)
            xsec, xsecerr = get_xsec_witherror(lhefile)
            print(coups, xsec)
            #xsec = 0
            coups_tuple = []
            for mm in range(len(coups)):
                coups_tuple.append(float(coups[mm]))
            X.append(tuple(coups_tuple))
            Z.append(float(xsec))
            ZERR.append(float(xsecerr))
            XSEC[tuple(coups_tuple)] = float(xsec)
            #print(X)
    return np.transpose(X), Z, ZERR, XSEC

def gen_coupbdasarray_dim_rand_range(coup_min, coup_max, nruns, randseed):
    random.seed(randseed)
    
    CouplingsArray_R = []
    CouplingsArrayF_R = []
    random_choice = 0
    # NOTE: legacy zeroes to comply with previous code! 
    while random_choice < nruns:
        coup1 = coup_min[0] + (coup_max[0] - coup_min[0]) * random.random()
        coup2 = coup_min[1] + (coup_max[1] - coup_min[1]) * random.random()
        coup3 = 0.0 * random.random()
        coup4 = 0.0 * random.random()
        coup5 = 0.0 * random.random()
        coup6 = 0.0 * random.random()
        coup7 = coup_min[2] + (coup_max[2] - coup_min[2]) * random.random()
        coup8 = 0.0 * random.random()
        coup9 = coup_min[3] + (coup_max[3] - coup_min[3]) * random.random()
        if MODEL == 'HEFT6':
            coup10 = coup_min[4] + (coup_max[4] - coup_min[4]) * random.random()
        else:
            coup10 = 0.0 * random.random()
        CouplingsArray = [str(round_sig(coup1,4)), str(round_sig(coup2,4)), str(round_sig(coup3,4)), str(round_sig(coup4,4)), str(round_sig(coup5,4)), str(round_sig(coup6,4)), str(round_sig(coup7,4)), str(round_sig(coup8,4)), str(round_sig(coup9,4)), str(round_sig(coup10,4))]
        CouplingsArrayF = tuple([round_sig(coup1,4), round_sig(coup2,4), round_sig(coup3,4), round_sig(coup4,4), round_sig(coup5,4), round_sig(coup6,4), round_sig(coup7,4), round_sig(coup8,4), round_sig(coup9,4), round_sig(coup10,4)])
        #print('CouplingsArray RANDOM=', CouplingsArray)
        CouplingsArray_R.append(CouplingsArray)
        CouplingsArrayF_R.append(CouplingsArrayF)
        random_choice = random_choice + 1
    print('Generated random arrays for Nruns=', nruns)
    return CouplingsArray_R, CouplingsArrayF_R


# function to read the mg5 event cross sections and compare to the fit              
def test_fit(runnum, mgloc, procloc, procname, CouplingsArray, ntotal, popt):
    X = []
    Z = []
    XSEC = {}
    ZERR = []
    func_CX_proc = partial(func_CX, procname=Process)
    fracdiff_avg = 0.
    for coups in CouplingsArray:
        lhe = 'run_' + procname + '_' + str(runnum) + '_' + '_'.join((coups)) + '/unweighted_events.lhe.gz'
        lhefile = mgloc + '/' + procloc + 'Events/' + lhe
        #TestBool = True
        #if TestBool is False:
        if os.path.exists(lhefile) is False:
            print('Error, lhe file or summary file:', lhefile, 'does not exist!')
            exit()
        else:
            #zgrepcommand = 'zgrep "Integrated weight" ' + lhefile
            #p = subprocess.Popen(zgrepcommand, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, cwd='.')
            #for line in iter(p.stdout.readline, b''):
            #    xsec = float(line.split()[5])
            xsec, xsecerr = get_xsec_witherror(lhefile)
            coups_tuple = []
            for mm in range(len(coups)):
                coups_tuple.append(float(coups[mm]))
            X.append(tuple(coups_tuple))
            Z.append(float(xsec))
            ZERR.append(float(xsecerr))
            # get the fitted XSEC
            xsec_fit = func_CX_proc(coups_tuple, *popt)
            fracdiff = abs(xsec-xsec_fit)/xsec
            if fracdiff > 0.2:
                print(coups, xsec)
                print('!!! lhefile=', lhefile)
                print('!!! xsec: real, fitted, frac diff =', xsec, xsec_fit, fracdiff)
            fracdiff_avg = fracdiff_avg + fracdiff
            XSEC[tuple(coups_tuple)] = float(xsec)
            #print(X)
    print('average fractional difference =', fracdiff_avg/ntotal)
    return np.transpose(X), Z, ZERR, XSEC




# 2D contour plot 
def contour_xsec(procname, plotname, plottitle, fit_coeffs, var1, var2, xlim, ylim, axext='', figext='', smtext=True, starsize=15, setxlabel=True, setylabel=True, nbins=100, savefig=True,variables=variables, variables_latex=variables_latex, labelsize=20, normalbar=True, contours=np.arange(0, 10, 0.5),norm_to_zeroth=True):
    output = procname + '_' + plotname + '_' + var1 + '_' + var2
    print('Plotting', output)
    nvar1 = [key for key, value in variables.items() if value == var1][0]
    nvar2 = [key for key, value in variables.items() if value == var2][0]
    #print(var1, var2)
    #print(nvar1, nvar2)
    # construct the axes for the plot
    # no need to modify this if you just need one plot
    gs = gridspec.GridSpec(4, 4)
    if figext == '':
        fig = plt.figure()
    else:
        fig = figext
    if axext == '':
        ax = fig.add_subplot(111)
    else:
        ax=axext
    ax.grid(False)
    ax.set_title(plottitle)
    # create legend and plot/font size
    #ax.legend()
    #ax.legend(loc="upper right", numpoints=1, frameon=False, prop={'size':8})
    # set the ticks, labels and limits etc.
    xlab = '$' + variables_latex[nvar1] + '$'
    ylab = '$' + variables_latex[nvar2] + '$'
    if setylabel == True:
        ax.set_ylabel(ylab, fontsize=labelsize)
    if setxlabel == True:
        ax.set_xlabel(xlab, fontsize=labelsize)
    
    # choose x and y log scales
    #if ylog:
    #    ax.set_yscale('log')
    #else:
    #    ax.set_yscale('linear')
    #if xlog:
    #    ax.set_xscale('log')
    #else:
    #    ax.set_xscale('linear')
    # set the limits on the x and y axes if required below:
    ymin = ylim[0]
    ymax = ylim[1]
    xmin = xlim[0]
    xmax = xlim[1]
    plt.xlim([xmin,xmax])
    plt.ylim([ymin,ymax])
    ctexts = []
    cvartexts = []
    for i in range(0, len(variables.keys())):
        if i != nvar1 and i != nvar2:
            ctext = variables[i] + '=0'
            ctexts.append(ctext)
        else:
            cvartexts.append(variables[i])
    #print(ctexts)
    fstr = 'partial(func_t_CX, ' + ','.join([ct for ct in ctexts]) + ', procname=Process)'
    global func_CX_partial
    func_CX_partial = eval(fstr)
    #print(func_CX_partial)
    #print(fit_coeffs)
    global fit_coeffs_g
    fit_coeffs_g = fit_coeffs
    #print(cvartexts[0], cvartexts[1])
    if norm_to_zeroth is True:
        feval = 'func_CX_partial(' + cvartexts[0] +'=x1,' + cvartexts[1] + '=x2,coeffs=fit_coeffs_g)/func_CX_partial(' + cvartexts[0] +'=0,' + cvartexts[1] + '=0,coeffs=fit_coeffs_g)'
    else:
        feval = 'func_CX_partial(' + cvartexts[0] +'=x1,' + cvartexts[1] + '=x2,coeffs=fit_coeffs_g)'
    func_fin = lambda x1, x2: eval(feval)
    #print(func_fin(0.05, -0.05))
    x = np.linspace(xlim[0], xlim[1], nbins)
    y = np.linspace(ylim[0], ylim[1], nbins)
    X, Y = np.meshgrid(x,y)
    Z = func_fin(X,Y)
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    cont = ax.contourf(X, Y, Z, contours, cmap='Spectral', extend='max')
    ax.plot(0,0,marker='*',ms=starsize, color='black')
    if smtext == True:
        ax.text(0.53, 0.53,"SM", transform=ax.transAxes)
    if normalbar == True:
        plt.colorbar(cont)
    if savefig == True:
        # save the figure
        print('saving the figure')
        # save the figure in PDF format
        infile = output + '.dat'
        print('---')
        print('output in', infile.replace('.dat','.pdf'))
        plt.savefig(infile.replace('.dat','.pdf'), bbox_inches='tight')
        plt.close(fig)
    return cont

# 1D plot of the xsec
def oned_xsec(procname, plotname, plottitle, fit_coeffs, var1, xlim, ylim, axext='', figext='', smtext=True, starsize=15, setxlabel=True, setylabel=True, nbins=100, savefig=True,variables=variables, variables_latex=variables_latex, labelsize=20, normalbar=True, contours=np.arange(0, 20, 0.5),norm_to_zeroth=True):
    output = procname + '_' + plotname + '_' + var1 
    print('Plotting', output)
    nvar1 = [key for key, value in variables.items() if value == var1][0]
    #print(var1, var2)
    #print(nvar1, nvar2)
    # construct the axes for the plot
    # no need to modify this if you just need one plot
    gs = gridspec.GridSpec(4, 4)
    if figext == '':
        fig = plt.figure()
    else:
        fig = figext
    if axext == '':
        ax = fig.add_subplot(111)
    else:
        ax=axext
    ax.grid(False)
    ax.set_title(plottitle)
    # create legend and plot/font size
    #ax.legend()
    #ax.legend(loc="upper right", numpoints=1, frameon=False, prop={'size':8})
    # set the ticks, labels and limits etc.
    xlab = '$' + variables_latex[nvar1] + '$'
    ylab = r'$\sigma/\sigma_\mathrm{SM}$'
    if setylabel == True:
        ax.set_ylabel(ylab, fontsize=labelsize)
    if setxlabel == True:
        ax.set_xlabel(xlab, fontsize=labelsize)
    
    # choose x and y log scales
    #if ylog:
    #    ax.set_yscale('log')
    #else:
    #    ax.set_yscale('linear')
    #if xlog:
    #    ax.set_xscale('log')
    #else:
    #    ax.set_xscale('linear')
    # set the limits on the x and y axes if required below:
    ymin = ylim[0]
    ymax = ylim[1]
    xmin = xlim[0]
    xmax = xlim[1]
    plt.xlim([xmin,xmax])
    plt.ylim([ymin,ymax])
    ctexts = []
    cvartexts = []
    for i in range(0, len(variables.keys())):
        if i != nvar1:
            ctext = variables[i] + '=0'
            ctexts.append(ctext)
        else:
            cvartexts.append(variables[i])
    #print(ctexts)
    fstr = 'partial(func_t_CX, ' + ','.join([ct for ct in ctexts]) + ', procname=Process)'
    global func_CX_partial
    func_CX_partial = eval(fstr)
    #print(func_CX_partial)
    #print(fit_coeffs)
    global fit_coeffs_g
    fit_coeffs_g = fit_coeffs
    #print(cvartexts[0], cvartexts[1])
    if norm_to_zeroth is True:
        feval = 'func_CX_partial(' + cvartexts[0] +'=x1,coeffs=fit_coeffs_g)/func_CX_partial(' + cvartexts[0] +'=0,coeffs=fit_coeffs_g)'
    else:
        feval = 'func_CX_partial(' + cvartexts[0] +'=x1,coeffs=fit_coeffs_g)'
    func_fin = lambda x1: eval(feval)
    #print(func_fin(0.05, -0.05))
    X = np.linspace(xlim[0], xlim[1], nbins)
    Z = func_fin(X)
    
    line = ax.plot(X, Z, marker='', ls='--', color='blue', lw=3)
    ax.axhline(y=1.0,  linewidth=0.5, color = 'k', ls='--')
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    if savefig == True:
        # save the figure
        print('saving the figure')
        # save the figure in PDF format
        infile = output + '.dat'
        print('---')
        print('output in', plot_dir + infile.replace('.dat','.pdf'))
        plt.savefig(plot_dir + infile.replace('.dat','.pdf'), bbox_inches='tight')
        plt.close(fig)
    return line


def correlation_plot(procname, plotname, popt, varnames,plottitle='',contours=np.arange(-2, 32, 2),norm_to_zeroth=True):
    ###################################################################################
    # correlation plots for cross section
    ###################################################################################
    print('---')
    print('plotting correlation plots for', procname, plotname)
    # plot settings ########
    output = procname + '_' + plotname + '_correlation'

    fig2 = plt.figure(figsize=(9,9))
    spec2 = gridspec.GridSpec(ncols=len(variables), nrows=len(variables),wspace=0, hspace=0, figure=fig2)

    f2_ax_array = []
    cc = 0
    for i in range(len(varnames)):
        for j in range(len(varnames)):
            if i > j:
                if procname == 'gg_hh' and (varnames[i] == 'd4' or varnames[j] == 'd4' or varnames[i] == 'ct3' or varnames[j] == 'ct3' or varnames[i] == 'cb3' or varnames[j] == 'cb3'):
                    continue
                f2_ax = fig2.add_subplot(spec2[i, j])
                f2_ax.set_box_aspect(1)
                f2_ax.xaxis.set_major_locator(MaxNLocator(nbins=4,prune='both'))
                f2_ax.yaxis.set_major_locator(MaxNLocator(nbins=4,prune='both'))
                f2_ax.tick_params(axis='both', labelsize=5)
                f2_ax_array.append(f2_ax)
                cc = cc+1
    spec2.update(wspace=0,hspace=0)

    nplots = len(varnames)**2
    cc = 0
    for i in range(len(varnames)):
        for j in range(len(varnames)):
            if i > j:
                if procname == 'gg_hh' and (varnames[i] == 'd4' or varnames[j] == 'd4' or varnames[i] == 'ct3' or varnames[j] == 'ct3' or varnames[i] == 'cb3' or varnames[j] == 'cb3'):
                    continue
                labelx=False
                labely=False
                if i == len(varnames)-1 or (procname=='gg_hh' and i==len(varnames)-4):
                    labelx=True
                else:
                    f2_ax_array[cc].set(xticks=[])
                if j == 0:
                    labely = True
                else:
                    f2_ax_array[cc].set(yticks=[])
                if varnames[j] == 'c3' and varnames[i] != 'd4':
                     cont = contour_xsec(Process, 'xsec', '', popt, varnames[j], varnames[i], [-10.0, 10.0], [-1.0, 1.0], smtext=False, starsize=2, setxlabel=labelx, setylabel=labely, figext=fig2, axext=f2_ax_array[cc], savefig=False, labelsize=15, normalbar=False,contours=contours, norm_to_zeroth=norm_to_zeroth)
                elif varnames[j] == 'c3' and varnames[i] == 'd4':
                    cont = contour_xsec(Process, 'xsec', '', popt, varnames[j], varnames[i], [-10.0, 10.0], [-40.0, 40.0], smtext=False, starsize=2, setxlabel=labelx, setylabel=labely, figext=fig2, axext=f2_ax_array[cc], savefig=False, labelsize=15, normalbar=False,contours=contours, norm_to_zeroth=norm_to_zeroth)
                elif varnames[j] != 'c3' and varnames[i] == 'd4':
                    cont = contour_xsec(Process, 'xsec', '', popt, varnames[j], varnames[i], [-1.0, 1.0], [-40.0, 40.0], smtext=False, starsize=2, setxlabel=labelx, setylabel=labely, figext=fig2, axext=f2_ax_array[cc], savefig=False, labelsize=15, normalbar=False,contours=contours, norm_to_zeroth=norm_to_zeroth)
                else:
                    cont = contour_xsec(Process, 'xsec', '', popt, varnames[j], varnames[i], [-1.0, 1.0], [-1.0, 1.0], smtext=False, starsize=2, setxlabel=labelx, setylabel=labely, figext=fig2, axext=f2_ax_array[cc], savefig=False, labelsize=15, normalbar=False,contours=contours, norm_to_zeroth=norm_to_zeroth)
                cc = cc + 1
    #fig2.tight_layout()
    #plt.subplots_adjust(wspace=0, hspace=0)
    #fig2.colorbar(cont, ax=f2_ax_array[-1])
    if procname == 'gg_hhh':
        axins = inset_axes(f2_ax_array[-1], # here using axis of the lowest plot
                width="20%",  # width = 5% of parent_bbox width
                height="280%",  # height : 340% good for a (4x4) Grid
                loc='lower left',
                    bbox_to_anchor=(1.08, 0.15, 1, 1),
                    bbox_transform=f2_ax_array[-1].transAxes,
                borderpad=0,
                )
    elif procname == 'gg_hh': 
        axins = inset_axes(f2_ax_array[-1], # here using axis of the lowest plot
                width="28%",  # width = 5% of parent_bbox width
                height="550%",  # height : 340% good for a (4x4) Grid
                loc='lower left',
                    bbox_to_anchor=(1.04, 0.1, 1, 1),
                    bbox_transform=f2_ax_array[-1].transAxes,
                borderpad=0,
                )
        
    cb = fig2.colorbar(cont, cax=axins)
    if procname == 'gg_hhh':
        fig2.suptitle(plottitle,y=0.72,fontsize=15)
    elif procname == 'gg_hh':
        fig2.suptitle(plottitle,x=0.4,y=0.8,fontsize=10)
    # save the figure
    print('saving the figure')
    # save the figure in PDF format
    infile = output + '.dat'
    print('---')
    print('output in', plot_dir + infile.replace('.dat','.pdf'))
    plt.savefig(plot_dir + infile.replace('.dat','.pdf'), bbox_inches='tight')
    plt.close(fig2)

    ####################


# function to save the fit for Process in the fit_dir for a specific RunNum:
def saveFit(popt, pcov, Process, RunNum):
    filename = fit_dir + 'fit_' + Process + '_run' + str(RunNum) + smearing_tag + '.dat'
    f = open(filename,'w')
    f.write('\t'.join((str(x) for x in popt)))
    f.write('\n')
    f.write('\t'.join((str(x) for x in pcov)))
    f.close()
# function to read the fit for Process in the fit_dir for a specific RunNum:
def readFit(Process, RunNum):
    filename = fit_dir + 'fit_' + Process + '_run' + str(RunNum) + smearing_tag + '.dat'
    print('Reading fit from', filename)
    f = open(filename, 'r')
    for i,line in enumerate(f):
        if i == 0:
            if len(line.split())!= NCoeffs[Process]:
                print('Error: the number of coefficients found is insufficient: expected:', NCoeffs[Process], 'got:', len(line.split()))
                exit()
            else:
                popt = [float(x) for x in line.split()]
    pcov = [] # WARNING: COVARIANCE IS EMPTY HERE!
    return popt, pcov


def drive_mg_proc(runnum, mgloc, procloc, procname, CouplingsArray, nevents, nruns, ecm=14):
    filename = mgloc + '/' + procname + '_coupvar_run' + str(runnum) + '.dcmd'
    print('generating mg5input:', filename)
    ebeam1 = ecm*1000/2
    ebeam2 = ebeam1
    counter = 0
    for coups in CouplingsArray:
        if counter > nruns:
            break
        lhe = 'run_' + procname + '_' + str(RunNum) + '_' + '_'.join((coups)) + '/unweighted_events.lhe.gz'
        lhefile = mgloc + '/' + procloc + 'Events/' + lhe
        if os.path.exists(lhefile) is False:
            filestream = open(filename,'w')
            filestream.write('launch run_' + procname + '_' + str(RunNum) + '_' + '_'.join((coups)) + ' --accuracy=0.25 --points=300 --iterations=1\n0\n')
            filestream.write('set ebeam1 ' + str(ebeam1) + '\n')
            filestream.write('set ebeam2 ' + str(ebeam2) + '\n')
            filestream.write('set d3 ' + str(coups[0]) + '\n')
            filestream.write('set d4 ' + str(coups[1]) + '\n')
            #filestream.write('set cg1 ' + str(coups[2]) + '\n')
            #filestream.write('set cg2 ' + str(coups[3]) + '\n')
            #filestream.write('set ct1 ' + str(coups[4]) + '\n')
            #filestream.write('set cb1 ' + str(coups[5]) + '\n')
            filestream.write('set ct2 ' + str(coups[6]) + '\n')
            #filestream.write('set cb2 ' + str(coups[7]) + '\n')
            filestream.write('set ct3 ' + str(coups[8]) + '\n')
            if MODEL == 'HEFT6':
                filestream.write('set ct1 ' + str(coups[9]) + '\n')
            filestream.write('set nevents ' + str(nevents) + '\n')
            filestream.write('0')
            filestream.close()
            # run mg5 with the file generated
            runcommand = 'cat ' + filename
            p = subprocess.run(runcommand, shell=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, cwd=mgloc + '/' + procloc)
            runcommand = mgloc + '/' + procloc + '/bin/madevent ' + filename
            p = subprocess.Popen(runcommand, shell=True, text=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, cwd=mgloc + '/' + procloc)
            for line in iter(p.stdout.readline, b''):
                print(line)
            print(p.stdout)
            print(p.stderr)
                
            counter = counter + 1
    return counter

# function that runs herwig for specific final states
def run_herwig_proc(runnum, mgloc, hwloc, procloc, procname, CouplingsArray, nevents, nruns, ecm=100):
    print('Running Herwig from the input files previously generated, for:', procname, 'at Energy=', Energy)
    for coups in CouplingsArray:
        #print(lams)
        lhe = 'run_' + procname + '_' + str(RunNum) + '_' + '_'.join((coups)) + '/unweighted_events.lhe.gz'
        lhefile = mgloc + '/' + procloc + 'Events/' + lhe
        if os.path.exists(lhefile) is False:
            print('File', lhefile, 'does not exist, cannot run Herwig!')
            exit()
        # get the template and write the input file:
        # Signal is LO
        HerwigInputTemplate = getTemplate(HW_template[0])
        processname = 'HW-' + str(RunNum) + '_' + '_'.join((coups)) + '_' + FinalState
        hwinputfile = processname + '.in'
        parmtextsubs = {
            'PROCESSNAME' : processname, 
            'LHEFILE' : lhefile,
            'OUTPUTLOCATION' : 'events/',
            'FatAnalysis' : '#',
            'HwSimLibrary' : 'HwSim',
            'FinalState6b' : FinalState6b,
            'FinalStatebtau' : FinalStatebtau,
            'FinalStatebgamma' : FinalStatebgamma
            
        }
        print('\t\twriting', hwinputfile)
        writeFile(HerwigLocation + hwinputfile, HerwigInputTemplate.substitute(parmtextsubs) )

        # check if the root file already exists. if it does, only run if ReRun is set to True
        hwrunfile = processname + '.run'
        outputlocation = HerwigOutputLocation
        rootfile = outputlocation + processname + '.root'
        print("Checking rootfile:", rootfile)
        
        if os.path.exists(rootfile) is True:
            print('File', rootfile, 'exists')
        if os.path.exists(rootfile) is False or (os.path.exists(rootfile) is True and ReRunHerwig is True): # if the root file exists, do not proceed except if ReRun is true
                if os.path.exists(rootfile) is True and ReRunHerwig is True:
                    print('File', rootfile, 'exists, but have chosen to re-run!')
                # get the number of events in the corresponding lhe file:
                zgrepcommand = 'zgrep "= nevents" ' + lhefile
                print(zgrepcommand)
                p = subprocess.Popen(zgrepcommand, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, cwd=HerwigLocation)
                for line in iter(p.stdout.readline, b''):
                    nevents = float(line.split()[0])
                print('\t\tHerwig reading:', hwinputfile)
                readcommand = 'Herwig read ' + hwinputfile
                print(readcommand)
                p = subprocess.Popen(readcommand, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, cwd=HerwigLocation)
                for line in iter(p.stdout.readline, b''):
                    print('\t\t', line, end=' ')
                out, err = p.communicate()
                #print out, err
                print('\t\tHerwig running:', hwrunfile, 'for', nevents, 'events')
                runcommand = 'Herwig run ' + hwrunfile + ' -N' + str(int(nevents*Reduction_Fac[0]))
                p = subprocess.Popen(runcommand, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, cwd=HerwigLocation)
                for line in iter(p.stdout.readline, b''):
                    print('\t\t', line, end=' ')
                out, err = p.communicate()
                #print out, err


def run_herwig_proc_parallel(runnum, mgloc, hwloc, procloc, procname, CouplingsArray, nevents, nruns, ecm=100):
    print('Running Herwig from the input files previously generated, for:', procname, 'at Energy=', ecm)
    
    def worker(coups):
        lhe = 'run_' + procname + '_' + str(runnum) + '_' + '_'.join((coups)) + '/unweighted_events.lhe.gz'
        lhefile = mgloc + '/' + procloc + 'Events/' + lhe
        if not os.path.exists(lhefile):
            print('File', lhefile, 'does not exist, cannot run Herwig!')
            return  # Skip this job

        HerwigInputTemplate = getTemplate(HW_template[0])
        processname = 'HW-' + str(runnum) + '_' + '_'.join((coups)) + '_' + FinalState
        hwinputfile = processname + '.in'
        parmtextsubs = {
            'PROCESSNAME' : processname, 
            'LHEFILE' : lhefile,
            'OUTPUTLOCATION' : 'events/',
            'FatAnalysis' : '#',
            'HwSimLibrary' : 'HwSim',
            'FinalState6b' : FinalState6b,
            'FinalStatebtau' : FinalStatebtau,
            'FinalStatebgamma' : FinalStatebgamma
        }
        print('\t\twriting', hwinputfile)
        writeFile(HerwigLocation + hwinputfile, HerwigInputTemplate.substitute(parmtextsubs))

        hwrunfile = processname + '.run'
        outputlocation = HerwigOutputLocation
        rootfile = outputlocation + processname + '.root'
        print("Checking rootfile:", rootfile)

        rerun = (not os.path.exists(rootfile)) or (os.path.exists(rootfile) and ReRunHerwig)
        if os.path.exists(rootfile):
            print('File', rootfile, 'exists')
        if rerun:
            if os.path.exists(rootfile) and ReRunHerwig:
                print('File', rootfile, 'exists, but have chosen to re-run!')
            zgrepcommand = f'zgrep "= nevents" {lhefile}'
            print(zgrepcommand)
            p = subprocess.Popen(zgrepcommand, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, cwd=HerwigLocation)
            output, _ = p.communicate()
            try:
                nevents_local = float(output.decode().split()[0])
            except Exception:
                print('Could not parse nevents from LHE file, skipping', lhefile)
                return
            print('\t\tHerwig reading:', hwinputfile)
            readcommand = f'Herwig read {hwinputfile}'
            print(readcommand)
            p = subprocess.Popen(readcommand, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, cwd=HerwigLocation)
            for line in iter(p.stdout.readline, b''):
                print('\t\t', line.decode(), end=' ')
            p.communicate()

            print('\t\tHerwig running:', hwrunfile, 'for', nevents_local, 'events')
            runcommand = f'Herwig run {hwrunfile} -N{int(nevents_local * Reduction_Fac[0])}'
            print(runcommand)
            p = subprocess.Popen(runcommand, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, cwd=HerwigLocation)
            for line in iter(p.stdout.readline, b''):
                print('\t\t', line.decode(), end=' ')
            p.communicate()

    # Launch all Herwig runs in parallel
    Parallel(n_jobs=-1, backend="loky")(
        delayed(worker)(coups) for coups in CouplingsArray
    )

# function to read the analysis results and test the fit          
def test_fit_analysis(runnum, mgloc, procloc, procname, CouplingsArray, ntotal, popt):
    X = []
    Z = []
    EFFICIENCY = {}
    ZERR = []
    func_CX_proc = partial(func_CX, procname=Process)
    fracdiff_avg = 0.
    for coups in CouplingsArray:
        outputlocation = HerwigOutputLocation
        processname = 'HW-' + str(RunNum) + '_' + '_'.join((coups))
        rootfile = outputlocation + processname + '_' + FinalState + '.root'
        print('rootfile=', rootfile)
        analysisOutputfile = outputlocation + processname + '.smear' + smearing_tag + '.dat'
        if os.path.exists(analysisOutputfile)is False:
            print('File', analysisOutputfile, 'does not exist!')
            exit()
        else:
            print('File', analysisOutputfile, ' exists, reading results')
            zgrepcommand = 'cat ' + analysisOutputfile
            p = subprocess.Popen(zgrepcommand, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, cwd='.')
            for line in iter(p.stdout.readline, b''):
                efficiency = float(line.split()[0])
            #print('efficiency=', efficiency)
            coups_tuple = []
            for mm in range(len(coups)):
                coups_tuple.append(float(coups[mm]))
            X.append(tuple(coups_tuple))
            Z.append(float(efficiency))
            EFFICIENCY[tuple(coups_tuple)] = float(efficiency)
            # get the fitted XSEC
            eff_fit = func_CX_proc(coups_tuple, *popt)
            if efficiency != 0:
                fracdiff = abs(efficiency-eff_fit)/efficiency
            else:
                fracdiff = 0
            if fracdiff > 0.5:
                print(coups, efficiency)
                print('!!! xsec: real, fitted, frac diff =', efficiency, eff_fit, fracdiff)
            fracdiff_avg = fracdiff_avg + fracdiff
            #print(X)
    print('average fractional difference =', fracdiff_avg/ntotal)

def run_analysis(command):
    p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, cwd='.')
    for line in iter(p.stdout.readline, b''):
        print('\t\t', line.decode(), end=' ')
    out, err = p.communicate()
    print('\n')

# run the xgboost analysis - chatgpt modification using joblib
def run_analysis_xgboost(runnum, mgloc, hwloc, procloc, procname, CouplingsArray, nevents, nruns,
                        model_file, Backgrounds, Background_files, Backgrounds_xsec, xsS, initial_S, 
                        sig_factors, initial_B, idB, bkg_factors, Luminosity, Energy, training_seed, ecm=14):
    print('Running Analysis on the root files, for:', procname, 'at Energy=', Energy)
    X = []
    Z = []
    EFFICIENCY = {}
    EFFICIENCY_BKG = {}
    format = "%(asctime)s: %(message)s"
    logging.basicConfig(format=format, level=logging.INFO, datefmt="%H:%M:%S")

    jobs = []
    for coups in CouplingsArray:
        outputlocation = HerwigOutputLocation
        processname = 'HW-' + str(runnum) + '_' + '_'.join((coups))
        rootfile = outputlocation + processname + '_' + FinalState + '.root'
        analysisOutputfile = outputlocation + processname + smearing_tag + '.XGBOOST.dat'
        analysisInputfile = outputlocation + processname + '_var.smear' + smearing_tag + '.root'
        #print("Checking analysis output:", analysisOutputfile)

        if os.path.exists(analysisOutputfile) and not ReRunAnalysisXGBOOST:
            #print('File', analysisOutputfile, ' already exists, reading results')
            p = subprocess.Popen(f"cat {analysisOutputfile}", shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            for line in iter(p.stdout.readline, b''):
                if not line:
                    break
                efficiency = float(line.split()[0])
            #print('efficiency=', efficiency)
            coups_tuple = tuple(float(c) for c in coups)
            X.append(coups_tuple)
            Z.append(float(efficiency))
            EFFICIENCY[coups_tuple] = float(efficiency)
        else:
            if os.path.exists(analysisOutputfile) and ReRunAnalysisXGBOOST:
                print('File', analysisOutputfile, 'exists, but have chosen to re-run analysis!')
            if not os.path.exists(rootfile):
                print('Error: ROOT file:', rootfile, 'does not exist!')
                exit()
            print('running the XGBOOST analysis on the input file', analysisInputfile)
            jobs.append((
                model_file, analysisInputfile, Backgrounds, Background_files, Backgrounds_xsec, xsS,
                initial_S, sig_factors, initial_B, idB, bkg_factors, Luminosity, Energy, training_seed, smearing_tag
            ))

    # Run all XGBoost analyses in parallel
    if jobs:
        Parallel(n_jobs=-1, backend="loky")(
            delayed(apply_xgboost_write)(*args) for args in jobs
        )

    for bkg in Backgrounds:  # background loop
        processname = 'HW-' + str(bkg) + '_' + str(Energy)
        rootfile = BackgroundLocation + processname + '.root'
        analysisOutputfile = BackgroundLocation + processname + smearing_tag + '.XGBOOST.dat'
        #print("Checking analysis output:", analysisOutputfile) 
        if os.path.exists(analysisOutputfile):
            print('File', analysisOutputfile, ' exists, reading results')
            p = subprocess.Popen(f"cat {analysisOutputfile}", shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            for line in iter(p.stdout.readline, b''):
                if not line:
                    break
                efficiency = float(line.split()[0])
            print(bkg, 'efficiency=', efficiency)
            EFFICIENCY_BKG[bkg] = float(efficiency)
            continue
        else:
            print('Error, analysis for bkg', bkg, 'does not exist!', analysisOutputfile)
            exit()
    return np.transpose(X), Z, EFFICIENCY, EFFICIENCY_BKG

# run the analysis on signal and background USING XGBOOST             
def run_analysis_xgboost_threads(runnum, mgloc, hwloc, procloc, procname, CouplingsArray, nevents, nruns, trained_model, Backgrounds, Background_files, Backgrounds_xsec, xsS, initial_S, sig_factors, initial_B, idB, bkg_factors, Luminosity, Energy, training_seed, ecm=14):
    print('Running Analysis on the root files, for:', procname, 'at Energy=', Energy)
    X = []
    Z = []
    EFFICIENCY = {}
    EFFICIENCY_BKG = {}
    format = "%(asctime)s: %(message)s"
    logging.basicConfig(format=format, level=logging.INFO,datefmt="%H:%M:%S")
    #print(Max_Jobs)
    threads = list()
    for coups in CouplingsArray:
        #  write the analysis input file:
        outputlocation = HerwigOutputLocation
        processname = 'HW-' + str(RunNum) + '_' + '_'.join((coups))
        rootfile = outputlocation + processname + '_' + FinalState + '.root'
        analysisOutputfile = outputlocation + processname + smearing_tag + '.XGBOOST.dat'
        analysisInputfile = outputlocation + processname + '_var.smear' + smearing_tag + '.root'
        #print("Checking analysis output:", analysisOutputfile)
        if os.path.exists(analysisOutputfile) is True and ReRunAnalysis is False:
            print('File', analysisOutputfile, ' already exists, reading results')
            zgrepcommand = 'cat ' + analysisOutputfile
            p = subprocess.Popen(zgrepcommand, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, cwd='.')
            for line in iter(p.stdout.readline, b''):
                efficiency = float(line.split()[0])
            print('efficiency=', efficiency)
            coups_tuple = []
            for mm in range(len(coups)):
                coups_tuple.append(float(coups[mm]))
            X.append(tuple(coups_tuple))
            Z.append(float(efficiency))
            EFFICIENCY[tuple(coups_tuple)] = float(efficiency)
        elif (os.path.exists(analysisOutputfile) is False) or (os.path.exists(analysisOutputfile) is True and ReRunAnalysisXGBOOST is True): # if the root file exists, do not proceed except if ReRun is true
                if os.path.exists(analysisOutputfile) is True and ReRunAnalysisXGBOOST is True:
                    print('File', analysisOutputfile, 'exists, but have chosen to re-run analysis!')
                if os.path.exists(rootfile) is False:
                    print('Error: ROOT file:', rootfile, 'does not exist!')
                    exit()
                print('running the XGBOOST analysis on the input file', analysisInputfile)
                print('Launching: apply_xgboost_write with:', trained_model, analysisInputfile, Backgrounds, Background_files, Backgrounds_xsec, xsS, initial_S, sig_factors, initial_B, idB, bkg_factors, Luminosity, Energy, training_seed, smearing_tag)
                
                x = threading.Thread(target=apply_xgboost_write, args=(trained_model, analysisInputfile, Backgrounds, Background_files, Backgrounds_xsec, xsS, initial_S, sig_factors, initial_B, idB, bkg_factors, Luminosity, Energy, training_seed, smearing_tag))
                #x = multiprocessing.Process(target=apply_xgboost_write, args=(trained_model, analysisInputfile, Backgrounds, Background_files, Backgrounds_xsec, xsS, initial_S, sig_factors, initial_B, idB, bkg_factors, Luminosity, Energy, training_seed,))
                x.start()
                x.join()
                print(x.exitcode) 
                #threads.append(x)
    #for index, thread in enumerate(threads):
    #    logging.info("Main    : before joining thread %d.", index)
    #    thread.join()
    #    logging.info("Main    : thread %d done", index)
    for bkg in Backgrounds: # background loop
        processname = 'HW-' + str(bkg) + '_' + str(Energy)
        rootfile = BackgroundLocation + processname + '.root'
        analysisOutputfile = BackgroundLocation + processname + smearing_tag + '.XGBOOST.dat'
        print("Checking analysis output:", analysisOutputfile) 
        if os.path.exists(analysisOutputfile) is True:
            print('File', analysisOutputfile, ' exists, reading results')
            zgrepcommand = 'cat ' + analysisOutputfile
            p = subprocess.Popen(zgrepcommand, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, cwd='.')
            for line in iter(p.stdout.readline, b''):
                efficiency = float(line.split()[0])
            print(bkg, 'efficiency=', efficiency)
            EFFICIENCY_BKG[bkg] = float(efficiency)
            continue
        if os.path.exists(analysisOutputfile) is False: # if the root file exists, do not proceed except if ReRun is true
                print('Error, analysis for bkg', bkg, 'does not exist!')
                exit()
    return np.transpose(X), Z, EFFICIENCY, EFFICIENCY_BKG


    
# run the analysis on signal and background               
def run_analysis_proc(runnum, mgloc, hwloc, procloc, procname, CouplingsArray, nevents, nruns, ecm=14):
    print('Running Analysis on the root files, for:', procname, 'at Energy=', Energy)
    X = []
    Z = []
    EFFICIENCY = {}
    EFFICIENCY_BKG = {}
    format = "%(asctime)s: %(message)s"
    logging.basicConfig(format=format, level=logging.INFO,datefmt="%H:%M:%S")
    #print(Max_Jobs)
    threads = list()
    for coups in CouplingsArray:
        #  write the analysis input file:
        outputlocation = HerwigOutputLocation
        processname = 'HW-' + str(RunNum) + '_' + '_'.join((coups))
        rootfile = outputlocation + processname + '_' + FinalState + '.root'
        analysisOutputfile = outputlocation + processname + '.smear' + smearing_tag + '.dat'
        analysisInputfile = outputlocation + processname + '.input'
        analysisInputstream = open(analysisInputfile,'w') 
        print("Checking analysis output:", analysisOutputfile)
        if os.path.exists(analysisOutputfile) is True and ReRunAnalysis is False:
            print('File', analysisOutputfile, ' already exists, reading results')
            zgrepcommand = 'cat ' + analysisOutputfile
            p = subprocess.Popen(zgrepcommand, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, cwd='.')
            for line in iter(p.stdout.readline, b''):
                efficiency = float(line.split()[0])
            print('efficiency=', efficiency)
            coups_tuple = []
            for mm in range(len(coups)):
                coups_tuple.append(float(coups[mm]))
            X.append(tuple(coups_tuple))
            Z.append(float(efficiency))
            EFFICIENCY[tuple(coups_tuple)] = float(efficiency)
        if os.path.exists(analysisOutputfile) is False or (os.path.exists(analysisOutputfile) is True and ReRunAnalysis is True): # if the root file exists, do not proceed except if ReRun is true
                if os.path.exists(analysisOutputfile) is True and ReRunAnalysis is True:
                    print('File', analysisOutputfile, 'exists, but have chosen to re-run analysis!')
                if os.path.exists(rootfile) is False:
                    print('Error: ROOT file:', rootfile, 'does not exist!')
                    exit()
                elif os.path.exists(rootfile) is True:
                    analysisInputstream.write(rootfile + '\n')
                    analysisInputstream.close()
                print('running the analysis', ExecutableSmear[Energy], 'on the input file', analysisInputfile)
                analysiscommand = ExecutableSmear[Energy] + ' ' + analysisInputfile
                print('Launching:', analysiscommand)

                #p = subprocess.Popen(analysiscommand, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, cwd='.')
                #for line in iter(p.stdout.readline, b''):
                #print('\t\t', line, end=' ')
                #out, err = p.communicate()
                #print('\n')
                x = threading.Thread(target=run_analysis, args=(analysiscommand,))
                threads.append(x)
                x.start()
    for index, thread in enumerate(threads):
        #logging.info("Main    : before joining thread %d.", index)
        thread.join()
        logging.info("Main    : thread %d done", index)
    for bkg in Backgrounds: # background loop
        processname = 'HW-' + str(bkg) + '_' + str(Energy)
        rootfile = BackgroundLocation + processname + '.root'
        analysisOutputfile = BackgroundLocation + processname + '.smear' + smearing_tag + '.dat'
        analysisInputfile = BackgroundLocation + processname + '.input'
        analysisInputstream = open(analysisInputfile,'w') 
        print("Checking analysis output:", analysisOutputfile) 
        if os.path.exists(analysisOutputfile) is True and ReRunAnalysis is False:
            print('File', analysisOutputfile, ' already exists, reading results')
            zgrepcommand = 'cat ' + analysisOutputfile
            p = subprocess.Popen(zgrepcommand, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, cwd='.')
            for line in iter(p.stdout.readline, b''):
                efficiency = float(line.split()[0])
            print(bkg, 'efficiency=', efficiency)
            EFFICIENCY_BKG[bkg] = float(efficiency)
            continue
        if os.path.exists(analysisOutputfile) is False or (os.path.exists(analysisOutputfile) is True and ReRunAnalysis is True): # if the root file exists, do not proceed except if ReRun is true
                if os.path.exists(analysisOutputfile) is True and ReRunAnalysis is True:
                    print('File', analysisOutputfile, 'exists, but have chosen to re-run analysis!')
                if os.path.exists(rootfile) is False:
                    print('Error: ROOT file:', rootfile, 'does not exist!')
                    exit()
                elif os.path.exists(rootfile) is True:
                    analysisInputstream.write(rootfile + '\n')
                    analysisInputstream.close()
                print('running the analysis', ExecutableSmear[Energy], 'on the input file', analysisInputfile)
                analysiscommand = ExecutableSmear[Energy] + ' ' + analysisInputfile
                p = subprocess.Popen(analysiscommand, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, cwd='.')
                for line in iter(p.stdout.readline, b''):
                        print('\t\t', line, end=' ')
                out, err = p.communicate()
                print('\n')
    return np.transpose(X), Z, EFFICIENCY, EFFICIENCY_BKG


def contour_pvalue_ct3d4_marginalized(procname, plotname, plottitle, fit_coeffs_xsec, fit_coeffs_eff, sigma_bkg, var1, var2, xlim, ylim, axext='', figext='', smtext=True, starsize=15, setxlabel=True, setylabel=True, nbins=200, savefig=True,variables=variables, variables_latex=variables_latex, labelsize=20, normalbar=True, contours=np.arange(0, 10, 0.5),norm_to_zeroth=True, lumi=Luminosity):
    output = procname + '_' + plotname + '_' + var1 + '_' + var2 + '_including_constraints_marginalized'
    print('Plotting', output)
    nvar1 = [key for key, value in variables.items() if value == var1][0]
    nvar2 = [key for key, value in variables.items() if value == var2][0]
    #print(var1, var2)
    #print(nvar1, nvar2)
    # construct the axes for the plot
    # no need to modify this if you just need one plot
    gs = gridspec.GridSpec(4, 4)
    if figext == '':
        fig = plt.figure()
    else:
        fig = figext
    if axext == '':
        ax = fig.add_subplot(111)
    else:
        ax=axext
    ax.grid(False)
    ax.set_title(plottitle)
    
    # set the ticks, labels and limits etc.
    xlab = '$' + variables_latex[nvar1] + '$'
    ylab = '$' + variables_latex[nvar2] + '$'
    if setylabel == True:
        ax.set_ylabel(ylab, fontsize=labelsize)
    if setxlabel == True:
        ax.set_xlabel(xlab, fontsize=labelsize)
        
    # set the limits on the x and y axes if required below:
    ymin = ylim[0]
    ymax = ylim[1]
    xmin = xlim[0]
    xmax = xlim[1]
    plt.xlim([xmin,xmax])
    plt.ylim([ymin,ymax])
    ctexts = []
    cvartexts = []
    for i in range(0, len(variables.keys())):
        #if i != nvar1 and i != nvar2:
        #   ctext = variables[i] + '=0'
        #    ctexts.append(ctext)
        #else:
        cvartexts.append(variables[i])
    #print(ctexts)
    #fstr = 'partial(func_t_CX, ' + ','.join([ct for ct in ctexts]) + ', procname=Process)'
    fstr = 'partial(func_t_CX, procname=Process)'
    global func_CX_partial
    func_CX_partial = eval(fstr)
    #print('func_CX_partial=', func_CX_partial)
    #print(fit_coeffs)
    global fit_coeffs_g_xsec
    fit_coeffs_g_xsec = fit_coeffs_xsec
    global fit_coeffs_g_eff
    fit_coeffs_g_eff = fit_coeffs_eff
    #print(cvartexts[0], cvartexts[1])

    # functions for xsec and significance:
    print(cvartexts)
    # 0     1      2      3     
    #['c3', 'ct2', 'ct3', 'd4']
    feval_xsec = 'func_CX_partial(' + cvartexts[0] +'=0,' + cvartexts[1] + '=0,' + cvartexts[2] + '=x1,' + cvartexts[3] + '=x2,coeffs=fit_coeffs_g_xsec)'
    feval_eff = 'func_CX_partial(' + cvartexts[0] +'=0,' + cvartexts[1] + '=0,'  + cvartexts[2] + '=x1,' + cvartexts[3] + '=x2,coeffs=fit_coeffs_g_eff)'

    feval_xsec_g = 'func_CX_partial(' + cvartexts[0] +'=x3,' + cvartexts[1] + '=x4,' + cvartexts[2] + '=x1,' + cvartexts[3] + '=x2,coeffs=fit_coeffs_g_xsec)'
    feval_eff_g  = 'func_CX_partial(' + cvartexts[0] +'=x3,' + cvartexts[1] + '=x4,' + cvartexts[2] + '=x1,' + cvartexts[3] + '=x2,coeffs=fit_coeffs_g_eff)'
  
    
    func_fin = lambda x1, x2: eval(feval_xsec) * sig_factors * eval(feval_eff) * lumi * 1000. / math.sqrt(sigma_bkg * bkg_factors * lumi + Systematics**2 * (sigma_bkg*bkg_factors)**2 * lumi**2)
    #func_fin = lambda x1, x2: significance(eval(feval_xsec) * eval(feval_eff) * lumi * 1000., sigma_bkg * lumi, Systematics)

    pfunc_fin = lambda x1, x2: (1 - 2 * scipy.special.ndtr(-(eval(feval_xsec) * sig_factors * eval(feval_eff) * lumi * 1000. / math.sqrt(sigma_bkg * bkg_factors * lumi + Systematics**2 * (sigma_bkg*bkg_factors)**2 * lumi**2))))

    pfunc_fin_gaussrw = lambda x1, x2, x3, x4: (1 - 2 * scipy.special.ndtr(-(eval(feval_xsec_g) * sig_factors * eval(feval_eff_g) * lumi * 1000. / math.sqrt(sigma_bkg * bkg_factors * lumi + Systematics**2 * (sigma_bkg*bkg_factors)**2 * lumi**2)))) * gaussian(x3, 0, constraints[Energy][0]) * gaussian(x4, 0, constraints[Energy][1])

    # SM significance: 
    feval_xsec_sm = 'func_CX_partial(' + cvartexts[0] +'=0,' + cvartexts[1] + '=0,' + cvartexts[2] + '=0,' + cvartexts[3] + '=0,coeffs=fit_coeffs_g_xsec)'
    feval_eff_sm = 'func_CX_partial(' + cvartexts[0] +'=0,' + cvartexts[1] + '=0,' + cvartexts[2] + '=0,' + cvartexts[3] + '=0,coeffs=fit_coeffs_g_eff)'
    
    sm_signif= eval(feval_xsec_sm) * sig_factors * eval(feval_eff_sm) * 1000. * lumi / math.sqrt(sigma_bkg * bkg_factors * lumi + Systematics**2 * (sigma_bkg*bkg_factors)**2 * lumi**2)

    # if SM is the "null" hypothesis:
    # SM number of events:
    S_SM = eval(feval_xsec_sm) * sig_factors * eval(feval_eff_sm) * 1000. * lumi
    # SM total uncertainty, including the background uncertainty:
    delta_SM =  math.sqrt(S_SM + sigma_bkg * bkg_factors * lumi + Systematics**2 * (sigma_bkg*bkg_factors)**2 * lumi**2) 
    # {c_i} number of events in the 4D model:
    S_i_4D = lambda x1, x2, x3, x4: eval(feval_xsec_g) * sig_factors * eval(feval_eff_g) * lumi * 1000.
    # {c_i} number of events in 2D model:
    S_i_2D = lambda x1, x2: eval(feval_xsec) * sig_factors * eval(feval_eff) * lumi * 1000.
    # significance versus the SM in the 4D model:
    func_fin_SM_4D = lambda x1, x2, x3, x4: np.power( (S_SM - S_i_4D(x1, x2, x3, x4))/delta_SM, 2)
    # significance versus the SM in the 2D mode: 
    func_fin_SM_2D = lambda x1, x2: np.power( (S_SM - S_i_2D(x1, x2))/delta_SM,2)
    
    # p-value in the 4D model (NO gaussian RW):
    pfunc_fin_SM_4D = lambda x1, x2, x3, x4: 1/(np.sqrt(2.*np.pi)*delta_SM)*np.exp(-func_fin_SM_4D(x1, x2, x3, x4)/2)
    # p-value in the 4D model (WITH gaussian RW):
    pfunc_fin_SM_4D_g = lambda x1, x2, x3, x4: 1/(np.sqrt(2.*np.pi)*delta_SM)*np.exp(-func_fin_SM_4D(x1, x2, x3, x4)/2) * gaussian(x3, 0, constraints[Energy][0]) * gaussian(x4, 0, constraints[Energy][1])
    # p-value in the 2D model:
    pfunc_fin_SM_2D = lambda x1, x2:  1/(np.sqrt(2.*np.pi)*delta_SM)*np.exp(-func_fin_SM_2D(x1, x2)/2)

    print('pfunc_fin_SM_4D_g(0,0,0,0)=',pfunc_fin_SM_4D_g(0,0,0,0))
    print('pfunc_fin_SM_2D(0,0)=',pfunc_fin_SM_2D(0,0))
    print("sigma_sig before anal. [fb]=", eval(feval_xsec_sm)*1000*sig_factors)
    print("analysis eff. on signal=", eval(feval_eff_sm))
    print("sigma_bkg after anal. [fb]=", sigma_bkg * bkg_factors)
    print("sigma sig SM after anal. [fb]=",eval(feval_xsec_sm) * sig_factors * eval(feval_eff_sm) * 1000)
    print("N(bkg)@lumi=", sigma_bkg * bkg_factors * lumi)
    print("N(sig SM)@lumi=", eval(feval_xsec_sm) * sig_factors * eval(feval_eff_sm) * lumi * 1000.) 
    print("SM significance=", sm_signif)
    
    # The two-dimensional p-value (all other coefficients zero):
    x = np.linspace(xlim[0], xlim[1], nbins)
    y = np.linspace(ylim[0], ylim[1], nbins)
    X, Y = np.meshgrid(x,y)
    P = pfunc_fin_SM_2D(X,Y) #func_fin(X,Y)
    #P = P/pfunc_fin_SM_2D(0,0)
    #print(np.amax(P))
    # convert to chi-sq.:
    chisq = stats.chi2.isf(P,2)
    chisq_sub = chisq - np.amin(chisq)
    print('np.amin(chisq)=',np.amin(chisq))
    
    # The four-dimensional p-value: (ct3, d4, c3, ct2)
    x1 = np.linspace(xlim[0], xlim[1], nbins) # ct3
    x2 = np.linspace(ylim[0], ylim[1], nbins) # d4
    nsigma = 10 # number of standard deviations away from the central value
    x3 = np.linspace(-nsigma*constraints[Energy][0],nsigma*constraints[Energy][0], nbins) # c3 limits
    x4 = np.linspace(-nsigma*constraints[Energy][1],nsigma*constraints[Energy][1], nbins) # ct2 limits
    #x3 = np.zeros(nbins)
    #x4 = np.zeros(nbins)
    X1, X2, X3, X4 = np.meshgrid(x1,x2,x3,x4)
    P_g = pfunc_fin_SM_4D_g(X1,X2,X3,X4)
    P_g_marg = np.apply_over_axes(np.sum, P_g, [2,3])
    P_g_marg_s = P_g_marg.reshape(P_g_marg.shape[0], P_g_marg.shape[1])
    P_g_marg_bar = P_g_marg_s*2*nsigma*constraints[Energy][0]*2*nsigma*constraints[Energy][1]/nbins/nbins
    # convert to chi-sq.:
    chisq_marg = stats.chi2.isf(P_g_marg_bar,2)
    chisq_marg_sub = chisq_marg - np.amin(chisq_marg)
    print('np.amin(chisq_marg)=',np.amin(chisq_marg))

    
    #cont = ax.contourf(X, Y, P, contours, cmap='Spectral', extend='max')

    # do the one-dimensional marginalizations:
    P_g_marg_d4 = np.apply_over_axes(np.sum, P_g, [0, 2, 3])
    P_g_marg_ct3 = np.apply_over_axes(np.sum, P_g, [1, 2, 3])
    P_g_marg_d4_s = P_g_marg_d4.reshape(P_g_marg_d4.shape[1])
    P_g_marg_ct3_s = P_g_marg_ct3.reshape(P_g_marg_ct3.shape[0])
    P_g_marg_d4_bar = P_g_marg_d4_s*2*nsigma*constraints[Energy][0]*2*nsigma*constraints[Energy][1]*(xlim[1]-xlim[0])/nbins/nbins/nbins
    P_g_marg_ct3_bar = P_g_marg_d4_s*2*nsigma*constraints[Energy][0]*2*nsigma*constraints[Energy][1]*(ylim[1]-ylim[0])/nbins/nbins/nbins
    chisq_marg_d4 = stats.chi2.isf(P_g_marg_d4_bar,1)
    chisq_marg_d4_sub = chisq_marg_d4 - np.amin(chisq_marg_d4)
    chisq_marg_ct3 = stats.chi2.isf(P_g_marg_ct3_bar,1)
    chisq_marg_ct3_sub = chisq_marg_ct3 - np.amin(chisq_marg_ct3)

    # remove inf and nans if necessary:
    #chisq_marg_d4_sub[np.isinf(chisq_marg_d4_sub)] = np.nan
    #chisq_marg_d4_sub[np.isnan(chisq_marg_d4_sub)] = np.nanmax(chisq_marg_d4_sub, axis=0)
    #chisq_marg_ct3_sub[np.isinf(chisq_marg_ct3_sub)] = np.nan
    #chisq_marg_ct3_sub[np.isnan(chisq_marg_ct3_sub)] = np.nanmax(chisq_marg_ct3_sub, axis=0)
    #print(chisq_marg_d4_sub)
    #print(chisq_marg_ct3_sub)

    # interpolate the 1D functions: 
    func_chisq_1D_d4 = interp1d(x2,chisq_marg_d4_sub, fill_value="extrapolate")
    func_chisq_1D_ct3 = interp1d(x1,chisq_marg_ct3_sub, fill_value="extrapolate")
    #func_chisq_1D_d4 =  make_interp_spline(x2,chisq_marg_d4_sub, k=3)
    #func_chisq_1D_ct3 =  make_interp_spline(x1,chisq_marg_ct3_sub, k=3)
    
    # construction functions to find 1 and 2 sigma limits on d4 and ct3 (from chi-sq min). 
    def func_d4_1sigma(x): return (func_chisq_1D_d4(x) - 0.99)
    def func_d4_2sigma(x): return (func_chisq_1D_d4(x) - 3.84)
    def func_ct3_1sigma(x): return (func_chisq_1D_ct3(x) - 0.99)
    def func_ct3_2sigma(x): return (func_chisq_1D_ct3(x) - 3.84)
        
    # guesses for the locations of the solutions in 1D [change with energy]:
    d4_min_1 = {}
    d4_max_1 = {}
    d4_min_2 = {}
    d4_max_2 = {}
    d4_min_1[13.6] = -10
    d4_max_1[13.6] = 10
    d4_min_2[13.6] = -35 # triple-ins
    d4_max_2[13.6] = 80 # triple-ins
    #d4_min_2[13.6] = -50 # double-ins
    #d4_max_2[13.6] = 20 # double-ins
    
    d4_min_1[100] = -5
    d4_max_1[100] = 32
    d4_min_2[100] = -5
    d4_max_2[100] = 32

    ct3_min_1 = {}
    ct3_max_1 = {}
    ct3_min_2 = {}
    ct3_max_2 = {}
    ct3_min_1[13.6] = -1
    ct3_max_1[13.6] = 2
    ct3_min_2[13.6] = -2
    ct3_max_2[13.6] = 4
    
    ct3_min_1[100] = -0.1
    ct3_max_1[100] = 0.6
    ct3_min_2[100] = -0.8
    ct3_max_2[100] = 0.5

    # calculate and print out the solutions:
    #print('d4@68% CL:', fsolve(func_d4_1sigma, d4_min_1[Energy]), fsolve(func_d4_1sigma, d4_max_1[Energy]))
    #print('d4@95% CL:', fsolve(func_d4_2sigma, d4_max_2[Energy]), fsolve(func_d4_2sigma, d4_max_2[Energy]))
    #print('ct3@68% CL:', fsolve(func_ct3_1sigma, ct3_min_1[Energy]), fsolve(func_ct3_1sigma, ct3_max_1[Energy]))
    #print('ct3@95% CL:', fsolve(func_ct3_2sigma, ct3_max_2[Energy]), fsolve(func_ct3_2sigma, ct3_max_2[Energy]))
    print('d4@68% CL:', fsolve(func_d4_1sigma, [d4_min_1[Energy], d4_max_1[Energy]]))
    print('d4@95% CL:', fsolve(func_d4_2sigma, [d4_min_2[Energy], d4_max_2[Energy]]))
    print('ct3@68% CL:', fsolve(func_ct3_1sigma, [ct3_min_1[Energy], ct3_max_1[Energy]]))
    print('ct3@95% CL:', fsolve(func_ct3_2sigma, [ct3_min_2[Energy], ct3_max_2[Energy]]))

    # plot the contours:
    #ax.clabel(cont)#, inline=True)
    ax.plot(0,0,marker='*',ms=starsize, color='black')
    #cont2 = ax.contour(X, Y, P_g_marg_bar, contours, extend='max', colors=('black'), label='4D')
    #cont = ax.contour(X, Y, P, contours, extend='max', colors=('red'), linestyles=('--'), label='2D')
    cont = ax.contour(X, Y, chisq_marg_sub, contours, extend='max', colors=('black', 'red'), linestyles=('-','--'))
    labels = ['$1\\sigma$', '$2\\sigma$']
    for i in range(len(labels)):
        cont.collections[i].set_label(labels[i])

    #cont = ax.contour(X, Y, chisq_sub, contours, extend='max', colors=('red'), linestyles=('--'), label='2D')

    # add constraints:
    if constraints[Energy][nvar1] != -1:
        ax.axvline(x=constraints[Energy][nvar1],  linewidth=0.5, color = 'k', ls='--')
        ax.axvline(x=-constraints[Energy][nvar1],  linewidth=0.5, color = 'k', ls='--')
    if constraints[Energy][nvar2] != -1:
        ax.axhline(y=constraints[Energy][nvar2],  linewidth=0.5, color = 'k', ls='--')
        ax.axhline(y=-constraints[Energy][nvar2],  linewidth=0.5, color = 'k', ls='--')
    
    if smtext == True:
        ax.text(0.53, 0.53,"SM", transform=ax.transAxes)
    if normalbar == True:
        plt.colorbar(cont)
    #handles, labels = cs.legend_elements()

    # after you’ve done your contour call…
    black_line = mlines.Line2D([], [], color='black', linestyle='-',
                            label='$1\\sigma$')
    red_line   = mlines.Line2D([], [], color='red',   linestyle='--',
                           label='$2\\sigma$')

    ax.legend(handles=[black_line, red_line],
          loc="upper right", frameon=False, prop={'size':8})
        
    #ax.legend()
    #ax.legend(loc="upper right", numpoints=1, frameon=False, prop={'size':8}, handles=[cont, cont2])
    ax.yaxis.set_minor_locator(MultipleLocator(5))
    if Energy == 100:
        ax.xaxis.set_minor_locator(MultipleLocator(0.05))
    elif Energy == 13.6:
        ax.xaxis.set_minor_locator(MultipleLocator(0.2))
    if savefig == True:
        # save the figure
        print('saving the figure')
        # save the figure in PDF format
        infile = output + '.dat'
        print('---')
        print('output in', infile.replace('.dat','.pdf'))
        plt.savefig(plot_dir + infile.replace('.dat','.pdf'), bbox_inches='tight')
        plt.close(fig)
        
    return cont




def contour_pvalue_only_old(procname, plotname, plottitle, fit_coeffs_xsec, fit_coeffs_eff, sigma_bkg, var1, var2, plotlimits, searchlimits, deltac3=-1, axext='', figext='', smtext=True, starsize=15, setxlabel=True, setylabel=True, nbins=400, savefig=True,variables=variables, variables_latex=variables_latex, labelsize=20, normalbar=True, contours=np.arange(0, 10, 0.5),norm_to_zeroth=True, lumi=Luminosity):
    output = procname + '_' + plotname + '_' + var1 + '_' + var2
    print('Plotting', output)
    nvar1 = [key for key, value in variables.items() if value == var1][0]
    nvar2 = [key for key, value in variables.items() if value == var2][0]
    nvar3 = [key for key, value in variables.items() if value != var1 and value != var2][0]
    nvar4 = [key for key, value in variables.items() if value != var1 and value != var2][1]

    #print(var1, var2)
    print('nvar1, nvar2=', nvar1, nvar2)
    #print(nvar3, nvar4)
   
        
    # set the limits on the x and y axes if required below:
    ymin = plotlimits[Energy][nvar2][0]
    ymax = plotlimits[Energy][nvar2][1]
    xmin = plotlimits[Energy][nvar1][0]
    xmax = plotlimits[Energy][nvar1][1]
  
    ctexts = []
    cvartexts = []
    for i in range(0, len(variables.keys())):
        cvartexts.append(variables[i])
    fstr = 'partial(func_t_CX, procname=Process)'
    global func_CX_partial
    func_CX_partial = eval(fstr)
    #print('func_CX_partial=', func_CX_partial)
    #print(fit_coeffs)
    global fit_coeffs_g_xsec
    fit_coeffs_g_xsec = fit_coeffs_xsec
    global fit_coeffs_g_eff
    fit_coeffs_g_eff = fit_coeffs_eff
    #print(cvartexts[0], cvartexts[1])

    # functions for xsec and significance:
    print(cvartexts)
    # 0     1      2      3     
    #['c3', 'ct2', 'ct3', 'd4']
    print(cvartexts[nvar3], cvartexts[nvar4], cvartexts[nvar1], cvartexts[nvar2])
    feval_xsec = 'func_CX_partial(' + cvartexts[nvar3] +'=0,' + cvartexts[nvar4] + '=0,' + cvartexts[nvar1] + '=x1,' + cvartexts[nvar2] + '=x2,coeffs=fit_coeffs_g_xsec)'
    feval_eff = 'func_CX_partial(' + cvartexts[nvar3] +'=0,' + cvartexts[nvar4] + '=0,'  + cvartexts[nvar1] + '=x1,' + cvartexts[nvar2] + '=x2,coeffs=fit_coeffs_g_eff)'
      
    func_fin = lambda x1, x2: eval(feval_xsec) * sig_factors * eval(feval_eff) * lumi * 1000. / math.sqrt(sigma_bkg * lumi + Systematics**2 * (sigma_bkg)**2 * lumi**2)

    pfunc_fin = lambda x1, x2: (1 - 2 * scipy.special.ndtr(-(eval(feval_xsec) * sig_factors * eval(feval_eff) * lumi * 1000. / math.sqrt(sigma_bkg * lumi + Systematics**2 * (sigma_bkg)**2 * lumi**2))))

    # SM significance: 
    feval_xsec_sm = 'func_CX_partial(' + cvartexts[0] +'=0,' + cvartexts[1] + '=0,' + cvartexts[2] + '=0,' + cvartexts[3] + '=0,coeffs=fit_coeffs_g_xsec)'
    feval_eff_sm = 'func_CX_partial(' + cvartexts[0] +'=0,' + cvartexts[1] + '=0,' + cvartexts[2] + '=0,' + cvartexts[3] + '=0,coeffs=fit_coeffs_g_eff)'
    
    sm_signif= eval(feval_xsec_sm) * sig_factors * eval(feval_eff_sm) * 1000. * lumi / math.sqrt(sigma_bkg * lumi + Systematics**2 * (sigma_bkg)**2 * lumi**2)

    # if SM is the "null" hypothesis:
    # SM number of events:
    S_SM = eval(feval_xsec_sm) * sig_factors * eval(feval_eff_sm) * 1000. * lumi
    # SM total uncertainty, including the background uncertainty:
    delta_SM =  math.sqrt(S_SM + sigma_bkg * lumi + Systematics**2 * (sigma_bkg)**2 * lumi**2) 
    # {c_i} number of events in 2D model:
    S_i_2D = lambda x1, x2: eval(feval_xsec) * sig_factors * eval(feval_eff) * lumi * 1000.
    # significance versus the SM in the 2D mode: 
    func_fin_SM_2D = lambda x1, x2: np.power( (S_SM - S_i_2D(x1, x2))/delta_SM,2)

    # p-value in the 2D model:
    #pfunc_fin_SM_2D = lambda x1, x2:  1/(np.sqrt(2.*np.pi)*delta_SM)*np.exp(-func_fin_SM_2D(x1, x2)/2)
    # print('pfunc_fin_SM_2D(0,0)=',pfunc_fin_SM_2D(0,0))
    
    print("sigma_sig before anal. [fb]=", eval(feval_xsec_sm)*1000*sig_factors)
    print("analysis eff. on signal=", eval(feval_eff_sm))
    print("sigma_bkg after anal. [fb]=", sigma_bkg)
    print("sigma sig SM after anal. [fb]=",eval(feval_xsec_sm) * sig_factors * eval(feval_eff_sm) * 1000)
    print("N(bkg)@lumi=", sigma_bkg * lumi)
    print("N(sig SM)@lumi=", eval(feval_xsec_sm) * sig_factors * eval(feval_eff_sm) * lumi * 1000.) 
    print("SM significance=", sm_signif)
    

    
    x = np.linspace(xmin, xmax, nbins)
    y = np.linspace(ymin, ymax, nbins)
    dx = x[1] - x[0]
    dy = y[1] - y[0]
    X, Y = np.meshgrid(x,y)
    #P = pfunc_fin_SM_2D(X,Y)  #func_fin(X,Y)
    chisq_prior = 0
    if deltac3 > 0:
        #P=P*gaussian(X, 0, deltac3) # REWEIGH BY C3 PRIOR if deltac3 > 0
        chisq_prior = (X - 0)**2 / deltac3**2
    chisq = func_fin_SM_2D(X,Y) + chisq_prior

    #P = P/pfunc_fin_SM_2D(0,0)
    #print(np.amax(P))
    # convert to chi-sq.:
    #chisq = stats.chi2.isf(P,2)
        
    chisq_sub = chisq - np.amin(chisq)
    print('np.amin(chisq)=',np.amin(chisq))

    # MARGINALIZATION ATTEMPTS BELOW
    # do the one-dimensional marginalizations:
    # sum over the marginalized direction
    #P_marg_nvar1 = np.apply_over_axes(np.sum, P, [1]) * dy
    #P_marg_nvar2 = np.apply_over_axes(np.sum, P, [0]) * dx
    # change the shape:
    #P_marg_nvar1_s = P_marg_nvar1.reshape(P_marg_nvar1.shape[0])
    #P_marg_nvar2_s = P_marg_nvar2.reshape(P_marg_nvar2.shape[1])
    #print('P_marg_nvar1_s=',P_marg_nvar1_s)
    #print('P_marg_nvar2_s=',P_marg_nvar2_s)
    # convert each probability to chisq:
    #chisq_marg_nvar1 = stats.chi2.isf(P_marg_nvar1_s,1)
    #chisq_marg_nvar2 = stats.chi2.isf(P_marg_nvar2_s,1)
    #print('chisq_marg_nvar1=', chisq_marg_nvar1)
    #print('chisq_marg_nvar2=', chisq_marg_nvar2)
    # remove infinities and nans:
    #chisq_marg_nvar1 = np.nan_to_num(chisq_marg_nvar1)
    #chisq_marg_nvar2 = np.nan_to_num(chisq_marg_nvar2)
    # subtract the minimum of chisq
    #chisq_marg_nvar1_sub = chisq_marg_nvar1 - np.amin(chisq_marg_nvar1)
    #chisq_marg_nvar2_sub = chisq_marg_nvar2 - np.amin(chisq_marg_nvar2)
    #print('chisq_marg_nvar1_sub=', chisq_marg_nvar1_sub)
    #print('chisq_marg_nvar2_sub=', chisq_marg_nvar2_sub)
    # remove infinities and nans:
    #chisq_marg_nvar1_sub = np.nan_to_num(chisq_marg_nvar1_sub)
    #chisq_marg_nvar2_sub = np.nan_to_num(chisq_marg_nvar2_sub)
    #print('chisq_marg_nvar1_sub=', chisq_marg_nvar1_sub)
    #print('chisq_marg_nvar2_sub=', chisq_marg_nvar2_sub)

    # profiling attempt:
    chisq_marg_nvar1_sub = np.min(chisq_sub, axis=1) 
    chisq_marg_nvar2_sub = np.min(chisq_sub, axis=0)

    x1 = np.linspace(xmin, xmax, nbins) # nvar1
    x2 = np.linspace(ymin, ymax, nbins) # nvar2
    
    # interpolate the 1D functions: 
    func_chisq_1D_nvar1 = interp1d(x1,chisq_marg_nvar1_sub, fill_value="extrapolate")
    func_chisq_1D_nvar2 = interp1d(x2,chisq_marg_nvar2_sub, fill_value="extrapolate")
    
    # construction functions to find 1 and 2 sigma limits on nvar1 and nvar2 (from chi-sq min). 
    def func_nvar1_1sigma(x): return (func_chisq_1D_nvar1(x) - 0.99)
    def func_nvar1_2sigma(x): return (func_chisq_1D_nvar1(x) - 3.84)
    def func_nvar2_1sigma(x): return (func_chisq_1D_nvar2(x) - 0.99)
    def func_nvar2_2sigma(x): return (func_chisq_1D_nvar2(x) - 3.84)
        
    # guesses for the locations of the solutions in 1D [change with energy]:
    nvar1_min_1 = {}
    nvar1_max_1 = {}
    nvar1_min_2 = {}
    nvar1_max_2 = {}

    nvar1_min_1[100] = searchlimits[Energy][nvar1][0]
    nvar1_max_1[100] = searchlimits[Energy][nvar1][1]
    nvar1_min_2[100] = searchlimits[Energy][nvar1][0]
    nvar1_max_2[100] = searchlimits[Energy][nvar1][1]

    nvar2_min_1 = {}
    nvar2_max_1 = {}
    nvar2_min_2 = {}
    nvar2_max_2 = {}
    
    nvar2_min_1[100] = searchlimits[Energy][nvar2][0]
    nvar2_max_1[100] = searchlimits[Energy][nvar2][1]
    nvar2_min_2[100] = searchlimits[Energy][nvar2][0]
    nvar2_max_2[100] = searchlimits[Energy][nvar2][1]

    CL_threshold = 3.84  # 95% CL
    allowed = np.where(chisq_marg_nvar1_sub <= CL_threshold)[0]
    x_limits = x1[allowed]
    x_lower, x_upper = x_limits[0], x_limits[-1]
    print(f"95% CL for c3: {x_lower:.3f} to {x_upper:.3f} (c3)")
    allowed = np.where(chisq_marg_nvar2_sub <= CL_threshold)[0]
    x_limits = x2[allowed]
    x_lower, x_upper = x_limits[0], x_limits[-1]
    print(f"95% CL for d4: {x_lower:.3f} to {x_upper:.3f} (d4)")
    
    
    # calculate and print out the solutions:
    #print(variables[nvar1] + '@68% CL:', fsolve(func_nvar1_1sigma, [nvar1_min_1[Energy], nvar1_max_1[Energy]]))
    print(variables[nvar1] + '@95% CL:', fsolve(func_nvar1_2sigma, [nvar1_min_2[Energy], nvar1_max_2[Energy]]))
    #print(variables[nvar2] + '@68% CL:', fsolve(func_nvar2_1sigma, [nvar2_min_1[Energy], nvar2_max_1[Energy]]))
    print(variables[nvar2] + '@95% CL:', fsolve(func_nvar2_2sigma, [nvar2_min_2[Energy], nvar2_max_2[Energy]]))

    # TEST FUNCTIONS TO SOLVE HERE:
    plt.clf()
    x2 = np.linspace(ymin, ymax, nbins) # nvar2
    y2 = func_nvar2_2sigma(x2)
    y1 = func_nvar2_1sigma(x2)
    plt.plot(x2, y1)
    plt.plot(x2, y2)
    plt.axhline(0, color='k', linestyle='--')
    plt.savefig(plot_dir + output + 'test_d4.pdf', bbox_inches='tight')
    plt.clf()
    x2 = np.linspace(xmin, xmax, nbins) # nvar2
    y2 = func_nvar1_2sigma(x2)
    y1 = func_nvar1_1sigma(x2)
    plt.plot(x2, y1)
    plt.plot(x2, y2)
    plt.axhline(0, color='k', linestyle='--')
    plt.savefig(plot_dir + output + 'test_c3.pdf', bbox_inches='tight')
    plt.clf()
    # END OF TEST FUNCTIONS TO SOLVE

    # construct the axes for the plot
    # no need to modify this if you just need one plot
    gs = gridspec.GridSpec(4, 4)
    if figext == '':
        fig = plt.figure()
    else:
        fig = figext
    if axext == '':
        ax = fig.add_subplot(111)
    else:
        ax=axext
    ax.grid(False)
    ax.set_title(plottitle, fontsize=10)
    
    # set the ticks, labels and limits etc.
    xlab = '$' + variables_latex[nvar1] + '$'
    ylab = '$' + variables_latex[nvar2] + '$'
    if setylabel == True:
        ax.set_ylabel(ylab, fontsize=labelsize)
    if setxlabel == True:
        ax.set_xlabel(xlab, fontsize=labelsize)
    # plot the contours:
    ax.plot(0,0,marker='*',ms=starsize, color='black')
    
    #cont = ax.contour(X, Y, chisq_sub, contours, extend='max', colors=('black', 'red'), linestyles=('-','--'))
    #labels = ['$1\\sigma$', '$2\\sigma$']
    #for i in range(len(labels)):
    #    cont.collections[i].set_label(labels[i])
    plt.xlim([xmin,xmax])
    plt.ylim([ymin,ymax])
    cont = ax.contour(
    X, Y, chisq_sub, contours,
    extend='max',
    colors=('black', 'red'),
    linestyles=('-', '--')
        )
    labels = ['$1\\sigma$', '$2\\sigma$']

    # Create a dictionary mapping each level to its label
    label_dict = {level: label for level, label in zip(cont.levels, labels)}
    
    # Add the legend labels via ax.clabel with a formatter
    ax.clabel(cont, fmt=label_dict, manual=False)  # Remove manual=... for automatic labeling

    
    if smtext == True:
        ax.text(0.53, 0.40,"SM", transform=ax.transAxes)
    if normalbar == True:
        plt.colorbar(cont)
    #handles, labels = cs.legend_elements()

    # after you’ve done your contour call…
    black_line = mlines.Line2D([], [], color='black', linestyle='-',
                            label='$1\\sigma$')
    red_line   = mlines.Line2D([], [], color='red',   linestyle='--',
                           label='$2\\sigma$')

    ax.legend(handles=[black_line, red_line],
          loc="upper right", frameon=False, prop={'size':8})
        
    #ax.legend()
    #ax.legend(loc="upper right", numpoints=1, frameon=False, prop={'size':8}, handles=[cont, cont2])
    ax.yaxis.set_major_locator(ticker.AutoLocator())
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.xaxis.set_major_locator(ticker.AutoLocator())
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())

    if savefig == True:
        # save the figure
        print('saving the figure')
        # save the figure in PDF format
        infile = output + '.dat'
        print('---')
        print('output in', infile.replace('.dat','.pdf'))
        plt.savefig(plot_dir + infile.replace('.dat','.pdf'), bbox_inches='tight')
        plt.close(fig)
        
    return cont, X, Y, chisq_sub


def contour_pvalue_only(procname, plotname, plottitle, fit_coeffs_xsec, fit_coeffs_eff, sigma_bkg, var1, var2, plotlimits, searchlimits, deltac3=-1, axext='', figext='', smtext=True, starsize=15, setxlabel=True, setylabel=True, nbins=400, savefig=True,variables=variables, variables_latex=variables_latex, labelsize=20, normalbar=True, contours=np.arange(0, 10, 0.5),norm_to_zeroth=True, lumi=Luminosity):
    output = procname + '_' + plotname + '_' + var1 + '_' + var2
    print('Plotting', output)
    nvar1 = [key for key, value in variables.items() if value == var1][0]
    nvar2 = [key for key, value in variables.items() if value == var2][0]
    nvar3 = [key for key, value in variables.items() if value != var1 and value != var2][0]
    nvar4 = [key for key, value in variables.items() if value != var1 and value != var2][1]

    #print(var1, var2)
    print('nvar1, nvar2=', nvar1, nvar2)
    #print(nvar3, nvar4)
   
        
    # set the limits on the x and y axes if required below:
    ymin = plotlimits[Energy][nvar2][0]
    ymax = plotlimits[Energy][nvar2][1]
    xmin = plotlimits[Energy][nvar1][0]
    xmax = plotlimits[Energy][nvar1][1]
  
    ctexts = []
    cvartexts = []
    for i in range(0, len(variables.keys())):
        cvartexts.append(variables[i])
    fstr = 'partial(func_t_CX, procname=Process)'
    global func_CX_partial
    func_CX_partial = eval(fstr)
    #print('func_CX_partial=', func_CX_partial)
    #print(fit_coeffs)
    global fit_coeffs_g_xsec
    fit_coeffs_g_xsec = fit_coeffs_xsec
    global fit_coeffs_g_eff
    fit_coeffs_g_eff = fit_coeffs_eff
    #print(cvartexts[0], cvartexts[1])

    # functions for xsec and significance:
    print(cvartexts)
    # 0     1      2      3     
    #['c3', 'ct2', 'ct3', 'd4']
    print(cvartexts[nvar3], cvartexts[nvar4], cvartexts[nvar1], cvartexts[nvar2])
    feval_xsec = 'func_CX_partial(' + cvartexts[nvar3] +'=0,' + cvartexts[nvar4] + '=0,' + cvartexts[nvar1] + '=x1,' + cvartexts[nvar2] + '=x2,coeffs=fit_coeffs_g_xsec)'
    feval_eff = 'func_CX_partial(' + cvartexts[nvar3] +'=0,' + cvartexts[nvar4] + '=0,'  + cvartexts[nvar1] + '=x1,' + cvartexts[nvar2] + '=x2,coeffs=fit_coeffs_g_eff)'
      
    func_fin = lambda x1, x2: eval(feval_xsec) * sig_factors * eval(feval_eff) * lumi * 1000. / math.sqrt(sigma_bkg * lumi + Systematics**2 * (sigma_bkg)**2 * lumi**2)

    pfunc_fin = lambda x1, x2: (1 - 2 * scipy.special.ndtr(-(eval(feval_xsec) * sig_factors * eval(feval_eff) * lumi * 1000. / math.sqrt(sigma_bkg * lumi + Systematics**2 * (sigma_bkg)**2 * lumi**2))))

    # SM significance: 
    feval_xsec_sm = 'func_CX_partial(' + cvartexts[0] +'=0,' + cvartexts[1] + '=0,' + cvartexts[2] + '=0,' + cvartexts[3] + '=0,coeffs=fit_coeffs_g_xsec)'
    feval_eff_sm = 'func_CX_partial(' + cvartexts[0] +'=0,' + cvartexts[1] + '=0,' + cvartexts[2] + '=0,' + cvartexts[3] + '=0,coeffs=fit_coeffs_g_eff)'
    
    sm_signif= eval(feval_xsec_sm) * sig_factors * eval(feval_eff_sm) * 1000. * lumi / math.sqrt(sigma_bkg * lumi + Systematics**2 * (sigma_bkg)**2 * lumi**2)

    # if SM is the "null" hypothesis:
    # SM number of events:
    S_SM = eval(feval_xsec_sm) * sig_factors * eval(feval_eff_sm) * 1000. * lumi
    # SM total uncertainty, including the background uncertainty:
    delta_SM =  math.sqrt(S_SM + sigma_bkg * lumi + Systematics**2 * (sigma_bkg)**2 * lumi**2) 
    # {c_i} number of events in 2D model:
    S_i_2D = lambda x1, x2: eval(feval_xsec) * sig_factors * eval(feval_eff) * lumi * 1000.
    # significance versus the SM in the 2D mode: 
    func_fin_SM_2D = lambda x1, x2: np.power( (S_SM - S_i_2D(x1, x2))/delta_SM,2)

    # p-value in the 2D model:
    #pfunc_fin_SM_2D = lambda x1, x2:  1/(np.sqrt(2.*np.pi)*delta_SM)*np.exp(-func_fin_SM_2D(x1, x2)/2)
    # print('pfunc_fin_SM_2D(0,0)=',pfunc_fin_SM_2D(0,0))
    
    print("sigma_sig before anal. [fb]=", eval(feval_xsec_sm)*1000*sig_factors)
    print("analysis eff. on signal=", eval(feval_eff_sm))
    print("sigma_bkg after anal. [fb]=", sigma_bkg)
    print("sigma sig SM after anal. [fb]=",eval(feval_xsec_sm) * sig_factors * eval(feval_eff_sm) * 1000)
    print("N(bkg)@lumi=", sigma_bkg * lumi)
    print("N(sig SM)@lumi=", eval(feval_xsec_sm) * sig_factors * eval(feval_eff_sm) * lumi * 1000.) 
    print("SM significance=", sm_signif)
    

    
    x = np.linspace(xmin, xmax, nbins)
    y = np.linspace(ymin, ymax, nbins)
    dx = x[1] - x[0]
    dy = y[1] - y[0]
    X, Y = np.meshgrid(x,y)

    # ---- New Bayesian & Frequentist interval calculation clearly added here ----
    chisq_prior = (X)**2 / deltac3**2 if deltac3 > 0 else 0.0
        
    # Combine chi-squared with prior
    chisq_total = func_fin_SM_2D(X,Y) + chisq_prior
    chisq_sub = chisq_total - np.amin(chisq_total)
    
    # Bayesian posterior
    posterior = np.exp(-chisq_total / 2)
    dx = x[1] - x[0]
    dy = y[1] - y[0]
    posterior /= np.sum(posterior) * dx * dy  # normalized posterior
    
    # Bayesian marginalization over each parameter:
    posterior_nvar1 = np.sum(posterior, axis=0) * dy  # marginalized over Y (nvar2)
    posterior_nvar2 = np.sum(posterior, axis=1) * dx  # marginalized over X (nvar1)
    
    cdf_nvar1 = np.cumsum(posterior_nvar1) * dx
    cdf_nvar2 = np.cumsum(posterior_nvar2) * dy
    
    cdf_nvar1 /= cdf_nvar1[-1]
    cdf_nvar2 /= cdf_nvar2[-1]
    
    # 95% Bayesian credible intervals (central 95%)
    nvar1_low95 = np.interp(0.025, cdf_nvar1, x)
    nvar1_high95 = np.interp(0.975, cdf_nvar1, x)
    
    nvar2_low95 = np.interp(0.025, cdf_nvar2, y)
    nvar2_high95 = np.interp(0.975, cdf_nvar2, y)
    
    print("Bayesian 95% credible interval for", variables[nvar1], f": {nvar1_low95:.3f} to {nvar1_high95:.3f}")
    print("Bayesian 95% credible interval for", variables[nvar2], f": {nvar2_low95:.3f} to {nvar2_high95:.3f}")

    # --- Frequentist profiling clearly included here ---
    chisq_profile_nvar1 = np.min(chisq_sub, axis=0)  # profile over nvar2 (y)
    chisq_profile_nvar2 = np.min(chisq_sub, axis=1)  # profile over nvar1 (x)
    
    # Frequentist 95% confidence intervals (Δχ²=3.84 for 1 parameter)
    allowed_nvar1 = x[chisq_profile_nvar1 <= 3.84]
    allowed_nvar2 = y[chisq_profile_nvar2 <= 3.84]

    freq_nvar1_low, freq_nvar1_high = allowed_nvar1[0], allowed_nvar1[-1]
    freq_nvar2_low, freq_nvar2_high = allowed_nvar2[0], allowed_nvar2[-1]

    print("Frequentist (profile) 95% CL for", variables[nvar1], f": {freq_nvar1_low:.3f} \t {freq_nvar1_high:.3f}")
    print("Frequentist (profile) 95% CL for", variables[nvar2], f": {freq_nvar2_low:.3f} \t {freq_nvar2_high:.3f}")
    
    # write frequentist results to files:
    filewrite_frequentist_c3 = ConstraintsDir + output + 'frequentist_c3.out'
    filewrite_frequentist_d4 = ConstraintsDir + output + 'frequentist_d4.out'
    with open(filewrite_frequentist_c3,'w') as f:
        f.write(str(f"{freq_nvar1_low:.3f} \t {freq_nvar1_high:.3f}"))
    with open(filewrite_frequentist_d4,'w') as f:
        f.write(str(f"{freq_nvar2_low:.3f} \t {freq_nvar2_high:.3f}"))
    
    # interpolate the 1D functions: 
    func_chisq_1D_nvar1 = interp1d(x,chisq_profile_nvar1, fill_value="extrapolate")
    func_chisq_1D_nvar2 = interp1d(y,chisq_profile_nvar2, fill_value="extrapolate")
    
    # construction functions to find 1 and 2 sigma limits on nvar1 and nvar2 (from chi-sq min). 
    def func_nvar1_1sigma(x): return (func_chisq_1D_nvar1(x) - 0.99)
    def func_nvar1_2sigma(x): return (func_chisq_1D_nvar1(x) - 3.84)
    def func_nvar2_1sigma(x): return (func_chisq_1D_nvar2(x) - 0.99)
    def func_nvar2_2sigma(x): return (func_chisq_1D_nvar2(x) - 3.84)
    
    # TEST FUNCTIONS TO SOLVE HERE:
    plt.clf()
    x2 = np.linspace(ymin, ymax, nbins) # nvar2
    y2 = func_nvar2_2sigma(y)
    y1 = func_nvar2_1sigma(y)
    plt.plot(y, y1)
    plt.plot(y, y2)
    plt.axhline(0, color='k', linestyle='--')
    plt.savefig(plot_dir + output + 'test_d4.pdf', bbox_inches='tight')
    plt.clf()
    x2 = np.linspace(xmin, xmax, nbins) # nvar2
    y2 = func_nvar1_2sigma(x)
    y1 = func_nvar1_1sigma(x)
    plt.plot(x, y1)
    plt.plot(x, y2)
    plt.axhline(0, color='k', linestyle='--')
    plt.savefig(plot_dir + output + 'test_c3.pdf', bbox_inches='tight')
    plt.clf()
    # END OF TEST FUNCTIONS TO SOLVE

    # construct the axes for the plot
    # no need to modify this if you just need one plot
    gs = gridspec.GridSpec(4, 4)
    if figext == '':
        fig = plt.figure()
    else:
        fig = figext
    if axext == '':
        ax = fig.add_subplot(111)
    else:
        ax=axext
    ax.grid(False)
    ax.set_title(plottitle, fontsize=10)
    
    # set the ticks, labels and limits etc.
    xlab = '$' + variables_latex[nvar1] + '$'
    ylab = '$' + variables_latex[nvar2] + '$'
    if setylabel == True:
        ax.set_ylabel(ylab, fontsize=labelsize)
    if setxlabel == True:
        ax.set_xlabel(xlab, fontsize=labelsize)
    # plot the contours:
    ax.plot(0,0,marker='*',ms=starsize, color='black')
    
    #cont = ax.contour(X, Y, chisq_sub, contours, extend='max', colors=('black', 'red'), linestyles=('-','--'))
    #labels = ['$1\\sigma$', '$2\\sigma$']
    #for i in range(len(labels)):
    #    cont.collections[i].set_label(labels[i])
    plt.xlim([xmin,xmax])
    plt.ylim([ymin,ymax])
    cont = ax.contour(
    X, Y, chisq_sub, contours,
    extend='max',
    colors=('black', 'red'),
    linestyles=('-', '--')
        )
    labels = ['$1\\sigma$', '$2\\sigma$']

    # Create a dictionary mapping each level to its label
    label_dict = {level: label for level, label in zip(cont.levels, labels)}
    
    # Add the legend labels via ax.clabel with a formatter
    ax.clabel(cont, fmt=label_dict, manual=False)  # Remove manual=... for automatic labeling

    
    if smtext == True:
        ax.text(0.53, 0.40,"SM", transform=ax.transAxes)
    if normalbar == True:
        plt.colorbar(cont)
    #handles, labels = cs.legend_elements()

    # after you’ve done your contour call…
    black_line = mlines.Line2D([], [], color='black', linestyle='-',
                            label='$1\\sigma$')
    red_line   = mlines.Line2D([], [], color='red',   linestyle='--',
                           label='$2\\sigma$')

    ax.legend(handles=[black_line, red_line],
          loc="upper right", frameon=False, prop={'size':8})
        
    #ax.legend()
    #ax.legend(loc="upper right", numpoints=1, frameon=False, prop={'size':8}, handles=[cont, cont2])
    ax.yaxis.set_major_locator(ticker.AutoLocator())
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.xaxis.set_major_locator(ticker.AutoLocator())
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())

    if savefig == True:
        # save the figure
        print('saving the figure')
        # save the figure in PDF format
        infile = output + '.dat'
        print('---')
        print('output in', infile.replace('.dat','.pdf'))
        plt.savefig(plot_dir + infile.replace('.dat','.pdf'), bbox_inches='tight')
        plt.close(fig)
        
    return cont, X, Y, chisq_sub


# save the data:
def save_data(data, filename):
    with open(filename,'wb') as f:
        pickle.dump(data,f)

# load the data:
def load_data(filename):
    with open(filename, 'rb') as f:
        data = pickle.load(filename)
    return data


################################c#########################
# RUN THE CODE HERE                                     # 
#########################################################


############################################
# GENERATE MG5 LHE FILES FOR SIGNAL:       #
############################################

if MODEL != 'HEFT4C3D4' and MODEL != 'C3D4ONLY':
    # reduced set for Global HHH (100 TeV runs):
    # [c3, d4, ct2, ct3 ] -> 4 couplings
    couplings_min = [-5.0, -5.0, -0.1, -0.1]
    couplings_max = [5.0, 5.0, 0.1, 0.1]
    if MODEL == 'HEFT6':
        couplings_min = couplings_min + [-0.1]
        couplings_max = couplings_max + [0.1]
elif MODEL == 'HEFT4C3D4' or MODEL== 'C3D4ONLY':
    couplings_min = [-5.0, -5.0, 0.0, 0.0]
    couplings_max = [5.0, 5.0, 0.0, 0.0]

# generate random coupling arrays:
randseed = 999
CouplingsArray_R, CouplingsArrayF_R = gen_coupbdasarray_dim_rand_range(couplings_min, couplings_max, Nruns, randseed)

# additional set:
Nadditional=560
couplings_min = [-0.5, -0.5, -0.001, -0.001]
couplings_max = [0.5, 0.5, 0.001, 0.001]
if MODEL == 'HEFT6':
    couplings_min = couplings_min + [-0.1]
    couplings_max = couplings_max + [0.1]
CouplingsArray_R_add, CouplingsArrayF_R_add = gen_coupbdasarray_dim_rand_range(couplings_min, couplings_max, Nadditional, randseed+31)

# additional set 2:
Nadditional2=200
couplings_min = [-5.0, -50, 0, 0]
couplings_max = [5.0, 50, 0, 0]
if MODEL == 'HEFT6':
    couplings_min = couplings_min + [-0.1]
    couplings_max = couplings_max + [0.1]
CouplingsArray_R_add2, CouplingsArrayF_R_add2 = gen_coupbdasarray_dim_rand_range(couplings_min, couplings_max, Nadditional2, randseed+29)

# additional set 3:
Nadditional3=300
couplings_min = [-100.0, -100, 0, 0]
couplings_max = [100.0, 100, 0, 0]
if MODEL == 'HEFT6':
    couplings_min = couplings_min + [-0.1]
    couplings_max = couplings_max + [0.1]
CouplingsArray_R_add3, CouplingsArrayF_R_add3 = gen_coupbdasarray_dim_rand_range(couplings_min, couplings_max, Nadditional3, randseed+27)

# additional set 4:
Nadditional4=290
couplings_min = [-10.0, -100, -2, -2]
couplings_max = [10.0, 100, 2, 2]
if MODEL == 'HEFT6':
    couplings_min = couplings_min + [-2]
    couplings_max = couplings_max + [2]
CouplingsArray_R_add4, CouplingsArrayF_R_add4 = gen_coupbdasarray_dim_rand_range(couplings_min, couplings_max, Nadditional4, randseed+33)

# additional set 5:
Nadditional5=500
couplings_min = [-20.0, -100, -4, -4]
couplings_max = [20.0, 100, 4, 4]
if MODEL == 'HEFT6':
    couplings_min = couplings_min + [-5]
    couplings_max = couplings_max + [5]
CouplingsArray_R_add5, CouplingsArrayF_R_add5 = gen_coupbdasarray_dim_rand_range(couplings_min, couplings_max, Nadditional5, randseed+99)

# additional set 6:
Nadditional6=100
couplings_min = [-40.0, -100, 0, 0]
couplings_max = [40.0, 100, 0, 0]
if MODEL == 'HEFT6':
    couplings_min = couplings_min + [-0.1]
    couplings_max = couplings_max + [0.1]
CouplingsArray_R_add6, CouplingsArrayF_R_add6 = gen_coupbdasarray_dim_rand_range(couplings_min, couplings_max, Nadditional6, randseed+4)

# additional set 7:
#Nadditional7=0
#CouplingsArray_R_add7, CouplingsArrayF_R_add7 = [], []
Nadditional7=200
couplings_min = [-50.0, -800, -2, -2] 
couplings_max = [50.0, 800, 2, 2]
if MODEL == 'HEFT6':
    couplings_min = couplings_min + [-0.1]
    couplings_max = couplings_max + [0.1]
CouplingsArray_R_add7, CouplingsArrayF_R_add7 = gen_coupbdasarray_dim_rand_range(couplings_min, couplings_max, Nadditional7, randseed+3)

# additional set 8: 
#Nadditional8=0
#CouplingsArray_R_add8, CouplingsArrayF_R_add8 = [], []
Nadditional8=200
couplings_min = [-20.0, -600, -2, -2]
couplings_max = [20.0, 600, 2, 2]
if MODEL == 'HEFT6':
    couplings_min = couplings_min + [-1.0]
    couplings_max = couplings_max + [1.0]
CouplingsArray_R_add8, CouplingsArrayF_R_add8 = gen_coupbdasarray_dim_rand_range(couplings_min, couplings_max, Nadditional8, randseed+2)

# for testing: reset some to zero:
#CouplingsArray_R_add, CouplingsArrayF_R_add = [], []
#CouplingsArray_R_add2, CouplingsArrayF_R_add2 = [], []
#CouplingsArray_R_add3, CouplingsArrayF_R_add3 = [], []
#CouplingsArray_R_add4, CouplingsArrayF_R_add4 = [], []
#CouplingsArray_R_add5, CouplingsArrayF_R_add5 = [], []
#CouplingsArray_R_add6, CouplingsArrayF_R_add6 = [], []
#CouplingsArray_R_add7, CouplingsArrayF_R_add7 = [], []
#CouplingsArray_R_add8, CouplingsArrayF_R_add8 = [], []
#Nadditional, Nadditional2, Nadditional3, Nadditional4,
#Nadditional5, Nadditional6, Nadditional7, Nadditional8 = 0,0,0,0

# concatenate
CouplingsArray_R+=CouplingsArray_R_add+CouplingsArray_R_add2+CouplingsArray_R_add3+CouplingsArray_R_add4+CouplingsArray_R_add5+CouplingsArray_R_add6+CouplingsArray_R_add7+CouplingsArray_R_add8
CouplingsArrayF_R+=CouplingsArrayF_R_add+CouplingsArrayF_R_add2+CouplingsArrayF_R_add3+CouplingsArrayF_R_add4+CouplingsArrayF_R_add5+CouplingsArrayF_R_add6+CouplingsArrayF_R_add7+CouplingsArrayF_R_add8

# Launch MG5 event generation
nevents=1
drive_mg_proc(RunNum, MGLocation, ProcLocations[Process], Process, CouplingsArray_R, nevents, Nruns+Nadditional+Nadditional2+Nadditional3+Nadditional4+Nadditional5+Nadditional6+Nadditional7+Nadditional8, ecm=Energy)

###################################
# PERFORM THE FIT OR READ THE FIT #
###################################

# read the generated MG5 files:
if DoFit is True:
    print('reading in generated files')
    print('CouplingsArray_R=',CouplingsArray_R)
    X, Z, ZERR, XSEC = read_files(RunNum, MGLocation, ProcLocations[Process], Process, CouplingsArray_R, Nruns+Nadditional+Nadditional2+Nadditional3+Nadditional4+Nadditional5+Nadditional6+Nadditional7+Nadditional8)
    print(X)
else:
    print('Not reading in files, will read fit!')


# generate the list of initial guesses:
p0_i = []
p0_iE = []
for i in range(0,NCoeffs[Process]):
    p0_i.append(0.01)
    p0_iE.append(0.1)
# get the partial function with the process fixed:
func_CX_proc = partial(func_CX, procname=Process)
# perform the fit:

if DoFit is True:
    popt[Process], pcov[Process] = curve_fit(func_CX_proc, tuple(X) , Z, sigma=ZERR, method='lm', maxfev=2000, p0=p0_i)
    saveFit(popt[Process], pcov[Process], Process, RunNum)
    # test the fit:
    test_fit(RunNum, MGLocation, ProcLocations[Process], Process, CouplingsArray_R,  Nruns+Nadditional+Nadditional2+Nadditional3+Nadditional4+Nadditional5+Nadditional6+Nadditional7+Nadditional8, popt[Process])
else:
    popt[Process], pcov[Process] = readFit(Process, RunNum)

if debug:
    print('fitted parameters:')
    print(popt[Process])
    
if RunHerwig is False:
    print('Fit coefficients for MODEL=', MODEL, '=',  popt[Process]/popt[Process][-1])
    print('Errors=', np.sqrt(np.diag(pcov[Process]))/popt[Process][-1])
    print('RunHerwig is False: Not running Herwig or analysis, exiting')
    exit()

####################################
# RUN HERWIG ON LHE FILES          #
# AND PERFORM THE ANALYSIS         #
# FIT THE EFFICIENCY               # 
####################################


print('Running Herwig on generated MG5 LHEs')
run_herwig_proc_parallel(RunNum, MGLocation, HerwigOutputDirectory, ProcLocations[Process], Process, CouplingsArray_R, nevents,  Nruns+Nadditional+Nadditional2+Nadditional3+Nadditional4+Nadditional5+Nadditional6, ecm=Energy)

print('Running analysis on signal and background')
XE, ZE, EFFICIENCY, EFFICIENCY_BKG = run_analysis_proc(RunNum, MGLocation, HerwigOutputDirectory, ProcLocations[Process], Process, CouplingsArray_R, nevents,  Nruns+Nadditional+Nadditional2+Nadditional3+Nadditional4+Nadditional5+Nadditional6, ecm=Energy)
# to fix the issue of not reading it the first time ReRunAnalysis is True
#if ReRunAnalysis is True:
#    ReRunAnalysis = False
#    XE, ZE, EFFICIENCY, EFFICIENCY_BKG = run_analysis_proc(RunNum, MGLocation, HerwigOutputDirectory, ProcLocations[Process], Process, CouplingsArray_R, nevents,  Nruns+Nadditional+Nadditional2+Nadditional3+Nadditional4+Nadditional5+Nadditional6, ecm=Energy)

print(ZE)
popt_eff, pcov = curve_fit(func_CX_proc, tuple(XE), ZE, method='lm', maxfev=10000, p0=p0_iE)
test_fit_analysis(RunNum, MGLocation, ProcLocations[Process], Process, CouplingsArray_R,  Nruns+Nadditional+Nadditional2+Nadditional3+Nadditional4+Nadditional5+Nadditional6, popt_eff)

# get the SM efficiency:
analysisInputfile = './Herwig/events/HW-8_SM_6b.root' 
print('running the analysis', ExecutableSmear[Energy], 'on the input file', analysisInputfile)
analysiscommand = ExecutableSmear[Energy] + ' ' + analysisInputfile
print('Launching:', analysiscommand)
p = subprocess.Popen(analysiscommand, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, cwd='.')
for line in iter(p.stdout.readline, b''):
    print('\t\t', line, end=' ')
out, err = p.communicate()
analysisOutputfile = analysisInputfile.replace('.root', '.smear' + smearing_tag + '.dat')
if os.path.exists(analysisOutputfile)is False:
    print('File', analysisOutputfile, 'does not exist!')
    exit()
else:
    print('File', analysisOutputfile, ' exists, reading results')
    zgrepcommand = 'cat ' + analysisOutputfile
    p = subprocess.Popen(zgrepcommand, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, cwd='.')
    for line in iter(p.stdout.readline, b''):
        SM_efficiency = float(line.split()[0])
        
# print the SM
print('SM RESULTS:')
print('Signal cross section BEFORE cuts (no b-tagging/BRs/k-factors)=', xsS)
print('Signal cross section BEFORE cuts (WITH b-tagging/BRs/k-factors)=', xsS*sig_factors)
print('Signal cut efficiency=', SM_efficiency)
sigma_SM_after_cuts = xsS * sig_factors * SM_efficiency
print('Signal cross section AFTRER cuts=', sigma_SM_after_cuts)
print('NSM(EVENTS)=', Luminosity*sigma_SM_after_cuts)

# Print the backgrounds
sigma_bkg = 0
print('Background cross sections BEFORE cuts:') 
for bkg in Backgrounds:
    print(bkg, 'sigma=', Backgrounds_xsec[(Energy, bkg)])
print('Background cut efficiency:') 
for bkg in Backgrounds:
    print(bkg, 'eff=', EFFICIENCY_BKG[bkg])
print('Background cross sections AFTER cuts:') 
for bkg in Backgrounds:
    sigma_bkg = sigma_bkg + EFFICIENCY_BKG[bkg] * Backgrounds_xsec[(Energy, bkg)]*bkg_factors
    print(bkg, 'sigma=', EFFICIENCY_BKG[bkg] * Backgrounds_xsec[(Energy, bkg)])
print('Background EXPECTED NUMBER OF EVENTS AFTER cuts:') 
for bkg in Backgrounds:
    print(bkg, 'N(EVENTS)=', bkg_factors*Luminosity*EFFICIENCY_BKG[bkg] * Backgrounds_xsec[(Energy, bkg)])
print('sigma_bkg total (fb) = ', sigma_bkg)

print("EXPECTED SM SIGNIFICANCE (CUTS)=", Luminosity*sigma_SM_after_cuts/np.sqrt(sigma_bkg*Luminosity + (Systematics*sigma_bkg*Luminosity)**2))

########################################
# XGBOOST ANALYSIS HERE
########################################
# do the training on the SM
training_seed = 12345
if DoTraining is True:
    trained_model = train_xgboost(signal_SM_file, Backgrounds, Background_files, Backgrounds_xsec, xsS, initial_S_SM, sig_factors, initial_B, idB, bkg_factors, Luminosity, Energy, training_seed)
    trained_model_file = 'trained_model' + str(RunNum) + smearing_tag + '.json'
    save_model(trained_model, trained_model_file)
else:
    trained_model_file = 'trained_model' + str(RunNum) + smearing_tag + '.json'
    trained_model = load_model(trained_model_file)
    
# apply the model on the SM (testing):
apply_xgboost(trained_model, signal_SM_file, Backgrounds, Background_files, Backgrounds_xsec, xsS, initial_S_SM, sig_factors, initial_B, idB, bkg_factors, Luminosity, Energy, training_seed)
time.sleep(10)
# apply to all points, get the efficiencies for signal and backgrounds
print('Running XGBOOST on all points')
XE, ZE, EFFICIENCY, EFFICIENCY_BKG = run_analysis_xgboost(RunNum, MGLocation, HerwigOutputDirectory, ProcLocations[Process], Process, CouplingsArray_R, nevents,  Nruns+Nadditional+Nadditional2+Nadditional3+Nadditional4+Nadditional5+Nadditional6, trained_model_file, Backgrounds, Background_files, Backgrounds_xsec, xsS, initial_S, sig_factors, initial_B, idB, bkg_factors, Luminosity, Energy, training_seed, ecm=Energy)
popt_eff_XGBOOST, pcov_XGBOOST = curve_fit(func_CX_proc, tuple(XE), ZE, method='lm', maxfev=1000000, p0=p0_iE)
# calculate the background cross section after the xgboost analysis: 
sigma_bkg_xgboost = 0
for bkg in Backgrounds:
    sigma_bkg_xgboost = sigma_bkg_xgboost + EFFICIENCY_BKG[bkg] * Backgrounds_xsec[(Energy, bkg)]*bkg_factors

########################################
# PLOTTING STARTS HERE                 #
########################################
    
# Plot "correlation" plot of the cross section
if MODEL != 'C3D4ONLY' and MODEL !='HEFT3':
    correlation_plot(Process, 'xsec'+ str(Energy), popt[Process], variables, plottitle='$\\sigma(gg\\rightarrow hhh)$@' + str(Energy) + ' TeV, normalized to SM value')

# plot 1D plots of the variation of the cross section with coefficient
oned_xsec(Process, 'xsec' + str(Energy), r'$\sigma(gg\rightarrow hhh)$@' + str(Energy) + ' TeV, normalized to SM value', popt[Process], 'c3',[-5.0, 5.0], [0.5, 10])
oned_xsec(Process, 'xsec' + str(Energy), r'$\sigma(gg\rightarrow hhh)$@' + str(Energy) + ' TeV, normalized to SM value', popt[Process], 'd4',[-40.0, 40.0], [0.5, 10])

if MODEL != 'C3D4ONLY' and MODEL !='HEFT3':
    oned_xsec(Process, 'xsec' + str(Energy), r'$\sigma(gg\rightarrow hhh)$@' + str(Energy) + ' TeV, normalized to SM value', popt[Process], 'ct2',[-1.0, 1.0], [0.5, 25.0])
    oned_xsec(Process, 'xsec' + str(Energy), r'$\sigma(gg\rightarrow hhh)$@' + str(Energy) + ' TeV, normalized to SM value', popt[Process], 'ct3',[-1.0, 1.0], [0.5, 25.0])

# plot the "correlation plot" of the efficiency
if MODEL != 'C3D4ONLY' and MODEL !='HEFT3':
    correlation_plot(Process, 'eff'+ str(Energy), popt_eff, variables, plottitle='$\\epsilon(gg\\rightarrow hhh)$@' + str(Energy) + ' TeV', contours=np.arange(0.005, 0.02, 0.0005), norm_to_zeroth=False)

# for the XGBOOST case:
if MODEL != 'C3D4ONLY' and MODEL !='HEFT3':
    correlation_plot(Process, 'eff_XGBOOST'+ str(Energy), popt_eff_XGBOOST, variables, plottitle='$\\epsilon_\\mathrm{XG}(gg\\rightarrow hhh)$@' + str(Energy) + ' TeV', contours=np.arange(0.005, 0.02, 0.0005), norm_to_zeroth=False)


#########################################
# p-value contours and calculations
#########################################
nbinsdist=5000

# limits on the plots (exclusion)
plotlimits = {}
plotlimits[100] = {}
searchlimits = {}
searchlimits[100] = {}
if Systematics == 0.0:
    plotlimits[100][0] = [-5.0, 5.0]
    plotlimits[100][1] = [-1.0, 1.0]
    plotlimits[100][2] = [-0.5, 0.5]
    plotlimits[100][3] = [-30.0, 40.0]
    # search limits for the exclusion
    searchlimits[100][0] = [-8.0, 8.0]
    searchlimits[100][1] = [-1.0, 1.0]
    searchlimits[100][2] = [-0.5, 0.5]
    searchlimits[100][3] = [-30.0, 40.0]    
else:
    plotlimits[100][0] = [-10.0, 12.0]
    plotlimits[100][1] = [-1.0, 1.0]
    plotlimits[100][2] = [-0.5, 0.5]
    plotlimits[100][3] = [-180.0, 110.0]
    searchlimits[100][0] = [-10.0, 12.0]
    searchlimits[100][1] = [-1.0, 1.0]
    searchlimits[100][2] = [-0.5, 0.5]
    searchlimits[100][3] = [-180.0, 90.0]

if EnergyToRescale == 13 and DoRescaling is True:
    plotlimits[100][0] = [-20.0, 20.0]
    plotlimits[100][1] = [-1.0, 1.0]
    plotlimits[100][2] = [-0.5, 0.5]
    plotlimits[100][3] = [-100.0, 100.0]
    # search limits for the exclusion
    searchlimits[100][0] = [-20.0, 20.0]
    searchlimits[100][1] = [-1.0, 1.0]
    searchlimits[100][2] = [-0.5, 0.5]
    searchlimits[100][3] = [-100.0, 100.0]

if Luminosity < 1000:
    plotlimits[100][0] = [-10.0, 10.0]
    plotlimits[100][1] = [-1.0, 1.0]
    plotlimits[100][2] = [-0.5, 0.5]
    plotlimits[100][3] = [-60.0, 60.0]
    # search limits for the exclusion
    searchlimits[100][0] = [-8.0, 10.0]
    searchlimits[100][1] = [-1.0, 1.0]
    searchlimits[100][2] = [-0.5, 0.5]
    searchlimits[100][3] = [-60.0, 60.0]    
    

#contour_pvalue_ct3d4_marginalized(Process, 'pvalue'+ str(Energy) + '_L' + str(Luminosity), '$gg\\rightarrow hhh$@' + str(Energy) + ' TeV, L=' + str(Luminosity) + ' fb$^{-1}$, $\\alpha_\\mathrm{syst.} = ' + str(100*Systematics) +  '\%$', popt[Process], popt_eff, sigma_bkg, 'ct3', 'd4', plotlimits[Energy][2], plotlimits[Energy][3], contours=[2.278868566376729, 5.99], nbins=nbinsdist, normalbar=False)


# the tag for the output PDFs:
fulltag = str(Energy) + '_L' + str(Luminosity) + '_Syst' + str(Systematics) + '_pb' + str(btagging) + smearing_tag + '_' + MODEL + RESCALETAG + KFACTAG

# ct3, d4 (all others zero)
if MODEL != 'C3D4ONLY' and MODEL !='HEFT3':
    cont, X, Y, chisq_sub = contour_pvalue_only(Process, 'pvalue'+ fulltag, '$gg\\rightarrow hhh$@' + str(Energy) + ' TeV, L=' + str(Luminosity) + ' fb$^{-1}$, $\\mathcal{P}(b \\rightarrow b ) =' + str(btagging) + ' $' + ', $\\alpha_\\mathrm{syst.} = ' + str(100*Systematics) +  '\%$', popt[Process], popt_eff, sigma_bkg, 'ct3', 'd4', plotlimits, searchlimits,contours=[onesigma, twosigma], nbins=nbinsdist, normalbar=False)
    save_data([cont, X, Y, chisq_sub], ResultsDir + 'contourdata'+ fulltag + '_ct3_d4.pkl')
    

# c3, d4 (all others zero) 
cont, X, Y, chisq_sub = contour_pvalue_only(Process, 'pvalue'+ fulltag, '$gg\\rightarrow hhh$@' + str(Energy) + ' TeV, L=' + str(Luminosity) + ' fb$^{-1}$, $\\mathcal{P}(b \\rightarrow b ) =' + str(btagging) + ' $' + ', $\\alpha_\\mathrm{syst.} = ' + str(100*Systematics) +  '\%$', popt[Process], popt_eff, sigma_bkg, 'c3', 'd4', plotlimits, searchlimits,contours=[onesigma, twosigma], nbins=nbinsdist, normalbar=False)
save_data([cont, X, Y, chisq_sub], ResultsDir + 'contourdata'+ fulltag + '_c3_d4.pkl')


# XGBOOST:
#searchlimits[100][0] = [-3.0, 4.0]
#searchlimits[100][1] = [-1.0, 1.0]
#searchlimits[100][2] = [-0.5, 0.5]
#ssearchlimits[100][3] = [-10.0, 21.0]
# ct3, d4 (all others zero)
if MODEL != 'C3D4ONLY' and MODEL !='HEFT3':
    cont, X, Y, chisq_sub = contour_pvalue_only(Process, 'pvalueXGBOOST' + fulltag, '$gg\\rightarrow hhh$@' + str(Energy) + ' TeV, L=' + str(Luminosity) + ' fb$^{-1}$, $\\mathcal{P}(b \\rightarrow b ) =' + str(btagging) + ' $' + ', $\\alpha_\\mathrm{syst.} = ' + str(100*Systematics) +  '\%$', popt[Process], popt_eff_XGBOOST, sigma_bkg_xgboost, 'ct3', 'd4', plotlimits, searchlimits,contours=[onesigma, twosigma], nbins=nbinsdist, normalbar=False)
    save_data([cont, X, Y, chisq_sub], ResultsDir + 'contourdataXGBOOST'+ fulltag + '_ct3_d4.pkl')


# c3, d4 (all others zero) 
cont, X, Y, chisq_sub = contour_pvalue_only(Process, 'pvalueXGBOOST'+ fulltag, '$gg\\rightarrow hhh$@' + str(Energy) + ' TeV, L=' + str(Luminosity) + ' fb$^{-1}$, $\\mathcal{P}(b \\rightarrow b ) =' + str(btagging) + ' $' + ', $\\alpha_\\mathrm{syst.} = ' + str(100*Systematics) +  '\%$', popt[Process], popt_eff_XGBOOST, sigma_bkg_xgboost, 'c3', 'd4', plotlimits, searchlimits, contours=[onesigma, twosigma], nbins=nbinsdist, normalbar=False)
save_data([cont, X, Y, chisq_sub], ResultsDir + 'contourdataXGBOOST'+ fulltag + '_c3_d4.pkl')

########################
# MARGINALIZE OVER C3
########################
# c3, d4 (all others zero)
deltac3 = 0.05

plotlimits[100][0] = [-2.0, 2.0]
plotlimits[100][1] = [-1.0, 1.0]
plotlimits[100][2] = [-0.5, 0.5]
plotlimits[100][3] = [-60.0, 80.0]
searchlimits[100][0] = [-0.3, 0.3]
searchlimits[100][1] = [-1.0, 1.0]
searchlimits[100][2] = [-0.5, 0.5]
searchlimits[100][3] = [-8.0, 8.0]

# CUTS
cont, X, Y, chisq_sub = contour_pvalue_only(Process, 'pvalue_deltac3' + str(deltac3) + '_' + fulltag, '$gg\\rightarrow hhh$@' + str(Energy) + ' TeV, L=' + str(Luminosity) + ' fb$^{-1}$, $\\mathcal{P}(b \\rightarrow b ) =' + str(btagging) + ' $' + ', $\\alpha_\\mathrm{syst.} = ' + str(100*Systematics) +  '\%$', popt[Process], popt_eff, sigma_bkg, 'c3', 'd4', plotlimits, searchlimits,contours=[onesigma, twosigma], nbins=nbinsdist, normalbar=False, deltac3=deltac3)
save_data([cont, X, Y, chisq_sub], ResultsDir + 'contourdata_deltac3' + str(deltac3) + fulltag + '_c3_d4.pkl')

# XGBOOST 
cont, X, Y, chisq_sub = contour_pvalue_only(Process, 'pvalueXGBOOST_deltac3' + str(deltac3) + '_' + fulltag, '$gg\\rightarrow hhh$@' + str(Energy) + ' TeV, L=' + str(Luminosity) + ' fb$^{-1}$, $\\mathcal{P}(b \\rightarrow b ) =' + str(btagging) + ' $' + ', $\\alpha_\\mathrm{syst.} = ' + str(100*Systematics) +  '\%$', popt[Process], popt_eff_XGBOOST, sigma_bkg_xgboost, 'c3', 'd4', plotlimits, searchlimits, contours=[onesigma, twosigma], nbins=nbinsdist, normalbar=False, deltac3=deltac3)
save_data([cont, X, Y, chisq_sub], ResultsDir + 'contourdataXGBOOST_deltac3_' + str(deltac3) + fulltag + '_c3_d4.pkl')

####################################
# PRINT COEFFICIENTS OF XSEC FIT:
####################################

print('Fit coefficients for MODEL=', MODEL, '=',  popt[Process])


# ChargePairing
b-jet charge pairing in multi-Higgs boson production

## Repository Setup

1. First clone the repository:
```
git clone https://github.com/apapaefs/ChargePairing.git
```

2. In the directory of the repository pull the lfs files as well:
```
git lfs install
git lfs pull
git lfs checkout
```

## Analysis 
The main analyses are: 

```Code/HwSimPostAnalysis_smear_13.6_variables_CMS.cc``` (no charge pairing)
and 
```Code/HwSimPostAnalysis_smear_13.6_variables_CMS.cc``` (with charge pairing)

The current signal ($pp \rightarrow hhh \rightarrow 6b$) and background (QCD $pp\rightarrow 6b$) samples (at 13.6 TeV) are in:

```Signals/events/HW-gg_hhh_SM.root
Backgrounds/events/HW-all_events_6b_13.6_new1.root
```

each has 100k events. The corresponding MG5 total cross sections, cross sections after $h \to b\bar{b}$ (=0.5824) and b-tagging with $p=0.85$, and total number of events at 3/ab are:
MG5 XSEC [fb]	| XSEC * BR * BTAG * KFACS	| Nevents, 3000/fb
4.07E-02	| 6.07E-03	| 18.20
1.06E+03	| 7.98E+02	| 2.39E+06

#!/usr/bin/env bash
# for most figures: figure 2, 3d,e, 4, 6b-j, S5, S7, s8, s9, s10, s11, s12, s13, s16, s17 and data for: figure 5, S6, s14
python SVZ_molar_model_phytotypes.py 3 0 mesotrophic 0 0
python SVZ_molar_model_phytotypes.py 3 0 mesotrophic 1.4 0
python SVZ_molar_model_phytotypes.py 3 0 mesotrophic 0 1800
python SVZ_molar_model_phytotypes.py 3 0 mesotrophic 1.4 1800

python SIVZ_molar_model_phytotypes.py 3 0 mesotrophic 0 0 1 0 0
python SIVZ_molar_model_phytotypes.py 3 0 mesotrophic 1.4 0 1 0 0
python SIVZ_molar_model_phytotypes.py 3 0 mesotrophic 0 1800 1 0 0
python SIVZ_molar_model_phytotypes.py 3 0 mesotrophic 1.4 1800 1 0 0

python SVRZ_molar_model_phytotypes.py 3 mesotrophic 0.8 0 0 phir 5
python SVRZ_molar_model_phytotypes.py 3 mesotrophic 0.8 1.4 0 phir 5
python SVRZ_molar_model_phytotypes.py 3 mesotrophic 0.8 0 1800 phir 5
python SVRZ_molar_model_phytotypes.py 3 mesotrophic 0.8 1.4 1800 phir 5

python SIVRZ_molar_model_phytotypes.py 3 0 0.8 SIVRZ mesotrophic 0 0 lp_phir 5
python SIVRZ_molar_model_phytotypes.py 3 0 0.8 SIVRZ mesotrophic 1.4 0 lp_phir 5
python SIVRZ_molar_model_phytotypes.py 3 0 0.8 SIVRZ mesotrophic 0 1800 lp_phir 5
python SIVRZ_molar_model_phytotypes.py 3 0 0.8 SIVRZ mesotrophic 1.4 1800 lp_phir 5

python SVRZ_molar_model_phytotypes.py 3 mesotrophic 0.8 0 0 phir 0
python SVRZ_molar_model_phytotypes.py 3 mesotrophic 0.8 1.4 0 phir 0
python SVRZ_molar_model_phytotypes.py 3 mesotrophic 0.8 0 1800 phir 0
python SVRZ_molar_model_phytotypes.py 3 mesotrophic 0.8 1.4 1800 phir 0

python SIVRZ_molar_model_phytotypes.py 3 0 0.8 SIVRZ mesotrophic 0 0 lp_phir 0
python SIVRZ_molar_model_phytotypes.py 3 0 0.8 SIVRZ mesotrophic 1.4 0 lp_phir 0
python SIVRZ_molar_model_phytotypes.py 3 0 0.8 SIVRZ mesotrophic 0 1800 lp_phir 0
python SIVRZ_molar_model_phytotypes.py 3 0 0.8 SIVRZ mesotrophic 1.4 1800 lp_phir 0


# figure 2a-c
python SIVZ_MCT.py 3 0 mesotrophic 0 1 0 

# figure s8, figure 4 and figure s13 (time series)
python SVRZ_molar_model_phytotypes.py 3 mesotrophic 0.8 0 0 phir no
python SVRZ_molar_model_phytotypes.py 3 mesotrophic 0.8 1.4 0 phir no
python SVRZ_molar_model_phytotypes.py 3 mesotrophic 0.8 0 1800 phir no
python SVRZ_molar_model_phytotypes.py 3 mesotrophic 0.8 1.4 1800 phir no

python SIVRZ_molar_model_phytotypes.py 3 0 0.8 SIVRZ mesotrophic 0 0 phir 0
python SIVRZ_molar_model_phytotypes.py 3 0 0.8 SIVRZ mesotrophic 1.4 0 phir 0
python SIVRZ_molar_model_phytotypes.py 3 0 0.8 SIVRZ mesotrophic 0 1800 phir 0
python SIVRZ_molar_model_phytotypes.py 3 0 0.8 SIVRZ mesotrophic 1.4 1800 phir 0

python SVRZ_molar_model_phytotypes.py 3 mesotrophic 0.8 0 0 epsr no
python SVRZ_molar_model_phytotypes.py 3 mesotrophic 0.8 1.4 0 epsr no
python SVRZ_molar_model_phytotypes.py 3 mesotrophic 0.8 0 1800 epsr no
python SVRZ_molar_model_phytotypes.py 3 mesotrophic 0.8 1.4 1800 epsr no

python SIVRZ_molar_model_phytotypes.py 3 0 0.8 SIVRZ mesotrophic 0 0 epsr 0
python SIVRZ_molar_model_phytotypes.py 3 0 0.8 SIVRZ mesotrophic 1.4 0 epsr 0
python SIVRZ_molar_model_phytotypes.py 3 0 0.8 SIVRZ mesotrophic 0 1800 epsr 0
python SIVRZ_molar_model_phytotypes.py 3 0 0.8 SIVRZ mesotrophic 1.4 1800 epsr 0

# figure s9
python SIVZ_molar_model_phytotypes.python2 0 mesotrophic 1.4 0 1 0 0
python SIVZ_molar_model_phytotypes.python2 0 mesotrophic 1.4 1800 1 0 0
python SIVZ_molar_model_phytotypes.python1 0 mesotrophic 1.4 0 1 0 0
python SIVZ_molar_model_phytotypes.python1 0 mesotrophic 1.4 1800 1 0 0
python SIVZ_molar_model_phytotypes.python0 0 mesotrophic 1.4 0 1 0 0
python SIVZ_molar_model_phytotypes.python0 0 mesotrophic 1.4 1800 1 0 0

# figure s11
python SIVRZ_molar_model_phytotypes.python2 0 0.8 SIVRZ mesotrophic 1.4 0 lp_phir 5
python SIVRZ_molar_model_phytotypes.python2 0 0.8 SIVRZ mesotrophic 1.4 1800 lp_phir 5
python SIVRZ_molar_model_phytotypes.python1 0 0.8 SIVRZ mesotrophic 1.4 0 lp_phir 5
python SIVRZ_molar_model_phytotypes.python1 0 0.8 SIVRZ mesotrophic 1.4 1800 lp_phir 5
python SIVRZ_molar_model_phytotypes.python0 0 0.8 SIVRZ mesotrophic 1.4 0 lp_phir 5
python SIVRZ_molar_model_phytotypes.python0 0 0.8 SIVRZ mesotrophic 1.4 1800 lp_phir 5


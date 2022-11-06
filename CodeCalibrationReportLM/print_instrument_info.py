  #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 19 21:00:16 2022

@author: lindamegner


"""


from mats_l1_processing.experimental_utils import plot_CCDimage,read_all_files_in_protocol
from  mats_l1_processing.experimental_utils import readprotocol
from mats_l1_processing.instrument import Instrument





calibration_file='/Users/lindamegner/MATS/retrieval/git/MATS-L1-processing/scripts/calibration_data_linda.toml'



#CCDitems = read_CCDitems(directory)  # read in data

instrument = Instrument(calibration_file)

for CCDunit in enumerate([instrument.IR1, instrument.IR2, instrument.IR3, instrument.IR4, instrument.UV1, instrument.UV2, instrument.NADIR]):
    CCDunit=CCDunit[1] #Dont know why this is read in as tuple
    print(CCDunit.channel+" ro_avr_HSM "+str(CCDunit.ro_avr_HSM))
    print(CCDunit.channel+" ro_std_HSM "+str(CCDunit.ro_std_HSM))
    print(CCDunit.channel+" alpha_avr_HSM "+str(CCDunit.alpha_avr_HSM))
    print(CCDunit.channel+" alpha_std_HSM "+str(CCDunit.alpha_std_HSM))

    print(CCDunit.channel+" ro_avr_LSM "+str(CCDunit.ro_avr_LSM))
    print(CCDunit.channel+" ro_std_LSM "+str(CCDunit.ro_std_LSM))
    print(CCDunit.channel+" alpha_avr_LSM "+str(CCDunit.alpha_avr_LSM))
    print(CCDunit.channel+" alpha_std_LSM "+str(CCDunit.alpha_std_LSM))
        
    # 0D dark current subtraction stuff
    print(CCDunit.channel+" log_a_avr_HSM "+str(CCDunit.log_a_avr_HSM))
    print(CCDunit.channel+" log_a_std_HSM "+str(CCDunit.log_a_std_HSM))
    print(CCDunit.channel+" log_b_avr_HSM "+str(CCDunit.log_b_avr_HSM))
    print(CCDunit.channel+" log_b_std_HSM "+str(CCDunit.log_b_std_HSM))

    print(CCDunit.channel+" log_a_avr_LSM "+str(CCDunit.log_a_avr_LSM))
    print(CCDunit.channel+" log_a_std_LSM "+str(CCDunit.log_a_std_LSM))
    print(CCDunit.channel+" log_b_avr_LSM "+str(CCDunit.log_b_avr_LSM))
    print(CCDunit.channel+" log_b_std_LSM "+str(CCDunit.log_b_std_LSM))

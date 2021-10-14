import numpy as np
import math
from abc import ABC, abstractmethod
from random import shuffle, random

class Constants():

    #### 30 MINUTES STEP SIZE #####

    D = 16

    TIME_STEP_SIZE = 2 # minutes
    HOUR = 60 # minutes
    UNIT_T = TIME_STEP_SIZE/30

    ##################################################################
    ITER_TO_SWELLING = int(4 * (HOUR/TIME_STEP_SIZE) - 2) #6 # 8???
    PR_ASPERGILLUS_CHANGE = -math.log(0.5)/(6*(HOUR/TIME_STEP_SIZE))
    ITER_TO_GERMINATE = int(2 * (HOUR/TIME_STEP_SIZE) - 2) #2 # $4???
    ##################################################################
    ITER_TO_CHANGE_STATE = int(2 * (HOUR/TIME_STEP_SIZE)) #4
    ITER_TO_REST = 30*6
    ITER_TO_GROW = (int(HOUR/TIME_STEP_SIZE) - 1) #1
    ITER_TO_NEUTROPHIL_DIE = 12*int(HOUR/TIME_STEP_SIZE)


    VOXEL_VOL = 6.4e-11 # 0.26  # L
    HYPHAE_VOL = 1.056832e-12  # L
    CONIDIA_VOL = 4.844e-14 #L #5.324E-015 #4.844e-14
    MA_VOL = 4.849048e-12  # 4.849048e-15 #  --- WRONG!!! (4.849048e-12)
    SPACE_VOL = 6.4e-8  ## FOR NOW!!!!

    #ITER_TO_HEP_CHANGE = 10
    PR_BRANCH = 0.25


    TURNOVER_RATE = 1-math.log(1.2)/int(30/TIME_STEP_SIZE)#0.2 # 10.1124/jpet.118.250134 (approx) 0.2/h CHANGE!!!!
    #ERYTHROCYTE_TURNOVER_RATE = 1-math.log(1.2)/int(30/TIME_STEP_SIZE)

    LAC_QTTY = 4.3680464554587733e-17 * UNIT_T/2 # REF 62  HALF-HOUR????

    TAFC_UP = (1e-12/VOXEL_VOL) * UNIT_T/2 # (cell^-1.h^-1) Ref 218, 224 (1e-12/VOXEL_VOL) -- NOT UNIT_T

    TAFC_QTTY =  1e-15 * UNIT_T #2e-15 * UNIT_T # REF 15  (Ref 2.7e-13 219)   1.74e-15#1.74e-12 #TAFC secretion rate  3.0659e-15
    HEMOLYSIN_QTTY = 2.08e-10 * UNIT_T#5e-9 * UNIT_T # Arbitrary Unities
    MA_IL6_QTTY = 2.192281e-20 * UNIT_T / 3 # THIS IS CORRECT IT IS 3!!!!!!!
    MA_IL8_QTTY = 7.609164e-19 * UNIT_T / 3
    MA_MCP1_QTTY = 2.635868e-20 * UNIT_T / 3
    MA_MIP1B_QTTY = 2.684384e-20 * UNIT_T / 3
    MA_MIP2_QTTY = 1.659086e-19 * UNIT_T / 3
    MA_IL10_QTTY = 1.046032e-21 * UNIT_T / 3
    MA_TNF_QTTY = 4.826850e-20 * UNIT_T / 3
    MA_TGF_QTTY = 1.517909e-21 * UNIT_T / 3

    N_IL6_QTTY = 0.005875191*MA_IL6_QTTY
    N_IL8_QTTY = 0.005875191*MA_IL8_QTTY
    N_MIP2_QTTY = 0.005875191*MA_MIP2_QTTY
    N_TNF_QTTY = 0.005875191*MA_TNF_QTTY

    P_IL6_QTTY = 1*MA_IL6_QTTY
    P_IL8_QTTY = 1*MA_IL8_QTTY
    P_MCP1_QTTY = 1*MA_MCP1_QTTY
    P_MIP1B_QTTY = 1*MA_MIP1B_QTTY
    P_MIP2_QTTY = 1*MA_MIP2_QTTY
    P_TNF_QTTY = 1*MA_TNF_QTTY

    IL6_HALF_LIFE = 1+math.log(0.5)/(1 * HOUR/TIME_STEP_SIZE)#1
    IL8_HALF_LIFE = 1+math.log(0.5)/(1 * HOUR/TIME_STEP_SIZE)#4
    MCP1_HALF_LIFE = 1+math.log(0.5)/(1 * HOUR/TIME_STEP_SIZE)#2
    MIP1B_HALF_LIFE = 1+math.log(0.5)/(1 * HOUR/TIME_STEP_SIZE)#2
    MIP2_HALF_LIFE = 1+math.log(0.5)/(1 * HOUR/TIME_STEP_SIZE)#2
    IL10_HALF_LIFE = 1+math.log(0.5)/(1 * HOUR/TIME_STEP_SIZE)#3
    TNF_HALF_LIFE = 1+math.log(0.5)/(1 * HOUR/TIME_STEP_SIZE)#0.5
    TGF_HALF_LIFE = 1+math.log(0.5)/(1 * HOUR/TIME_STEP_SIZE)#0.1
    ANTI_TNF_HALF_LIFE = 1+math.log(0.5)/(5 * 24 * HOUR/TIME_STEP_SIZE)# five days -- 10.1002/eji.1830180221
    HEMOPEXIN_HALF_LIFE = 1+math.log(0.5)/(1.2 * 24 * HOUR/TIME_STEP_SIZE) # 1.2 days -- 10.1006/abbi.1993.1014

    Kd_IL6 = 3.3e-10
    Kd_IL8 = 1.045e-9
    Kd_MCP1 = 5.686549e-10
    Kd_MIP1B = 1.8e-10
    Kd_MIP2 = 9.166739837e-11
    Kd_IL10 = 1.4e-10
    Kd_TNF = 3.26e-10
    Kd_TGF = 2.65e-11
    Kd_KC = 1.0e-9 #doi/10.1074/jbc.M112.443762
    Kd_Hep = 8.55e-07 # REF 223
    Kd_LIP = 7.90456448805514E-05 #2.762975e-05
    Kd_HEMO = 19.987634208144584#185.18518518518516/2 # Arbitrary Unitis x h^-1 10.1016/j.ijmm.2011.04.016
    Kd_MA_IRON = 0.0020799959084752307

    MA_IRON_IMPORT_RATE = 5.3333e-12/VOXEL_VOL # L*cell^-1*h^-1. REF 27
    MA_IRON_EXPORT_RATE = 1367.3051298168639/VOXEL_VOL #calculations made based on an MA iron concentration of 2.08 mM and 6.5 uM of apo-Tf (see Pedro's paper)
    MA_HEME_IMPORT_RATE = 0.0833328125

    #TIME_STEP = 30
    CYT_BIND_T = 35 # SEE MANUAL (REF 199)
    IRON_AFNT_T = 60
    IRON_IMP_EXP_T = 60
    PHAG_AFNT_T = 1
    HYPHAE_INT_T = 60
    REL_CYT_BIND_UNIT_T = TIME_STEP_SIZE/CYT_BIND_T
    REL_IRON_AFNT_UNIT_T = TIME_STEP_SIZE/IRON_AFNT_T
    REL_PHAG_AFNT_UNIT_T = TIME_STEP_SIZE/PHAG_AFNT_T
    REL_IRON_IMP_EXP_UNIT_T = TIME_STEP_SIZE/IRON_IMP_EXP_T
    REL_HYPHAEL_INT_UNIT_T = TIME_STEP_SIZE/PHAG_AFNT_T
    REL_N_HYPHAE_INT_UNIT_T = TIME_STEP_SIZE/HOUR
    STD_UNIT_T = TIME_STEP_SIZE/HOUR

    DRIFT_BIAS = 1e-100 #0.9 # REF 156 (Est.)
    PR_MOVE_ASP = 0.75 # DUMMY
    MA_MOVE_RATE_REST = 1.44*TIME_STEP_SIZE/40 # 1.44--doi:10.1371/journal.pone.0004693 -NOT
    MA_MOVE_RATE_ACT = 1.44*TIME_STEP_SIZE/40
    REC_BIAS = 0.9995 # DUMMY VALUE CREATED TO AVOID INFINTY LOOP!

    MAX_N = 522.5#15 # for an alveoli sac with 6 alveoli (Review)
    MAX_MA = 209
    MIN_MA = 15
    #MA_MAX = 360 # for an alveoli sac with 6 alveoli (Review) using MIP-1b

    #N_REC = 1.0884e-11
    #MA_REC = 5*N_REC # the bigger the "worst"
    NUM_ALVEOLI = 3.26e+06 #REF 242, 184

    #IRON_EXP_RATE = 7.035e4 # mol*cell^-1*h^-1 (It is based on mol not molar). Total Iron inside a macrophage: 8.713E-21 mols.
    #IRON_IMP_RATE = 0.322   # mol*cell^-1*h^-1 (It is based on mol not molar). Total Iron in a voxel (40x40x40): 1.904E-15 mols.

    PR_MA_PHAG =  1 - math.exp(-(1/VOXEL_VOL)*REL_PHAG_AFNT_UNIT_T/1.32489230813214E+10) # 30 min --> 1 - exp(-cells*t/Kd) --> Kd = 1.32489230813214E+10
    PR_N_PHAG =   1 - math.exp(-(1/VOXEL_VOL)*REL_PHAG_AFNT_UNIT_T/1.96177129709014E+11) # 30 min --> Kd = 1.96177129709014E+11
    PR_E_PHAG =   1 - math.exp(-(1/VOXEL_VOL)*REL_PHAG_AFNT_UNIT_T/6.80269538729756E+11) # 30 min --> Kd = 6.80269538729756E+11
    PR_N_HYPHAE = 1 - math.exp(-(1/VOXEL_VOL)*REL_N_HYPHAE_INT_UNIT_T/2.02201143330207E+09) # 0.5 h --> Kd = 2.02201143330207E+09
    PR_MA_HYPHAE = 1 - math.exp(-(1/VOXEL_VOL)*REL_N_HYPHAE_INT_UNIT_T/5.02201143330207E+9)  #Kd ~10x Neut. (REF 71)
    #PR_MA_HYPHAE = 1 - math.exp(-(1 / VOXEL_VOL) * REL_N_HYPHAE_INT_UNIT_T / 3.3333e+10)  # 0.5 h --> Kd = 2.02201143330207E+09
    PR_MA_PHAG_ERYT = 1 - math.exp(-(1/VOXEL_VOL)*REL_N_HYPHAE_INT_UNIT_T/73873893586.4061) #https://doi.org/10.1002/(SICI)1097-0320(19971015)30:5<269::AID-CYTO8>3.0.CO;2-C
    PR_P_INT = PR_E_PHAG #DUMMY VALUE???


    MA_MAX_CONIDIA = 18
    N_MAX_CONIDIA = 3
    E_MAX_CONIDIA = 18

    PR_KILL = -math.log(0.01)/(12*HOUR/TIME_STEP_SIZE) #10.1080/13693780400029247

    K_M_TF_TAFC = 2.514985e-3
    K_M_TF_LAC = 2.5052031141601793e-3 ### REF 62 (2.5-10% uptake in 1-4h)


    MA_INTERNAL_IRON = 1.0086e-14 # 2.08e-3*4.849048e-12 (CORRECT VALUE)
    CONIDIA_INIT_IRON = Kd_LIP * CONIDIA_VOL

    NEUTROPHIL_HALF_LIFE = - math.log(0.5) / (6 * (HOUR/TIME_STEP_SIZE)) # 12h (half-life infection) doi.org/10.1189/jlb.1112571
    MA_HALF_LIFE = -math.log(0.5) / ( 1 * 24 * (HOUR/TIME_STEP_SIZE))

    #MA_OCCUPANCY = 6 # ~6 times bigger then neutrophils

    #SECRETING_HALF_LIFE = 0.9

    #PR_INT_RESTING = 0.15
    #PR_INT_ACTIVE = 0.3
    #MAX_INTERNALIZED_CONIDIA = 50

    #TF_ENHANCE = 2

    #PHAGOSOME_UP_RATE = 0.025

    P1 =  0.2734
    P2 = -1.1292
    P3 =  0.8552

    IL6_THRESHOLD = 1.372243e-10#1e-11 # 5e-14 #Molar (threshold to activate liver)
    HEP_INTERCEPT = -0.3141 # [IL6] > 1e-11
    HEP_SLOPE = 0.7793 # [IL6] > 1e-11

    TF_INTERCEPT = -1.194e-05
    TF_SLOPE = -5.523e-06
    THRESHOLD_LOG_HEP = -8
    THRESHOLD_HEP = math.pow(10, THRESHOLD_LOG_HEP)

    DEFAULT_APOTF_REL_CONCENTRATION = 0.4
    DEFAULT_TFFE_REL_CONCENTRATION = 0.1657
    DEFAULT_TFFE2_REL_CONCENTRATION = 0.4343

    DEFAULT_TF_CONCENTRATION = (TF_INTERCEPT + TF_SLOPE * THRESHOLD_LOG_HEP) * VOXEL_VOL #1.3490016161826808e-16
    DEFAULT_APOTF_CONCENTRATION = DEFAULT_APOTF_REL_CONCENTRATION * DEFAULT_TF_CONCENTRATION
    DEFAULT_TFFE_CONCENTRATION = DEFAULT_TFFE_REL_CONCENTRATION * DEFAULT_TF_CONCENTRATION
    DEFAULT_TFFE2_CONCENTRATION = DEFAULT_TFFE2_REL_CONCENTRATION * DEFAULT_TF_CONCENTRATION

    #MA_DEFAULT_IRON_CONCENTRATION = 2.08e-3

    RECRUITMENT_RATE = 2.0#2.9333567478400004#0.5e18#3.0e18 * STD_UNIT_T #1.5e+19  # 0.04 #0.04 ATTENTION HALVED!

    N_REC_MUL = 1.0
    MA_REC_MUL = 1.0

    N_FRAC = 1.0
    #N_REC_RATE = 1.2054830189634535e+19  # 0.2


    ANTI_TNFA_REACT_TIME_UNIT = 120 #sec
    K_M_ANTI_TNFA = 6.9737e-07  # (M^-1.sec^-1)  http://www.jimmunol.org/content/162/10/6040.full#ref-list-1
    ANTI_TNFA_SYSTEM_CONCENTRATION_REF = 2e-8*VOXEL_VOL
    ANTI_TNFA_SYSTEM_CONCENTRATION = ANTI_TNFA_SYSTEM_CONCENTRATION_REF # https://doi.org/10.1016/0022-1759(95)00278-2


    HEMOGLOBIN_UPTAKE_RATE = (1e-13/VOXEL_VOL) * UNIT_T # TAFC_UP
    ERYTROCYTE_HEMOGLOBIN_CONCENTRATION = 4.6875e-16 #mol
    MAX_ERYTHROCYTE_VOXEL = 180


    ESTB_KM = 4e-4 #10.1128/EC.00066-07
    ESTB_HALF_LIFE = 1+math.log(0.5)/(24 * HOUR/TIME_STEP_SIZE)#2.9 - Advances in Enzymology and related Areas of Molecular Biology, Alton Meister, v39, 1973, p227
    ESTB_SYSTEM_CONCENTRATION = 0
    ESTB_KCAT = 3096 # 10.1128/AEM.65.8.3470-3472.1999 (GENERIC FUNGI ESTERASE)

    HEMOPEXIN_SYSTEM_CONCENTRATION = 0 # NOT CORRECT
    HEMOPEXIN_SYSTEM_CONCENTRATION_REF = 3.1e-7 * VOXEL_VOL # (Luis SV and 10.3181/00379727-75-18083)
    HEMOPEXIN_KM = 1e-9 # Kd < 1pM -- 10.1006/abbi.1993.1014 (Kd is not Km!!!)
    HEMOPEXIN_KCAT = 100

    #HEMOLYSIS_RATE = 0.0005
    HEMOLYSIN = False

    ITERATION_TO_KILL = 15

class Id():
    ID = 0

    @staticmethod
    def get_id():
        Id.ID = Id.ID + 1
        return Id.ID


class Util():

    cosTheta = math.cos(math.pi/4.0)
    sinTheta = math.sin(math.pi / 4.0);

    @staticmethod
    def michaelianKinetics(substract1, substract2, km, h, Kcat = 1, v=Constants.VOXEL_VOL): ### CAUTION CHANGED!!!!
        substract1 = substract1 / v # transform into M
        substract2 = substract2 / v
        ######### CHANGE START HERE!!!!!
        substract = substract1
        enzime = substract2
        """
        if substract1 > substract2:
            enzime = substract2
            substract = substract1
        else:
            enzime = substract1;
            substract = substract2;
        """
        ######### END HERE!!!!

        return h * Kcat * enzime * substract * v  / (substract + km) # (*v) transform back into mol

    @staticmethod
    def activation_function(x, kd, h, v=Constants.VOXEL_VOL, b=1):
        x = x/v #CONVERT MOL TO MOLAR
        return h * (1 - b*math.exp(-(x/kd)))

    @staticmethod
    def turnover_rate(x, x_system=0, k=Constants.TURNOVER_RATE, t=Constants.REL_CYT_BIND_UNIT_T, v=Constants.VOXEL_VOL):
        if x == 0 and x_system == 0:
            return 0
        x = x / v  # CONVERT MOL TO MOLAR
        x_system = x_system / v
        y = (x - x_system)*math.exp(-k*t) + x_system
        return y / x

    @staticmethod
    def iron_tf_reaction(iron, Tf, TfFe):
        total_binding_site = 2*(Tf + TfFe) # That is wright 2*(Tf + TfFe)!
        total_iron = iron + TfFe # it does not count TfFe2
        if total_iron <= 0 or total_binding_site <= 0:
            return 0.0
        rel_total_iron = total_iron/total_binding_site
        rel_total_iron = rel_total_iron if rel_total_iron <= 1.0 else 1.0
        rel_TfFe = Constants.P1*rel_total_iron*rel_total_iron*rel_total_iron + Constants.P2*rel_total_iron*rel_total_iron + Constants.P3*rel_total_iron
        rel_TfFe = rel_TfFe if rel_TfFe > 0.0 else 0.0
        return rel_TfFe

    @staticmethod
    def rotatation_matrix(phi):
        return np.array(
            [
                [Util.cosTheta, 0, Util.sinTheta],\
                [math.sin(phi)*Util.sinTheta, math.cos(phi), -math.sin(phi)*Util.cosTheta],\
                [-math.cos(phi)*Util.sinTheta, 0, math.cos(phi)*Util.cosTheta]
            ]
        )

    @staticmethod
    def normalize(v):
        return v / np.linalg.norm(v)

    @staticmethod
    def gram_schmidt(V):
        n = len(V)
        for i in range(n):
            V[:,i] = Util.normalize(V[:,i])
            for j in range(i+1,n):
                V[:,j] = V[:,j] - ((V[:,j].dot(V[:,i]))/(V[:,i].dot(V[:,i])))*V[:,i]

        return(V)

    @staticmethod
    def random_sphere_point() -> np.ndarray:
        u: np.ndarray = 2 * np.random.random_sample(size=2) - 1
        while np.linalg.norm(u) > 1.0:
            u = 2 * np.random.random_sample(size=2) - 1

        norm_squared_u = float(np.dot(u, u))
        return np.array(
            [
                2 * u[0] * np.sqrt(1 - norm_squared_u),
                2 * u[1] * np.sqrt(1 - norm_squared_u),
                1 - 2 * norm_squared_u,
                ],
            dtype=np.float64,
        )


class SimRand(ABC):
    pass

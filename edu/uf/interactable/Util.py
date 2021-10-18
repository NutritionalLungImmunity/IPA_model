import numpy as np
import math
from abc import ABC, abstractmethod
from random import shuffle, random

class Constants():

    D = 16 # diffusion constant. Lead to 850 um^2/min

    TIME_STEP_SIZE = 2 # minutes
    HALF_HOUR = 30 # minutes
    CYT_BIND_T = 35 # 35 came from a manual with curves dose response; Their experiments is done in 35 minutes
    HOUR = 60 #minutes
    UNIT_T = TIME_STEP_SIZE/HALF_HOUR
    STD_UNIT_T = TIME_STEP_SIZE/HOUR
    REL_CYT_BIND_UNIT_T = TIME_STEP_SIZE/CYT_BIND_T
    REL_IRON_IMP_EXP_UNIT_T = STD_UNIT_T

    ITER_TO_SWELLING = 118
    PR_ASPERGILLUS_CHANGE = 0.003850817669777474 #probability to germinate in the next time-step.
    ITER_TO_GERMINATE = 58 # Iterations until germination start

    ITER_TO_CHANGE_STATE = 60
    ITER_TO_REST = 30*6
    ITER_TO_GROW = 29 # 40 um every 58 min

    VOXEL_VOL = 6.4e-11 # L
    VOXEL_LEN = 40 #um
    HYPHAE_VOL = 1.056832e-12  # L
    CONIDIA_VOL = 4.844e-14 # L
    MA_VOL = 4.849048e-12 # L
    SPACE_VOL = 6.4e-8 # L

    PR_BRANCH = 0.25 # 25%


    TURNOVER_RATE = 0.9878452295470697
    LAC_QTTY = 4.3680e-17 * STD_UNIT_T

    TAFC_UP = (1e-12/VOXEL_VOL) * STD_UNIT_T

    TAFC_QTTY = 1.0000e-15 * STD_UNIT_T
    MA_IL6_QTTY = 1.4615e-20 * STD_UNIT_T
    MA_MIP1B_QTTY = 1.7896e-20 * STD_UNIT_T
    MA_MIP2_QTTY = 1.1061e-19 * STD_UNIT_T
    MA_IL10_QTTY = 6.9735e-22 * STD_UNIT_T
    MA_TNF_QTTY = 3.2179e-20 * STD_UNIT_T
    MA_TGF_QTTY = 1.0119e-21 * STD_UNIT_T

    N_IL6_QTTY = 0.005875191*MA_IL6_QTTY
    N_MIP2_QTTY = 0.005875191*MA_MIP2_QTTY
    N_TNF_QTTY = 0.005875191*MA_TNF_QTTY

    P_IL6_QTTY = MA_IL6_QTTY
    P_MIP1B_QTTY = MA_MIP1B_QTTY
    P_MIP2_QTTY = MA_MIP2_QTTY
    P_TNF_QTTY = MA_TNF_QTTY

    #degradation rate: [Cytokine] = HALF_LIFE * [Cytokine]; HALF_LIFE=0.9768950939813351 -> 1h half-life
    IL6_HALF_LIFE = 0.9768950939813351
    MIP1B_HALF_LIFE = 0.9768950939813351
    MIP2_HALF_LIFE = 0.9768950939813351
    IL10_HALF_LIFE = 0.9768950939813351
    IL8_HALF_LIFE = 0.9768950939813351
    MCP1_HALF_LIFE = 0.9768950939813351
    TNF_HALF_LIFE = 0.9768950939813351
    TGF_HALF_LIFE = 0.9768950939813351
    ANTI_TNF_HALF_LIFE = 0.9998074591165111# five days -- 10.1002/eji.1830180221

    Kd_IL6 = 3.3e-10
    Kd_MIP1B = 1.8e-10
    Kd_MIP2 = 9.166739837e-11
    Kd_IL10 = 1.4e-10
    Kd_TNF = 3.26e-10
    Kd_TGF = 2.65e-11
    Kd_Hep = 8.55e-07
    Kd_LIP = 7.90456448805514E-05
    Kd_MA_IRON = 0.0020799959084752307

    MA_IRON_IMPORT_RATE = 5.3333e-12/VOXEL_VOL
    MA_IRON_EXPORT_RATE = 1367.3051298168639/VOXEL_VOL

    DRIFT_BIAS = 1e-100
    PR_MOVE_ASP = 0.75 # DUMMY

    #Average number of displaced voxels
    MA_MOVE_RATE_REST = 1.44*TIME_STEP_SIZE/VOXEL_LEN
    MA_MOVE_RATE_ACT = 1.44*TIME_STEP_SIZE/VOXEL_LEN


    REC_BIAS = 0.9995 # DUMMY VALUE CREATED TO AVOID INFINTY LOOP!

    MAX_N = 522
    MAX_MA = 209
    MIN_MA = 15

    PR_MA_PHAG = 0.9054552746554831
    PR_N_PHAG = 0.1472550818938948
    PR_N_HYPHAE = 0.2270815966867843
    PR_MA_HYPHAE = 0.09851340404770559
    PR_P_INT = 0.04489851158453906

    MA_MAX_CONIDIA = 18
    N_MAX_CONIDIA = 3

    PR_KILL = 0.012792139405522474

    K_M_TF_TAFC = 2.514985e-3
    K_M_TF_LAC = 2.5052031141601793e-3


    MA_INTERNAL_IRON = 1.0086e-14 # mols
    CONIDIA_INIT_IRON = Kd_LIP * CONIDIA_VOL # mols

    #probability leukocytes die in one time-step (2 min): 6 and 24 hours half-life respectively
    NEUTROPHIL_HALF_LIFE = 0.003850817669777474
    MA_HALF_LIFE = 0.0009627044174443685

    #These are phenomenological parameters of a polynome used to approximate an ODE.
    #That ODE computes the relative amounts of TfFe and TfFe2 after the reaction of Transferrin (Tf) and Iron (Fe).
    #That ODE is dimension-free and assume that both binding sites of Tf has the same afinity and that there is no cooperativity.
    #That is a simplification.
    P1 =  0.2734
    P2 = -1.1292
    P3 =  0.8552

    IL6_THRESHOLD = 1.372243e-10 # mols
    HEP_INTERCEPT = -0.3141
    HEP_SLOPE = 0.7793

    TF_INTERCEPT = -1.194e-05
    TF_SLOPE = -5.523e-06
    THRESHOLD_LOG_HEP = -8
    THRESHOLD_HEP = math.pow(10, THRESHOLD_LOG_HEP)

    DEFAULT_APOTF_REL_CONCENTRATION = 0.4 #40%
    DEFAULT_TFFE_REL_CONCENTRATION = 0.1657 #16.57%
    DEFAULT_TFFE2_REL_CONCENTRATION = 0.4343 #43.43%

    DEFAULT_TF_CONCENTRATION = (TF_INTERCEPT + TF_SLOPE * THRESHOLD_LOG_HEP) * VOXEL_VOL
    DEFAULT_APOTF_CONCENTRATION = DEFAULT_APOTF_REL_CONCENTRATION * DEFAULT_TF_CONCENTRATION
    DEFAULT_TFFE_CONCENTRATION = DEFAULT_TFFE_REL_CONCENTRATION * DEFAULT_TF_CONCENTRATION
    DEFAULT_TFFE2_CONCENTRATION = DEFAULT_TFFE2_REL_CONCENTRATION * DEFAULT_TF_CONCENTRATION

    RECRUITMENT_RATE = 2.0 #arbitrary unity

    KCAT_ANTI_TNF = 120
    K_M_ANTI_TNFA = 6.9737e-07
    ANTI_TNFA_SYSTEM_CONCENTRATION_REF = 2e-8*VOXEL_VOL
    ANTI_TNFA_SYSTEM_CONCENTRATION = ANTI_TNFA_SYSTEM_CONCENTRATION_REF # https://doi.org/10.1016/0022-1759(95)00278-2

    ITERATION_TO_KILL = 15

    MA_REC_MUL = 1

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

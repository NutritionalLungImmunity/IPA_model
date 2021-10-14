import math
from edu.uf.interactable.MIP1B import *
import numpy as np

class Liver():

    name = "Liver"

    def __init__(self):
        self.id = Id.get_id()
        self._log_hepcidin = None

        self.vinteract = np.vectorize(self.interact)

    def interact(self, interactable):
        itype = type(interactable)
        if itype is Liver:
            return False
        elif itype is Macrophage:
            return False
        elif itype is Neutrophil:
            return False
        elif itype is Pneumocytes:
            return False
        elif itype is Afumigatus:
            return False
        elif itype is Iron:
            return False
        elif itype is Transferrin:
            if self._log_hepcidin is None or self._log_hepcidin < Constants.THRESHOLD_LOG_HEP:
                tf = Constants.TF_INTERCEPT + Constants.TF_SLOPE * Constants.THRESHOLD_LOG_HEP
            else:
                tf = Constants.TF_INTERCEPT + Constants.TF_SLOPE * self._log_hepcidin
            tf = tf#*0.25142602860942986

            rate_tf    = Util.turnover_rate(interactable.get("Tf"), tf * Constants.DEFAULT_APOTF_REL_CONCENTRATION * Constants.VOXEL_VOL) - 1
            rate_tffe  = Util.turnover_rate(interactable.get("TfFe"), tf * Constants.DEFAULT_TFFE_REL_CONCENTRATION * Constants.VOXEL_VOL) - 1
            rate_tffe2 = Util.turnover_rate(interactable.get("TfFe2"), tf * Constants.DEFAULT_TFFE2_REL_CONCENTRATION * Constants.VOXEL_VOL) - 1

            interactable.pinc(rate_tf, "Tf")
            interactable.pinc(rate_tffe, "TfFe")
            interactable.pinc(rate_tffe2, "TfFe2")

            return True
        elif itype is TAFC:
            return False
        elif itype is Lactoferrin:
            return False
        elif itype is IL8:
            return False
        elif itype is IL10:
            return False
        elif itype is TNFa:
            return False
        elif itype is TGFb:
            return False
        elif itype is MCP1:
            return False
        elif itype is MIP2:
            return False
        elif itype is MIP1B:
            return False
        elif itype is ROS:
            return False
        elif itype is IL6:
            global_il6_concentration = interactable.total_molecule[0]/(2*Constants.SPACE_VOL) #div 2 : serum
            if global_il6_concentration > Constants.IL6_THRESHOLD:
                self._log_hepcidin = Constants.HEP_INTERCEPT + Constants.HEP_SLOPE*math.log(global_il6_concentration, 10)
            else:
                self._log_hepcidin = None
            return True
        elif itype is Hepcidin:
            if self._log_hepcidin == None or self._log_hepcidin > Constants.THRESHOLD_LOG_HEP:
                rate = Util.turnover_rate(interactable.get("Hepcidin"), Constants.THRESHOLD_HEP * Constants.VOXEL_VOL) - 1
            else:
                rate = Util.turnover_rate(interactable.get("Hepcidin"), math.pow(10, self._log_hepcidin) * Constants.VOXEL_VOL) - 1
            interactable.pinc(rate)
            return True
        return interactable.interact(self)
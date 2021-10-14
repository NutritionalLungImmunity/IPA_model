from edu.uf.interactable.MCP1 import *

class MIP2(Chemokine):
    name = "MIP2"

    total_molecule = [0]
    NUM_STATES = 1
    threshold = -1

    def _init(self):
        #from edu.uchc.interactable.Cells import Neutrophil
        Neutrophil.chemokine = MIP2

    def degrade(self, p=Constants.MIP2_HALF_LIFE):
        #pass
        super().degrade(p)
        super().degrade(Util.turnover_rate(x=1, x_system=0, k=Constants.TURNOVER_RATE, t=Constants.REL_CYT_BIND_UNIT_T, v=1))

    def compute_total_molecule(self):
        MIP2.total_molecule[0] = MIP2.total_molecule[0] + self.values[0]

    def get_index(self, str):
        return 0

    def inc_total_molecule(self, index, inc):
        MIP2.total_molecule[index] = MIP2.total_molecule[index] + inc

    def interact(self, interactable):
        #from edu.uchc.interactable.Cells import Macrophage, Afumigatus, Neutrophil, Pneumocytes, Phagocyte
        itype = type(interactable)
        if itype is MIP2:
            return False
        if itype is MCP1:
            return False
        if itype is TGFb:
            return False
        if itype is Hepcidin:
            return False
        if itype is IL10:
            return False
        if itype is Iron:
            return False
        if itype is Transferrin:
            return False
        if itype is TAFC:
            return False
        if itype is Lactoferrin:
            return False
        if itype is ROS:
            return False
        if itype is Afumigatus:
            return False
        if itype is IL6:
            return False
        if itype is TNFa:
            return False
        if itype is IL8:
            return False
        if itype is Neutrophil:
            if interactable.status == Phagocyte.RESTING:
                if Util.activation_function(self.get(0), Constants.Kd_MIP2, Constants.STD_UNIT_T) > random():
                    interactable.status = Phagocyte.ACTIVATING
            elif interactable.tnfa:#interactable.status == Phagocyte.ACTIVE and interactable.state == Neutrophil.INTERACTING:
                self.inc(Constants.N_MIP2_QTTY, 0)
                if Util.activation_function(self.get(0), Constants.Kd_MIP2, Constants.STD_UNIT_T) > random():
                    interactable.interaction = 0
            #if Util.activation_function(self.values[0], Constants.Kd_MIP2) > random():
            #    self.pdec(0.5)
            return True
        if itype is Pneumocytes:
            if interactable.tnfa:#interactable.status == Phagocyte.ACTIVE:
                self.inc(Constants.P_MIP2_QTTY, 0)
            return True
        #if type(interactable) is Hepatocytes:
        #    return False
        if itype is Macrophage:
            if interactable.tnfa:#interactable.status == Phagocyte.ACTIVE:# and interactable.state == Neutrophil.INTERACTING:
                self.inc(Constants.MA_MIP2_QTTY, 0)
            return True
        return interactable.interact(self)

    def chemoatract(self, drift_bias=Constants.DRIFT_BIAS):
        return Util.activation_function(self.values[0], Constants.Kd_MIP2, Constants.STD_UNIT_T) + drift_bias
        #return drift_bias*Constants.Kd_MIP1B + self.get(0)
from edu.uf.interactable.MIP2 import *

class MIP1B(Chemokine):
    name = "MIP1B"

    total_molecule = [0]
    NUM_STATES = 1
    threshold = -1

    def _init(self):
        #from edu.uchc.interactable.Cells import Macrophage
        Macrophage.chemokine = MIP1B

    def degrade(self, p=Constants.MIP1B_HALF_LIFE):
        super().degrade(p)
        super().degrade(Util.turnover_rate(x=1, x_system=0, k=Constants.TURNOVER_RATE, t=Constants.REL_CYT_BIND_UNIT_T, v=1))

    def get_index(self, str):
        return 0

    def inc_total_molecule(self, index, inc):
        MIP1B.total_molecule[index] = MIP1B.total_molecule[index] + inc

    def compute_total_molecule(self):
        MIP1B.total_molecule[0] = MIP1B.total_molecule[0] + self.values[0]

    def interact(self, interactable):
        #from edu.uchc.interactable.Cells import Macrophage, Afumigatus, Neutrophil, Pneumocytes, Phagocyte
        itype = type(interactable)
        if itype is MIP1B:
            return False
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
            return False
        if itype is Pneumocytes:
            if interactable.tnfa:#interactable.status == Phagocyte.ACTIVE:
                self.inc(Constants.P_MIP1B_QTTY, 0)
            return True
        #if type(interactable) is Hepatocytes:
        #    return False
        if itype is Macrophage:
            if interactable.tnfa:#interactable.status == Phagocyte.ACTIVE:# and interactable.state == Neutrophil.INTERACTING:
                self.inc(Constants.MA_MIP1B_QTTY, 0)
            #if Util.activation_function(self.values[0], Constants.Kd_MIP1B) > random():
            #    self.pdec(0.5)
            return True
        return interactable.interact(self)

    def chemoatract(self, drift_bias=Constants.DRIFT_BIAS):
        return Util.activation_function(self.values[0], Constants.Kd_MIP1B, Constants.STD_UNIT_T) + drift_bias
        #return drift_bias*Constants.Kd_MIP1B + self.get(0)
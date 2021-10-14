from edu.uf.interactable.TGFb import *

class MCP1(Chemokine):
    name = "MCP1"

    total_molecule = [0]
    NUM_STATES = 1
    threshold = -1

    def _init(self):
        pass

    def degrade(self, p=Constants.MCP1_HALF_LIFE):
        super().degrade(p)
        super().degrade(Util.turnover_rate(x=1, x_system=0, k=Constants.TURNOVER_RATE, t=Constants.REL_CYT_BIND_UNIT_T, v=1))

    def get_index(self, str):
        return 0

    def inc_total_molecule(self, index, inc):
        MCP1.total_molecule[index] = MCP1.total_molecule[index] + inc

    def compute_total_molecule(self):
        MCP1.total_molecule[0] = MCP1.total_molecule[0] + self.values[0]

    def interact(self, interactable):
        #from edu.uchc.interactable.Cells import Macrophage, Afumigatus, Neutrophil, Pneumocytes, Phagocyte
        itype = type(interactable)
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
                self.inc(Constants.P_MCP1_QTTY, 0)
            return True
        #if type(interactable) is Hepatocytes:
        #    return False
        if itype is Macrophage:
            if interactable.tnfa:#interactable.status == Phagocyte.ACTIVE and interactable.state == Neutrophil.INTERACTING:
                self.inc(Constants.MA_MCP1_QTTY, 0)
            return True
        return interactable.interact(self)

    def chemoatract(self, drift_bias=Constants.DRIFT_BIAS):
        return Util.activation_function(self.values[0], Constants.Kd_MCP1, Constants.STD_UNIT_T) + drift_bias
        #return drift_bias*Constants.Kd_MIP1B + self.get(0)
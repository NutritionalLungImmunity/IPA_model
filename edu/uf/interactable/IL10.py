from edu.uf.interactable.IL8 import *

class IL10(Molecule):
    name = "IL10"

    total_molecule = [0]
    NUM_STATES = 1
    threshold = -1

    def degrade(self, p=Constants.IL10_HALF_LIFE):
        super().degrade(p)
        super().degrade(Util.turnover_rate(x=1, x_system=0, k=Constants.TURNOVER_RATE, t=Constants.REL_CYT_BIND_UNIT_T, v=1))

    def get_index(self, str):
        return 0

    def inc_total_molecule(self, index, inc):
        IL10.total_molecule[index] = IL10.total_molecule[index] + inc

    def compute_total_molecule(self):
        IL10.total_molecule[0] = IL10.total_molecule[0] + self.values[0]

    def interact(self, interactable):
        #from edu.uchc.interactable.Cells import Macrophage, Afumigatus, Neutrophil, Phagocyte
        itype = type(interactable)
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
        if itype is Macrophage:# or type(interactable) is Neutrophil:
            if (interactable.status == Phagocyte.ACTIVE and interactable.state == Neutrophil.INTERACTING):
                self.inc(Constants.MA_IL10_QTTY, 0)
            if interactable.status == Macrophage.ANERGIC:
                pass
                #self.inc(10*Constants.MA_IL10_QTTY, 0)
            if interactable.status != Phagocyte.DEAD and interactable.status != Phagocyte.APOPTOTIC and interactable.status != Phagocyte.NECROTIC:
                if Util.activation_function(self.get(0), Constants.Kd_IL10, Constants.STD_UNIT_T) > random():
                    interactable.status = Phagocyte.INACTIVATING if interactable.status != Phagocyte.INACTIVE else Phagocyte.INACTIVE
                    interactable.interation = 0
            return True
        if itype is Neutrophil:
            return False
        return interactable.interact(self)
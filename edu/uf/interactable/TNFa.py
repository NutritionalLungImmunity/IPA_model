from edu.uf.interactable.IL6 import *

class TNFa(Molecule):

    name = "TNFa"

    total_molecule = [0]
    NUM_STATES = 1
    threshold = 0

    def degrade(self, p=Constants.TNF_HALF_LIFE):
        super().degrade(p)
        super().degrade(Util.turnover_rate(x=1, x_system=0, k=Constants.TURNOVER_RATE, t=Constants.REL_CYT_BIND_UNIT_T, v=1))

    def get_index(self, str):
        return 0

    def inc_total_molecule(self, index, inc):
        TNFa.total_molecule[index] = TNFa.total_molecule[index] + inc

    def compute_total_molecule(self):
        if self.values[0] < 0:
            print(self.values[0])
        TNFa.total_molecule[0] = TNFa.total_molecule[0] + self.values[0]

    def interact(self, interactable):
        #from edu.uchc.interactable.Cells import Macrophage, Afumigatus, Neutrophil, Phagocyte
        itype = type(interactable)
        if itype is TNFa:
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
        if itype is Macrophage:
            if interactable.status == Phagocyte.ACTIVE:# and interactable.state == Neutrophil.INTERACTING:
                self.inc(Constants.MA_TNF_QTTY, 0)
            if interactable.status == Phagocyte.RESTING or interactable.status == Phagocyte.ACTIVE:
                if Util.activation_function(self.get(0), Constants.Kd_TNF, Constants.STD_UNIT_T) > random():
                    interactable.status = Phagocyte.ACTIVATING if interactable.status == Phagocyte.RESTING else Phagocyte.ACTIVE
                    interactable.iteration = 0
                    interactable.tnfa = True
            return True
        if itype is Neutrophil:
            if interactable.status == Phagocyte.ACTIVE:# and interactable.state == Neutrophil.INTERACTING:
                self.inc(Constants.N_TNF_QTTY, 0)
            if interactable.status == Phagocyte.RESTING or interactable.status == Phagocyte.ACTIVE:
                if Util.activation_function(self.get(0), Constants.Kd_TNF, Constants.STD_UNIT_T) > random():
                    interactable.status = Phagocyte.ACTIVATING if interactable.status == Phagocyte.RESTING else Phagocyte.ACTIVE
                    interactable.iteration = 0
                    interactable.tnfa = True
            return True
        return interactable.interact(self)

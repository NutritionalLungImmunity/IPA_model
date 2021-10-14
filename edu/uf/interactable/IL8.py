from edu.uf.interactable.TNFa import *

class IL8(Molecule):
    name = "IL8"

    total_molecule = [0]
    NUM_STATES = 1
    threshold = -1

    def degrade(self, p=Constants.IL8_HALF_LIFE):
        super().degrade(p)
        super().degrade(Util.turnover_rate(x=1, x_system=0, k=Constants.TURNOVER_RATE, t=Constants.REL_CYT_BIND_UNIT_T, v=1))

    def get_index(self, str):
        return 0

    def inc_total_molecule(self, index, inc):
        IL8.total_molecule[index] = IL8.total_molecule[index] + inc

    def compute_total_molecule(self):
        IL8.total_molecule[0] = IL8.total_molecule[0] + self.values[0]

    def interact(self, interactable):
        #from edu.uchc.interactable.Cells import Macrophage, Afumigatus, Neutrophil, Phagocyte
        itype = type(interactable)
        if itype is IL8:
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
        if itype is Macrophage or itype is Neutrophil:
            if interactable.tnfa:#interactable.status == Phagocyte.ACTIVE and interactable.state == Neutrophil.INTERACTING:
                self.inc(0, 0)
            if (interactable.status == Phagocyte.RESTING or interactable.status == Phagocyte.ACTIVE) and type(interactable) is Neutrophil:
                if Util.activation_function(self.get(0), Constants.Kd_IL8, Constants.STD_UNIT_T) > random():
                    interactable.status = Phagocyte.ACTIVE
                    interactable.iteration = 0
            return True
        return interactable.interact(self)

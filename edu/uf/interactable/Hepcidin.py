from edu.uf.interactable.Pneumocytes import *


class Hepcidin(Molecule):

    name = "Hepcidin"

    total_molecule = [0]
    NUM_STATES = 1
    threshold = -1

    def degrade(self, p=0):
        pass
        #HEP turnover is done on "liver"

    def get_index(self, str):
        return 0

    def inc_total_molecule(self, index, inc):
        Hepcidin.total_molecule[index] = Hepcidin.total_molecule[index] + inc

    def compute_total_molecule(self):
        Hepcidin.total_molecule[0] = Hepcidin.total_molecule[0] + self.values[0]

    def interact(self, interactable):
        #from edu.uchc.interactable.Cells import Macrophage, Afumigatus, Neutrophil, Pneumocytes, Phagocyte
        itype = type(interactable)
        if itype is Hepcidin:
            return False
        if itype is IL10:
            return False
        if itype is Iron:
            return False
        if itype is Transferrin:
            return False
        if itype is TAFC:
            #if Afumigatus.total_cells > 100 or self.f:
            #    self.f = True
            #    interactable.pdec(0.9, "TAFC")
            #    interactable.pdec(0.9, "TAFCBI")
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
            return False
        #if type(interactable) is Hepatocytes:
        #    if interactable.status == Phagocyte.ACTIVE:
        #        self.inc(Constants.HEP_QTTY)
        #    return True
        if itype is Macrophage:
            #if Util.hillProbability(self.get(0)) > random():
            if Util.activation_function(self.get(0), Constants.Kd_Hep, Constants.STD_UNIT_T) > random():
                interactable.fpn = False
                interactable.fpn_iteration = 0
            return True
        return interactable.interact(self)
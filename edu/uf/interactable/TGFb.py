from edu.uf.interactable.Hepcidin import *

class TGFb(Molecule):
    name = "TGFb"

    total_molecule = [0]
    NUM_STATES = 1
    threshold = -1

    def degrade(self, p=Constants.TGF_HALF_LIFE):
        super().degrade(p)
        super().degrade(Util.turnover_rate(x=1, x_system=0, k=Constants.TURNOVER_RATE, t=Constants.REL_CYT_BIND_UNIT_T, v=1))

    def get_index(self, str):
        return 0

    def inc_total_molecule(self, index, inc):
        TGFb.total_molecule[index] = TGFb.total_molecule[index] + inc

    def compute_total_molecule(self):
        TGFb.total_molecule[0] = TGFb.total_molecule[0] + self.values[0]

    def interact(self, interactable):
        #from edu.uchc.interactable.Cells import Macrophage, Afumigatus, Neutrophil, Pneumocytes, Phagocyte
        itype = type(interactable)
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
            return False
        #if type(interactable) is Hepatocytes:
        #    return False
        if itype is Macrophage:
            if interactable.status == Phagocyte.INACTIVE:
                self.inc(Constants.MA_TGF_QTTY, 0)
                if Util.activation_function(self.get(0), Constants.Kd_TGF, Constants.STD_UNIT_T) > random():
                    interactable.iteration = 0
            elif interactable.status != Phagocyte.APOPTOTIC and interactable.status != Phagocyte.NECROTIC and interactable.status != Phagocyte.DEAD:
                if Util.activation_function(self.get(0), Constants.Kd_TGF, Constants.STD_UNIT_T) > random():
                    interactable.status = Phagocyte.INACTIVATING
            return True
        return interactable.interact(self)
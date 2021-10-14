from edu.uf.interactable.AntiTNFa import *

class Hemolysin(Molecule):
    name = "Hemolysin"

    total_molecule = [0]
    NUM_STATES = 1
    threshold = -1
    choose = None

    def degrade(self, p=Constants.ANTI_TNF_HALF_LIFE):
        super().degrade(Util.turnover_rate(x=self.values[0], x_system=0,k=Constants.TURNOVER_RATE, t=Constants.REL_CYT_BIND_UNIT_T, v=1))

    def get_index(self, str):
        return 0

    def inc_total_molecule(self, index, inc):
        Hemolysin.total_molecule[index] = Hemolysin.total_molecule[index] + inc

    def compute_total_molecule(self):
        if self.values[0] < 0:
            print(self.values[0])
        Hemolysin.total_molecule[0] = Hemolysin.total_molecule[0] + self.values[0]

    def interact(self, interactable):
        itype = type(interactable)
        if itype is Hemolysin:
            return False
        elif itype is AntiTNFa:
            return False
        elif itype is Liver:
            return False
        elif itype is Macrophage:
            return False
        elif itype is Neutrophil:
            return False
        elif itype is Pneumocytes:
            return False
        elif itype is Afumigatus:
            if interactable.status == Afumigatus.HYPHAE and Afumigatus.iter > Constants.HEMORRHAGE_DELAY:
                self.inc(Constants.HEMOLYSIN_QTTY)
            return True
        elif itype is Iron:
            return False
        elif itype is Transferrin:
            return False
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
            return False
        elif itype is Hepcidin:
            return False
        return interactable.interact(self)
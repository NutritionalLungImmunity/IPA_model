from edu.uf.interactable.Liver import *

class AntiTNFa(Molecule):
    name = "AntiTNFa"

    total_molecule = [0]
    NUM_STATES = 1
    threshold = 0
    choose = None

    def degrade(self, p=Constants.ANTI_TNF_HALF_LIFE):
        #super().degrade(p)
        if AntiTNFa.choose is None or self.id == AntiTNFa.choose:
            Constants.ANTI_TNFA_SYSTEM_CONCENTRATION = Constants.ANTI_TNFA_SYSTEM_CONCENTRATION*p
            AntiTNFa.choose = self.id
        super().degrade(
            Util.turnover_rate(x=self.values[0], x_system=Constants.ANTI_TNFA_SYSTEM_CONCENTRATION, k=Constants.TURNOVER_RATE, t=Constants.REL_CYT_BIND_UNIT_T, v=1)
        )

    def get_index(self, str):
        return 0

    def inc_total_molecule(self, index, inc):
        AntiTNFa.total_molecule[index] = AntiTNFa.total_molecule[index] + inc

    def compute_total_molecule(self):
        if self.values[0] < 0:
            print(self.values[0])
        AntiTNFa.total_molecule[0] = AntiTNFa.total_molecule[0] + self.values[0]

    def interact(self, interactable):
        itype = type(interactable)
        if itype is AntiTNFa:
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
            return False
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
            v = Util.michaelianKinetics(self.values[0], interactable.values[0], km=Constants.K_M_ANTI_TNFA, h=Constants.ANTI_TNFA_REACT_TIME_UNIT)
            #print(self.values[0], interactable.values[0], v)
            v = v if v < min(self.values[0], interactable.values[0]) else min(self.values[0], interactable.values[0])
            self.dec(v)
            interactable.dec(v)
            return True
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

    def min(self, a, b):
        return a if a < b else b
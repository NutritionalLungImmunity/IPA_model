from edu.uf.interactable.Erythrocyte import *

class EstB(Molecule):
    name = "EstB"

    total_molecule = [0]
    NUM_STATES = 1
    threshold = 0
    choose = None

    def _init(self):
        self.iron_buffer = 0

    def degrade(self, p=None):
        # super().degrade(p)
        if EstB.choose is None or self.id == EstB.choose:
            Constants.ESTB_SYSTEM_CONCENTRATION = Constants.ESTB_SYSTEM_CONCENTRATION * Constants.ESTB_HALF_LIFE
            EstB.choose = self.id
        super().degrade(
            Util.turnover_rate(x=self.values[0], x_system=Constants.ESTB_SYSTEM_CONCENTRATION,
                               k=Constants.TURNOVER_RATE, t=Constants.REL_CYT_BIND_UNIT_T, v=1)
        )

    def get_index(self, str):
        return 0

    def inc_total_molecule(self, index, inc):
        EstB.total_molecule[index] = EstB.total_molecule[index] + inc

    def compute_total_molecule(self):
        if self.values[0] < 0:
            print(self.values[0])
        EstB.total_molecule[0] = EstB.total_molecule[0] + self.values[0]

    def interact(self, interactable):
        itype = type(interactable)
        if itype is EstB:
            return False
        elif itype is Erythrocyte:
            return False
        elif itype is Hemoglobin:
            return False
        elif itype is Hemolysin:
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
            return False
        elif itype is Iron:
            interactable.inc(self.iron_buffer)
            self.iron_buffer = 0
            return True
        elif itype is Transferrin:
            return False
        elif itype is TAFC:
            v1 = Util.michaelianKinetics(interactable.values[0], self.values[0], km=Constants.ESTB_KM, Kcat=Constants.ESTB_KCAT, h = Constants.STD_UNIT_T)
            v2 = Util.michaelianKinetics(interactable.values[1], self.values[0], km=Constants.ESTB_KM, Kcat=Constants.ESTB_KCAT, h = Constants.STD_UNIT_T)
            #if self.values[0] != 0:
            #    print(v1, v2)
            interactable.dec(v1, 0)
            interactable.dec(v2, 1)
            self.iron_buffer = v2
            return True
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
from edu.uf.interactable.EstB import *

class Hemopexin(Molecule):
    name = "Hemopexin"

    total_molecule = [0]
    NUM_STATES = 1
    threshold = 0
    choose = None

    def _init(self):
        self.iron_buffer = 0

    def degrade(self, p=None):
        # super().degrade(p)
        if Hemopexin.choose is None or self.id == Hemopexin.choose:
            Constants.HEMOPEXIN_SYSTEM_CONCENTRATION = Constants.HEMOPEXIN_SYSTEM_CONCENTRATION * Constants.HEMOPEXIN_HALF_LIFE
            Hemopexin.choose = self.id
        super().degrade(
            Util.turnover_rate(x=self.values[0], x_system=Constants.HEMOPEXIN_SYSTEM_CONCENTRATION,
                               k=Constants.TURNOVER_RATE, t=Constants.REL_CYT_BIND_UNIT_T, v=1)
        )

    def get_index(self, str):
        return 0

    def inc_total_molecule(self, index, inc):
        Hemopexin.total_molecule[index] = Hemopexin.total_molecule[index] + inc

    def compute_total_molecule(self):
        if self.values[0] < 0:
            print(self.values[0])
        Hemopexin.total_molecule[0] = Hemopexin.total_molecule[0] + self.values[0]

    def interact(self, interactable):
        itype = type(interactable)
        if itype is Hemopexin:
            return False
        if itype is EstB:
            return False
        elif itype is Erythrocyte:
            return False
        elif itype is Hemoglobin:
            v = Util.michaelianKinetics(self.values[0], interactable.values[0], km=Constants.HEMOPEXIN_KM,
                                        h=Constants.STD_UNIT_T, Kcat=Constants.HEMOPEXIN_KCAT)

            v = v if v < min(self.values[0], interactable.values[0]) else 0.9999999*min(self.values[0], interactable.values[0])

            self.dec(v)
            interactable.dec(v)
            return True
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
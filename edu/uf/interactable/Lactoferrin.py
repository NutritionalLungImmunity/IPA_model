from edu.uf.interactable.Neutrophil import *

class Lactoferrin(Molecule):
    name = "Lactoferrin"

    total_molecule = [0, 0, 0]
    INDEXES = {"Lactoferrin": 0, "LactoferrinFe": 1, "LactoferrinFe2": 2}
    NUM_STATES = 3
    threshold = Constants.K_M_TF_LAC * Constants.VOXEL_VOL / 1.0e6

    def degrade(self, p=-1):
        super().degrade(Util.turnover_rate(x=1, x_system=0, k=Constants.TURNOVER_RATE, t=Constants.REL_CYT_BIND_UNIT_T, v=1))

    def get_index(self, str):
        return Lactoferrin.INDEXES.get(str)

    def inc_total_molecule(self, index, inc):
        Lactoferrin.total_molecule[index] = Lactoferrin.total_molecule[index] + inc

    def compute_total_molecule(self):
        Lactoferrin.total_molecule[0] = Lactoferrin.total_molecule[0] + self.values[0]
        Lactoferrin.total_molecule[1] = Lactoferrin.total_molecule[1] + self.values[1]
        Lactoferrin.total_molecule[2] = Lactoferrin.total_molecule[2] + self.values[2]

    def interact(self, interactable):
        #from edu.uchc.interactable.Cells import Macrophage, Afumigatus, Neutrophil
        itype = type(interactable)
        if itype is Lactoferrin:
            return False
        if itype is Macrophage: #ADD UPTAKE
            qttyFe2 = self.get("LactoferrinFe2") * Constants.MA_IRON_IMPORT_RATE * Constants.REL_IRON_IMP_EXP_UNIT_T
            qttyFe = self.get("LactoferrinFe") * Constants.MA_IRON_IMPORT_RATE * Constants.REL_IRON_IMP_EXP_UNIT_T

            qttyFe2 = qttyFe2 if qttyFe2 < self.get("LactoferrinFe2") else self.get("LactoferrinFe2")
            qttyFe = qttyFe if qttyFe < self.get("LactoferrinFe") else self.get("LactoferrinFe")

            self.dec(qttyFe2, "LactoferrinFe2")
            self.dec(qttyFe, "LactoferrinFe")
            interactable.inc_iron_pool(2 * qttyFe2 + qttyFe)
            return True
        if itype is Afumigatus:
            return False
        if itype is Neutrophil:
            if interactable.status == Neutrophil.ACTIVE and interactable.state == Neutrophil.INTERACTING:
                self.inc(Constants.LAC_QTTY, "Lactoferrin")
            return True
        if itype is TAFC:
            return False
        if itype is Transferrin:
            dfe2dt = Util.michaelianKinetics(interactable.get("TfFe2"), self.get("Lactoferrin"), Constants.K_M_TF_LAC, Constants.STD_UNIT_T)
            dfedt  = Util.michaelianKinetics(interactable.get("TfFe"), self.get("Lactoferrin"), Constants.K_M_TF_LAC, Constants.STD_UNIT_T)

            dfe2dt_fe = Util.michaelianKinetics(interactable.get("TfFe2"), self.get("LactoferrinFe"), Constants.K_M_TF_LAC, Constants.STD_UNIT_T)
            dfedt_fe = Util.michaelianKinetics(interactable.get("TfFe"), self.get("LactoferrinFe"), Constants.K_M_TF_LAC, Constants.STD_UNIT_T)

            if dfe2dt + dfedt > self.get("Lactoferrin"):
                rel = self.get("Lactoferrin") / (dfe2dt + dfedt)
                dfe2dt = dfe2dt * rel
                dfedt = dfedt * rel

            if dfe2dt_fe + dfedt_fe > self.get("LactoferrinFe"):
                rel = self.get("LactoferrinFe") / (dfe2dt_fe + dfedt_fe)
                dfe2dt_fe = dfe2dt_fe * rel
                dfedt_fe = dfedt_fe * rel

            interactable.dec(dfe2dt + dfe2dt_fe, "TfFe2")
            interactable.inc(dfe2dt + dfe2dt_fe, "TfFe")

            interactable.dec(dfedt + dfedt_fe, "TfFe")
            interactable.inc(dfedt + dfedt_fe, "Tf")

            self.dec(dfe2dt + dfedt, "Lactoferrin")
            self.inc(dfe2dt + dfedt, "LactoferrinFe")

            self.dec(dfe2dt_fe + dfedt_fe, "LactoferrinFe")
            self.inc(dfe2dt_fe + dfedt_fe, "LactoferrinFe2")

            return True
        if itype is Iron:
            qtty = interactable.get("Iron")
            qtty = qtty if qtty <= 2 * self.get("Lactoferrin") + self.get("LactoferrinFe") else 2 * self.get("Lactoferrin") + self.get("LactoferrinFe")
            rel_TfFe = Util.iron_tf_reaction(qtty, self.get("Lactoferrin"), self.get("LactoferrinFe"))
            tffe_qtty = rel_TfFe * qtty
            tffe2_qtty = (qtty - tffe_qtty) / 2
            self.dec(tffe_qtty + tffe2_qtty, "Lactoferrin")
            self.inc(tffe_qtty, "LactoferrinFe")
            self.inc(tffe2_qtty, "LactoferrinFe2")
            interactable.dec(qtty, "Iron")
            return True

        if itype is ROS:
            return False
        return interactable.interact(self)

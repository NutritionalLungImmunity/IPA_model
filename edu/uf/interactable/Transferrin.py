from edu.uf.interactable.Iron import *

class Transferrin(Molecule):
    name = "Transferrin"

    total_molecule = [0, 0, 0]
    #total_transferrinBI = 0

    INDEXES = {"Tf": 0, "TfFe": 1, "TfFe2": 2}
    NUM_STATES = 3
    threshold = Constants.K_M_TF_TAFC *Constants.VOXEL_VOL / 1.0e6

    def degrade(self, p=-1):
        pass #turnover done on "liver"

    def compute_total_molecule(self):
        Transferrin.total_molecule[0] = Transferrin.total_molecule[0] + self.values[0]
        Transferrin.total_molecule[1] = Transferrin.total_molecule[1] + self.values[1]
        Transferrin.total_molecule[2] = Transferrin.total_molecule[2] + self.values[2]

    def get_index(self, str):
        return Transferrin.INDEXES[str]

    def inc_total_molecule(self, index, inc):
        Transferrin.total_molecule[index] = Transferrin.total_molecule[index] + inc

    def interact(self, interactable):
        #from edu.uchc.interactable.Cells import Macrophage
        itype = type(interactable)
        if itype is Transferrin:
            return False
        elif itype is Macrophage:
            qttyFe2 = self.get("TfFe2") * Constants.MA_IRON_IMPORT_RATE * Constants.REL_IRON_IMP_EXP_UNIT_T
            qttyFe  = self.get("TfFe")  * Constants.MA_IRON_IMPORT_RATE * Constants.REL_IRON_IMP_EXP_UNIT_T

            qttyFe2 = qttyFe2 if qttyFe2 < self.get("TfFe2") else self.get("TfFe2")
            qttyFe = qttyFe if qttyFe < self.get("TfFe") else self.get("TfFe")

            self.dec(qttyFe2, "TfFe2")
            self.dec(qttyFe, "TfFe")
            self.inc(qttyFe2 + qttyFe, "Tf");
            interactable.inc_iron_pool(2 * qttyFe2 + qttyFe)
            if interactable.fpn and interactable.status != Macrophage.ACTIVE and interactable.status != Macrophage.ACTIVATING:
                qtty = interactable.iron_pool * self.get("Tf") * Constants.MA_IRON_EXPORT_RATE * Constants.REL_IRON_IMP_EXP_UNIT_T
                #qtty = qttyFe + 2 * qttyFe2
                qtty = qtty if qtty <= 2*self.get("Tf") else 2*self.get("Tf")
                rel_TfFe = Util.iron_tf_reaction(qtty, self.get("Tf"), self.get("TfFe"))
                tffe_qtty  = rel_TfFe*qtty
                tffe2_qtty = (qtty - tffe_qtty)/2
                self.dec(tffe_qtty + tffe2_qtty, "Tf")
                self.inc(tffe_qtty, "TfFe")
                self.inc(tffe2_qtty, "TfFe2")
                interactable.inc_iron_pool(-qtty)
            return True
        elif itype is Iron:
            qtty = interactable.get("Iron")
            qtty = qtty if qtty <= 2 * self.get("Tf") + self.get("TfFe") else  2 * self.get("Tf") + self.get("TfFe")
            rel_TfFe = Util.iron_tf_reaction(qtty, self.get("Tf"), self.get("TfFe"))
            tffe_qtty = rel_TfFe * qtty
            tffe2_qtty = (qtty - tffe_qtty) / 2
            self.dec(tffe_qtty + tffe2_qtty, "Tf")
            self.inc(tffe_qtty, "TfFe")
            self.inc(tffe2_qtty, "TfFe2")
            interactable.dec(qtty, "Iron")

            return True
        return interactable.interact(self)
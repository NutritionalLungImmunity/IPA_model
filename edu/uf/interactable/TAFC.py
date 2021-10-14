from edu.uf.interactable.Afumigatus import *

class TAFC(Molecule):
    name = "TAFC"

    def degrade(self, p=-1):
        super().degrade(Util.turnover_rate(x=1, x_system=0, k=Constants.TURNOVER_RATE, t=Constants.REL_CYT_BIND_UNIT_T, v=1))

    total_molecule = [0, 0]
    #total_tafcbi = 0
    INDEXES = {"TAFC": 0, "TAFCBI": 1}
    NUM_STATES = 2
    threshold = Constants.K_M_TF_TAFC * Constants.VOXEL_VOL / 1.0e6

    def get_index(self, str):
        return TAFC.INDEXES.get(str)

    def compute_total_molecule(self):
        TAFC.total_molecule[0] = TAFC.total_molecule[0] + self.values[0]
        TAFC.total_molecule[1] = TAFC.total_molecule[1] + self.values[1]

    def inc_total_molecule(self, index, inc):
        TAFC.total_molecule[index] = TAFC.total_molecule[index] + inc

    def interact(self, interactable):
        #from edu.uchc.interactable.Cells import Macrophage, Afumigatus
        #from random import random
        itype = type(interactable)
        if itype is TAFC:
            return False
        if itype is Macrophage:
            return False
        if itype is Transferrin:
            dfe2dt = Util.michaelianKinetics(interactable.get("TfFe2"), self.get("TAFC"), Constants.K_M_TF_TAFC, Constants.STD_UNIT_T)
            dfedt  = Util.michaelianKinetics(interactable.get("TfFe"), self.get("TAFC"), Constants.K_M_TF_TAFC, Constants.STD_UNIT_T)

            if dfe2dt + dfedt > self.get("TAFC"):
                rel = self.get("TAFC") / (dfe2dt + dfedt)
                dfe2dt = dfe2dt * rel
                dfedt = dfedt * rel

            interactable.dec(dfe2dt, "TfFe2")
            interactable.inc(dfe2dt, "TfFe")

            interactable.dec(dfedt, "TfFe")
            interactable.inc(dfedt, "Tf")

            self.inc(dfe2dt + dfedt, "TAFCBI")
            self.dec(dfe2dt + dfedt, "TAFC")

            v = dfe2dt + dfedt

            return True
        if itype is Afumigatus:
            if interactable.state == Afumigatus.FREE and interactable.status != Afumigatus.DYING and interactable.status != Afumigatus.DEAD:
                if interactable.boolean_network[Afumigatus.MirB] == 1 and interactable.boolean_network[Afumigatus.EstB] == 1:
                    qtty = self.get("TAFCBI") * Constants.TAFC_UP
                    qtty = qtty if qtty < self.get("TAFCBI") else self.get("TAFCBI")

                    self.dec(qtty, "TAFCBI")
                    interactable.inc_iron_pool(qtty)
                if interactable.boolean_network[Afumigatus.TAFC] == 1 and \
                        (interactable.status == Afumigatus.SWELLING_CONIDIA or \
                         interactable.status == Afumigatus.HYPHAE or
                         interactable.status == Afumigatus.GERM_TUBE) : # SECRETE TAFC
                    self.inc(Constants.TAFC_QTTY, "TAFC")
            return True
        if itype is Iron:
            qtty_iron = interactable.get("Iron")
            qtty_tafc = self.get("TAFC")
            qtty = qtty_tafc if qtty_tafc < qtty_iron else qtty_iron
            self.dec(qtty, "TAFC")
            self.inc(qtty, "TAFCBI")
            interactable.dec(qtty, "Iron")

            return True

        return interactable.interact(self)
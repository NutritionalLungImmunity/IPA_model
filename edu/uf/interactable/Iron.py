from edu.uf.interactable.Macrophage import *

class Iron(Molecule):

    name = "Iron"

    total_molecule = [0]
    NUM_STATES = 1
    threshold = 0

    def get_index(self, str):
        return 0

    def degrade(self, p=-1):
        pass # turnover done on "liver" (if done at all)

    def compute_total_molecule(self):
        Iron.total_molecule[0] = Iron.total_molecule[0] + self.values[0]

    def inc_total_molecule(self, index, inc):
        Iron.total_molecule[0] = Iron.total_molecule[0] + inc

    def interact(self, interactable):
        #from edu.uchc.interactable.Cells import Macrophage
        itype = type(interactable)
        if itype is Iron:
            return False
        if itype is Macrophage:
            if interactable.status == Macrophage.NECROTIC:# or interactable.status == Macrophage.APOPTOTIC or interactable.status == Macrophage.DEAD:
                self.inc(interactable.iron_pool, "Iron")
                interactable.inc_iron_pool(-interactable.iron_pool)
            return True
        return interactable.interact(self)
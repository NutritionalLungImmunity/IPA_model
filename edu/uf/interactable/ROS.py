from edu.uf.interactable.TAFC import *

class ROS(Molecule):
    name = "ROS"

    def degrade(self, p=0):
        super().degrade(p)
        #NO Turnover

    total_molecule = 0
    NUM_STATES = 1
    threshold = -1

    def get_index(self, str):
        return 0

    def inc_total_molecule(self, index, inc):
        ROS.total_molecule = ROS.total_molecule + inc

    def interact(self, interactable):
        #from edu.uchc.interactable.Cells import Afumigatus, Macrophage

        if type(interactable) is ROS:
            return False
        elif type(interactable) is Transferrin:
            return False
        elif type(interactable) is TAFC:
            return False
        elif type(interactable) is Iron:
            return False
        elif type(interactable) is Afumigatus:
        #    if interactable.state == Afumigatus.FREE and not (interactable.status == Afumigatus.DYING or interactable.status == Afumigatus.DEAD):
        #        if Util.hillProbability(self.get()) > random():
        #            interactable.status = Afumigatus.DYING
            return True
        elif type(interactable) is Macrophage:
            if interactable.status == Macrophage.ACTIVE:
                self.inc(0)
            return True
        return interactable.interact(self)
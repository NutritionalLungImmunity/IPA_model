from edu.uf.interactable.IL10 import *

class Pneumocytes(Cell):
    name = "Pneumocytes"

    total_cells = 0

    def __init__(self):
        self.id = Id.get_id()
        self.status = Phagocyte.RESTING
        self.iteration = 0
        self.tnfa = False
        Pneumocytes.total_cells = Pneumocytes.total_cells + 1

    def _init_boolean_network(self):
        pass

    def _get_max_conidia(self):
        return Constants.E_MAX_CONIDIA

    def _get_p_phagocytose(self):
        return Constants.PR_E_PHAG

    #def leave(self, qtty):
    #    return False

    def die(self):
        if self.status != Phagocyte.DEAD:
            self.status = Neutrophil.DEAD ##CAUTION!!!
            Pneumocytes.total_cells = Pneumocytes.total_cells - 1

    def interact(self, interactable):
        itype = type(interactable)
        if itype is Pneumocytes:
            return False
        if itype is Neutrophil:
            return False
        if itype is Macrophage:
            return False
        if itype is Afumigatus:
            if self.status != Phagocyte.APOPTOTIC and self.status != Phagocyte.NECROTIC and self.status != Phagocyte.DEAD:
                if interactable.status != Afumigatus.RESTING_CONIDIA:
                    if self.status != Phagocyte.ACTIVE:
                        if random() < Constants.PR_P_INT:
                            self.status = Phagocyte.ACTIVATING
                    else:
                        self.iteration = 0

                    #if self.status == Phagocyte.ACTIVE:
                    #    self.iteration = 0
                    #else:
                    #    self.status = Phagocyte.ACTIVATING
                    #Phagocyte.int_aspergillus(self, interactable)
            return True
        if itype is IL6:
            if self.status == Phagocyte.ACTIVE:
                interactable.inc(Constants.P_IL6_QTTY, 0)
            return True
        if itype is IL8:
            if self.tnfa:#self.status == Phagocyte.ACTIVE:
                interactable.inc(Constants.P_IL8_QTTY, 0)
            return True
        if itype is TNFa:
            if random() < Util.activation_function(interactable.values[0], Constants.Kd_TNF, h = Constants.UNIT_T):
                #print(Util.activation_function(interactable.values[0], Constants.Kd_TNF, h = Constants.UNIT_T))
                if self.status == Phagocyte.ACTIVE:
                    self.iteration = 0
                    self.tnfa = True
                #elif self.status == Phagocyte.RESTING:
                #    self.status = Phagocyte.ACTIVATING

            if self.status == Phagocyte.ACTIVE:
                interactable.inc(Constants.P_TNF_QTTY, 0)
            return True
        if itype is IL10:
            return False
        if itype is ROS:
            return False
        if itype is Iron:
            return False
        if itype is TAFC:
            return False
        if itype is Transferrin:
            return False
        if itype is Lactoferrin:
            return False

        return interactable.interact(self)

    def process_boolean_network(self):
        pass

    def is_dead(self):
        return self.status == Phagocyte.DEAD

    def inc_iron_pool(self, qtty):
        pass

    def update_status(self):
        if self.status == Phagocyte.ACTIVE:
            if self.iteration > Constants.ITER_TO_REST:#30*6:
                self.iteration = 0
                self.status = Phagocyte.RESTING
                self.tnfa = False
            self.iteration = self.iteration + 1
        elif self.status == Phagocyte.ACTIVATING:
            if self.iteration > Constants.ITER_TO_CHANGE_STATE:
                self.iteration = 0
                self.status = Phagocyte.ACTIVE
            self.iteration = self.iteration + 1

    #Pneumocytes don't move
    def move(self, oldVoxel, newVoxel):
        pass

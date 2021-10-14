from edu.uf.interactable.Hemoglobin import *
import numpy as np

class Erythrocyte(Cell):
    name = "Erythrocyte"

    total_cells = 0
    hemo = 0

    iter = 0

    def __init__(self, count):
        self.id = Id.get_id()
        Erythrocyte.total_cells = Erythrocyte.total_cells + count
        self.hemoglobin = 0
        self.count = count
        self.hemorrhage = False

        self.iter = 0

    def _init_boolean_network(self):
        pass

    def _get_max_conidia(self):
        return Constants.E_MAX_CONIDIA

    def _get_p_phagocytose(self):
        return Constants.PR_E_PHAG

    #def leave(self, qtty):
    #    return False

    def die(self):
        pass

    def interact(self, interactable):
        itype = type(interactable)
        if itype is Erythrocyte:
            return False
        elif itype is Hemoglobin:
            interactable.inc(self.hemoglobin)
            self.hemoglobin = 0
            return True
        elif itype is Hemolysin:
            avg = self.count*Util.activation_function(interactable.values[0], Constants.Kd_HEMO, Constants.STD_UNIT_T)
            if avg > 0:
                num = np.random.poisson(avg,1)[0]
                num = num if num < self.count else self.count
                #print(self.count, num, avg)
                self.hemoglobin = self.hemoglobin + Constants.ERYTROCYTE_HEMOGLOBIN_CONCENTRATION * num
                self.count = self.count - num
                Erythrocyte.total_cells = Erythrocyte.total_cells - num
            return True
        elif itype is AntiTNFa:
            return False
        elif itype is Liver:
            return False
        elif itype is Macrophage:
            #return False
            if self.hemorrhage and Erythrocyte.iter > Constants.HEMORRHAGE_DELAY:
                n = np.random.poisson(Constants.PR_MA_PHAG_ERYT * self.count, 1)[0] if self.count > 0 else 0
                interactable.inc_iron_pool(n * Constants.ERYTROCYTE_HEMOGLOBIN_CONCENTRATION * 4)
                self.count = self.count - n
                Erythrocyte.total_cells = Erythrocyte.total_cells - n
                interactable.status = Macrophage.ANERGIC if n > 0 else interactable.status #######
            return True
        elif itype is Neutrophil:
            return False
        elif itype is Pneumocytes:
            return False
        elif itype is Afumigatus:
            if interactable.status == Afumigatus.HYPHAE:
                if not self.hemorrhage:
                    Erythrocyte.hemo = Erythrocyte.hemo + 1
                self.hemorrhage = True
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

    def process_boolean_network(self):
        pass

    def is_dead(self):
        return False

    def inc_iron_pool(self, qtty):
        pass

    def update_status(self):
        #avg = Constants.MAX_ERYTHROCYTE_VOXEL - self.count
        avg = (1-Constants.TURNOVER_RATE)*(1-self.count/Constants.MAX_ERYTHROCYTE_VOXEL)
        if avg > 0:
            num = np.random.poisson(avg, 1)[0]
            self.count = self.count + num
            Erythrocyte.total_cells = Erythrocyte.total_cells + num

        if self.hemorrhage:
            self.iter += 1
            Erythrocyte.iter = self.iter if self.iter > Erythrocyte.iter else Erythrocyte.iter
        #if self.hemorrhage and not Constants.HEMOLYSIN:
        #    num = 0
        #    avg = Constants.HEMOLYSIS_RATE*self.count
        #    if avg > 0:
        #        num = np.random.poisson(avg, 1)[0] if num < self.count  else self.count
        #        self.count = self.count - num
        #        Erythrocyte.total_cells = Erythrocyte.total_cells - num
        #        self.hemoglobin = num * Constants.ERYTROCYTE_HEMOGLOBIN_CONCENTRATION


    def move(self, oldVoxel, newVoxel):
        pass

#from edu.uchc.interactable.Cells import *
from edu.uf.interactable.Macrophage import Macrophage
from edu.uf.interactable.Neutrophil import Neutrophil
from random import random, randint
from random import shuffle
from abc import ABC, abstractmethod
from edu.uf.interactable.Util import *


class Recruiter(ABC):

    def recruit(self, grid, quadrants):
        qtty = int(self.get_qtty())
        if qtty > 0:
            shuffle(quadrants)
            while True:
                for q in quadrants:
                    if random() < self.chemoatract(q):
                        xmin = q.x_min
                        ymin = q.y_min
                        zmin = q.z_min
                        xmax = q.x_max
                        ymax = q.y_max
                        zmax = q.z_max

                        x = randint(xmin, xmax) if xmin < xmax else xmin
                        y = randint(ymin, ymax) if ymin < ymax else ymin
                        z = randint(zmin, zmax) if zmin < zmax else zmin

                        grid[x][y][z].set_agent(self.create_cell())
                        qtty = qtty - 1
                        if qtty <= 0 :
                            return


    @abstractmethod
    def create_cell(self):
        pass

    @abstractmethod
    def get_qtty(self):
        pass

    @abstractmethod
    def chemoatract(self, voxel):
        pass

    @abstractmethod
    def leave(self):
        pass

    @abstractmethod
    def get_cell(self, voxel):
        pass

class MacrophageRecruiter(Recruiter):

    def create_cell(self):
        return Macrophage(Constants.MA_INTERNAL_IRON)

    def get_qtty(self):

        avg = Constants.MA_REC_MUL*Constants.RECRUITMENT_RATE * Macrophage.chemokine.total_molecule[0] * \
              (1 - (Macrophage.total_cells) / Constants.MAX_MA) / (Constants.Kd_MIP1B * Constants.SPACE_VOL)
        #print(Macrophage.chemokine.total_molecule[0])
        #avg = Constants.RECRUITMENT_RATE*Macrophage.chemokine.total_molecule[0]*\
        #      (1 - (Macrophage.total_cells)/Constants.MAX_MA)
        #avg = Macrophage.chemokine.total_molecule[0]/(Constants.SPACE_VOL*Constants.Kd_MIP1B)*(1 - Macrophage.total_cells/Constants.MAX_MA)
        if avg > 0:
            return np.random.poisson(avg)
        else:
            if Macrophage.total_cells < Constants.MIN_MA:
                return np.random.poisson(1)
            return 0

    def chemoatract(self, quadrant):
        return Util.activation_function(quadrant.chemokines[Macrophage.chemokine.name].get(), Constants.Kd_MIP1B, Constants.STD_UNIT_T, b=Constants.REC_BIAS)

    def leave(self):
        return True

    def get_cell(self, voxel):
        if len(voxel.cells) == 0:
            return None
        for _, agent in voxel.cells.items():
            if type(agent) is Macrophage:
                return agent
        return None

class NeutrophilRecruiter(Recruiter):

    def create_cell(self):
        return Neutrophil(0)

    def get_qtty(self):
        avg = Constants.N_REC_MUL*Constants.RECRUITMENT_RATE * Constants.N_FRAC * Neutrophil.chemokine.total_molecule[0] * \
              (1 - (Neutrophil.total_cells) / Constants.MAX_N) / (Constants.Kd_MIP2 * Constants.SPACE_VOL)
        #avg = Macrophage.chemokine.total_molecule[0] / (Constants.SPACE_VOL * Constants.Kd_MIP2) * (
        #            1 - Neutrophil.total_cells / Constants.MAX_N)
        if avg > 0:
            return np.random.poisson(avg)
        else:
            return 0

    def chemoatract(self, quadrant):
        return  Util.activation_function(quadrant.chemokines[Neutrophil.chemokine.name].get(), Constants.Kd_MIP2, Constants.STD_UNIT_T, b=Constants.REC_BIAS)

    def leave(self):
        return False

    def get_cell(self, voxel):
        pass

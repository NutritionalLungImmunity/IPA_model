from edu.uf.interactable.Util import *
from random import random
from abc import ABC, abstractmethod
import sys


class Cell(ABC):

    @abstractmethod
    def interact(self, interactable):
        pass

    @abstractmethod
    def process_boolean_network(self):
        pass

    @abstractmethod
    def update_status(self):
        pass

    @abstractmethod
    def is_dead(self):
        pass

    @abstractmethod
    def move(self, oldVoxel, newVoxel):
        pass

    #@abstractmethod
    #def leave(self, qtty):
    #    pass

    @abstractmethod
    def die(self, qtty):
        pass

    def attracted_by(self):
        return None




class Phagocyte(Cell):

    INACTIVE = 0
    INACTIVATING = 1
    RESTING = 2
    ACTIVATING = 3
    ACTIVE = 4
    APOPTOTIC = 5
    NECROTIC = 6
    DEAD = 7
    ANERGIC = 8

    FREE = 0
    INTERACTING = 1

    #total_phagocyte = 0

    def __init__(self, iron_pool = 0):
        self.id = Id.get_id()
        self._init(iron_pool) #self.boolean_network = Macrophage.InitMacrophageBooleanState;
        self.iron_pool = iron_pool
        self.status_iteration = 0
        self.state_iteration = 0
        self.phagosome = Phagosome()
        self.move_step = 0
        #Phagocyte.total_phagocyte = Phagocyte.total_phagocyte + 1

    @abstractmethod
    def _init(self):
        pass

    @abstractmethod
    def _get_max_conidia(self):
        pass

    @abstractmethod
    def get_max_move_steps(self):
        pass

    @abstractmethod
    def kill(self):
        pass

    def get_phagosome(self):
        return self.phagosome

    def add_aspergillus(self, aspergillus):
        self.phagosome.agents[aspergillus.id] = aspergillus

    def kill(self):
        if random() < Constants.PR_KILL:
            for _, a in self.phagosome.agents.items():
                self.inc_iron_pool(a.iron_pool)
                a.inc_iron_pool(-a.iron_pool)
                a.die()
                #a.status = Afumigatus.DEAD
            self.phagosome.agents.clear()

    def move(self, oldVoxel, newVoxel):
        if self.move_step < self.get_max_move_steps():
            oldVoxel.remove_cell(self.id)
            newVoxel.set_agent(self)
            self.move_step = self.move_step + 1
            for _, a in self.phagosome.agents.items():
                a.x = newVoxel.x + random()
                a.y = newVoxel.y + random()
                a.z = newVoxel.z + random()
            return True
        return False

    @staticmethod
    def int_aspergillus(phagocyte, aspergillus, phagocytize=False):
        if aspergillus.state == aspergillus.FREE:
            if (aspergillus.status == aspergillus.RESTING_CONIDIA or aspergillus.status == aspergillus.SWELLING_CONIDIA or
                aspergillus.status == aspergillus.STERILE_CONIDIA) or phagocytize:
                if phagocyte.status != Phagocyte.NECROTIC and phagocyte.status != Phagocyte.APOPTOTIC and phagocyte.status != Phagocyte.DEAD:
                    if len(phagocyte.phagosome.agents) < phagocyte._get_max_conidia():
                        phagocyte.phagosome.has_conidia = True
                        aspergillus.state = aspergillus.INTERNALIZING
                        aspergillus.engulfed = True
                        phagocyte.phagosome.agents[aspergillus.id] = aspergillus
            if aspergillus.status != aspergillus.RESTING_CONIDIA:
                phagocyte.state = Phagocyte.INTERACTING
                if phagocyte.status != Phagocyte.ACTIVE:
                    phagocyte.status = Phagocyte.ACTIVATING
                else:
                    phagocyte.status_iteration = 0

class Molecule(ABC):

    def __init__(self, qttys):
        self._init()
        self.id = Id.get_id()
        self.molecule = None
        self.f = False
        self.values = qttys #qttys should be a list
        for i in range(len(qttys)):
            self.inc_total_molecule(i, qttys[i])


    def _init(self):
        pass

    def inc(self, qtty, index=0):
        if type(index) is str:
            index = self.get_index(index)
        self.values[index] = self.values[index] + qtty
        #self.inc_total_molecule(index, qtty)
        return self.values[index]

    def dec(self, qtty, index=0):
        if type(index) is str:
            index = self.get_index(index)
        qtty = qtty if self.values[index] - qtty >= 0 else 0
        self.values[index] = self.values[index] - qtty
        #self.inc_total_molecule(index, -qtty)
        return self.values[index]

    def pdec(self, p, index=0):
        if type(index) is str:
            index = self.get_index(index)
        dec = self.values[index] * p
        self.values[index] = self.values[index] - dec
        #self.inc_total_molecule(index, -dec)
        return self.values[index]

    def pinc(self, p, index=0):
        if type(index) is str:
            index = self.get_index(index)
        inc = self.values[index] * p
        self.values[index] = self.values[index] + inc
        #self.inc_total_molecule(index, inc)
        return self.values[index]

    def set(self, qtty, index=0):
        if type(index) is str:
            index = self.get_index(index)
        self.inc_total_molecule(index, - self.values[index])
        #self.inc_total_molecule(index, qtty)
        self.values[index] = qtty

    def get(self, index=0):
        if type(index) is str:
            index = self.get_index(index)
        return self.values[index]

    def merge(self, molecule):
        if type(self) != type(molecule):
            print("molecule is not of the same type as self!", file=sys.stderr)
            sys.exit(1)

        for i in range(len(molecule._molecule)):
            self.values[i] = self.values[i] + molecule._molecule[i]

    def degrade(self, p):
        if p < 0:
            return
        for i in range(len(self.values)):
            self.pdec(1-p, i)
            #self._molecule[i] = self._molecule[i]*p

    @abstractmethod
    def compute_total_molecule(self):
        pass

    @abstractmethod
    def interact(self, interactable):
        pass

    @abstractmethod
    def get_index(self, str):
        pass

    @abstractmethod
    def inc_total_molecule(self, index, inc):
        pass


class Chemokine(Molecule):

    @abstractmethod
    def chemoatract(self):
        pass


class Phagosome():

    def __init__(self):
        self.agents = {}
        self.has_conidia = False

    def set_agent(self, agent):
        self.agents[agent.id] = agent

    def remove_agent(self, id):
        return self.agents.pop(id)

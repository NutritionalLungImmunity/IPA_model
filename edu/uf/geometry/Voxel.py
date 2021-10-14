from random import shuffle
#import random
from edu.uf.interactable.Agent import *
from edu.uf.interactable.Liver import Liver
from edu.uf.interactable.Afumigatus import Afumigatus
import timeit
import sys

class Voxel():
    AIR = 0
    EPITHELIUM = 1
    REGULAR_TISSUE = 2
    BLOOD = 3

    aspergillus_count = 0

    #xbin = None
    #zbin = None
    #ybin = None

    #Grid = None
    count_refusal = 0
    count_out_of_bounds = 0

    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z
        self.p = 0.0
        self.tissue_type = Voxel.EPITHELIUM
        self.interactables = {}
        #self.interactables = []
        self.molecules = {}
        self.cells = {}
        self.aspergillus = {}
        self.neighbors = []
        self.quadrant = None

        self.h = 0

    def set_quadrant(self, quadrant):
        quadrant.set_coordinates(self.x, self.y, self.z)
        self.quadrant = quadrant

    def set_agents(self, interactables):
        self.interactables = interactables

    def set_type(self, type):
        self.type = type

    def set_agent(self, interactable):
        if isinstance(interactable, Molecule):
            self.molecules[interactable.name] = interactable
            #    self.set_molecule(interactable.name, interactable)
        elif isinstance(interactable, Cell):
            self.cells[interactable.id] = interactable
            if type(interactable) is Afumigatus:
                self.aspergillus[interactable.id] = interactable

        elif isinstance(interactable, Liver):
            self.cells[interactable.id] = interactable
        #    self.molecules[interactable.name] = interactable
        #self.interactables[interactable.id] = interactable

        if interactable.name not in self.interactables:
            self.interactables[interactable.name] = list()
        self.interactables[interactable.name].append(interactable)


    #def get_agents(self):
    #    return self.interactables

    def get_cells(self):
        return self.cells

    def get_agent(self, id):
        return self.interactables.get(id)

    def remove_cell(self, id):
        agent = self.cells.pop(id)
        if type(agent) is Afumigatus:
            self.aspergillus.pop(id)
        return agent

    def set_molecule(self, mol_name, molecule):
        #self.set_agent(molecule)
        self.molecules[mol_name] = molecule

        if mol_name not in self.interactables:
            self.interactables[mol_name] = list()
        self.interactables[mol_name].append(molecule)

    def get_molecule(self, mol_name):
        return self.molecules[mol_name]

    def set_molecule_qtty(self, mol_name, index, qtty):
        self.molecules[mol_name].set(qtty, index)

    def set_neighbor(self, neighbor):
        self.neighbors.append(neighbor)

    def loop(self, xbin, ybin, zbin, grid):
        list = self.cells.copy()
        for _, agent in list.items():
            self._lupdate(agent)
            self.lmove(agent)
            self.lhyphal_growth(agent, xbin, ybin, zbin, grid)




        """"""

    def _lupdate(self, agent):
        if isinstance(agent, Cell):
            agent.update_status()
            agent.process_boolean_network()
            if isinstance(agent, Phagocyte):
                self._update(agent.get_phagosome().agents)
                agent.kill()

    def lmove(self, agent):
        if isinstance(agent, Cell):
            from edu.uf.interactable.Macrophage import Macrophage #CHANGE HERE!!!!
            if type(agent) is Macrophage and agent.engaged is not None: #CHANGE HERE!!!
                return #CHANGE HENRE!!!
            self._calc_drift_probability(agent)
            nvoxel = self._get_voxel(random())
            if agent.move(self, nvoxel):
                Voxel._lmove(nvoxel, agent)

    @staticmethod
    def _lmove(oldVoxel, agent):
        oldVoxel._calc_drift_probability(agent)
        nvoxel = oldVoxel._get_voxel(random())
        if agent.move(oldVoxel, nvoxel):
            Voxel._move(nvoxel, agent)

    def lhyphal_growth(self, agent, xbin, ybin, zbin, grid):
        #from edu.uchc.interactable.Cells import Phagocyte, Afumigatus
        if isinstance(agent, Phagocyte):
            tmp_phagossome = agent.get_phagosome().agents.copy()
            for _, phagent in tmp_phagossome.items():
                if isinstance(phagent, Afumigatus):
                    if phagent.status == Afumigatus.GERM_TUBE and phagent.boolean_network[
                        Afumigatus.LIP] == 1:
                        pass
                        #agent.status = Phagocyte.NECROTIC if agent.status != Phagocyte.DEAD else Phagocyte.DEAD
        elif isinstance(agent, Afumigatus):
            if agent.state == Afumigatus.FREE:
                voxel = self._find_voxel(agent.x_tip, agent.y_tip, agent.z_tip, xbin, ybin, zbin, grid)
                if voxel is not None and voxel.tissue_type != Voxel.AIR:
                    next_septa = agent.elongate()
                    if next_septa is not None:
                        voxel.set_agent(next_septa)
                    next_septa = agent.branch()
                    if next_septa is not None:
                        voxel.set_agent(next_septa)
                else:
                    Voxel.count_refusal = Voxel.count_refusal + 1
    """"""


    def interact_new(self):
        keys = list(self.interactables.keys())
        shuffle(keys)

        for i in range(len(keys)):
            for j in range(i, len(keys)):
                L1 = self.interactables.get(keys[i])
                L2 = self.interactables.get(keys[j])
                interact = True
                for l1 in L1:
                    for l2 in L2:
                        if not l1.interact(l2):
                            print(keys[i], keys[j])
                            interact = False
                            break
                    if not interact:
                        break

    def interact(self):
        cells = list(self.cells.values())
        mol = list(self.molecules.values())
        size = len(cells)
        if size > 2:
            shuffle(cells)
            shuffle(mol)
        mol_size = len(mol)
        cell_size = len(cells)
        l = [0,1,2]
        shuffle(l)
        for k in l:
            if k == 0:
                for i in range(mol_size):
                    if mol[i].threshold != -1 and mol[i].values[0] > mol[i].threshold:
                        for j in range(i, mol_size):
                            mol[i].interact(mol[j])
            elif k == 1:
                for i in range(cell_size):
                    for j in range(i, cell_size):
                        cells[i].interact(cells[j])
            elif k == 2:
                for i in range(cell_size):
                    for j in range(mol_size):
                        cells[i].interact(mol[j])

    """
    def _interact(self, agents):

        interactables = list(agents.values())

        shuffle(interactables)
        size = len(interactables)
        #if size == 5:
        #    return None
        for i in range(size):
            a1 = interactables[i]
            if isinstance(a1, Phagocyte) and a1.phagosome.is_assembled:
                self._interact(a1.phagosome.agents)
            for j in range(i, size):
                a2 = interactables[j]
                a1.interact(a2)
    """

    def update(self):
        self._update(self.interactables)
        if self.quadrant is not None:
            self.quadrant.update_chemokines(self.molecules)

    def _update(self, agents):
        for _, agent in agents.items():
            if isinstance(agent, Cell):
                agent.update_status()
                agent.process_boolean_network()
                if isinstance(agent, Phagocyte):
                    self._update(agent.get_phagosome().agents)
                    agent.kill()

    def move(self):
        tmp_interactable = self.interactables.copy()
        for _, agent in tmp_interactable.items():
            if isinstance(agent, Cell):
                self._calc_drift_probability(agent)
                nvoxel = self._get_voxel(random())
                if agent.move(self, nvoxel):
                    Voxel._move(nvoxel, agent)

    @staticmethod
    def _move(oldVoxel, agent):
        oldVoxel._calc_drift_probability(agent)
        nvoxel = oldVoxel._get_voxel(random())
        if agent.move(oldVoxel, nvoxel):
            Voxel._move(nvoxel, agent)


    def hyphal_growth(self, xbin, ybin, zbin, grid):
        #from edu.uchc.interactable.Cells import Phagocyte, Afumigatus
        tmp_interactable = self.interactables.copy()
        for _, agent in tmp_interactable.items():
            if isinstance(agent, Phagocyte):
                tmp_phagossome = agent.get_phagosome().agents.copy()
                for _, phagent in tmp_phagossome.items():
                    if isinstance(phagent, Afumigatus):
                        if phagent.status == Afumigatus.GERM_TUBE and phagent.boolean_network[Afumigatus.LIP] == 1:
                            agent.status = Phagocyte.NECROTIC
            elif isinstance(agent, Afumigatus):
                if agent.state == Afumigatus.FREE:
                    voxel =  self._find_voxel(agent.x_tip, agent.y_tip, agent.z_tip, xbin, ybin, zbin, grid)
                    if voxel is not None and voxel.tissue_type != Voxel.AIR:
                        next_septa = agent.elongate()
                        if next_septa is not None:
                            voxel.set_agent(next_septa)
                        next_septa = agent.branch()
                        if next_septa is not None:
                            voxel.set_agent(next_septa)

    def degrade(self):
        for _, m in self.molecules.items():
            m.degrade()
            m.compute_total_molecule()

    def _calc_drift_probability(self, agent):
        from edu.uf.interactable.Neutrophil import Neutrophil
        if agent.attracted_by() == None:
            return

        chemokine = self.molecules[agent.attracted_by()]
        chemoattraction = chemokine.chemoatract()

        self.p = chemoattraction + (0.0 if len(self.aspergillus) > 0 else 0.0)


        cumP = self.p
        for v in self.neighbors:
            chemokine = v.molecules[agent.attracted_by()]
            chemoattraction = chemokine.chemoatract()
            #print(chemokine.get(), chemoattraction)
            v.p = (chemoattraction  + (0.0 if len(self.aspergillus) > 0 else 0.0))if v.tissue_type != Voxel.AIR else 0.0
            cumP = cumP + v.p
        self.p = self.p / cumP
        for v in self.neighbors:
            v.p = v.p / cumP

    def _get_voxel_(self, P):
        p=0
        vx = None
        for v in self.neighbors:
            if v.p > p:
                p = v.p
                vx = v
        return vx

    def _get_voxel(self, P):
        cumP = self.p
        if P <= cumP:
            return self
        for v in self.neighbors:
            cumP = cumP + v.p
            if P <= cumP:
                return v

    #def _find_voxel(self, x, y, z, xbin, ybin, zbin, grid):
    #    if x > xbin or x < 0 or y > ybin or y < 0 or z > zbin or z < 0:
    #        Voxel.count_out_of_bounds = Voxel.count_out_of_bounds + 1
    #        return None
    #
    #    return grid[int(x)][int(y)][int(z)]

    def _find_voxel(self, x, y, z, xbin, ybin, zbin, grid):
        return grid[int(x%xbin)][int(y%ybin)][int(z%zbin)]

class Quadrant():

    def __init__(self):
        self.chemokines = {}
        self.x_min = sys.maxsize
        self.y_min = sys.maxsize
        self.z_min = sys.maxsize
        self.x_max = 0
        self.y_max = 0
        self.z_max = 0

    def update_chemokines(self, molecules):
        #from edu.uchc.interactable.Molecules import Chemokine
        for _, chemokine in molecules.items():
            if isinstance(chemokine, Chemokine):
                if chemokine.name not in self.chemokines or chemokine.values[0] > self.chemokines[chemokine.name].values[0]:
                    self.chemokines[chemokine.name] = chemokine

    def set_coordinates(self, x, y, z):
        if x < self.x_min:
            self.x_min = x
        if x > self.x_max:
            self.x_max = x
        if y < self.y_min:
            self.y_min = y
        if y > self.y_max:
            self.y_max = y
        if z < self.z_min:
            self.z_min = z
        if z > self.z_max:
            self.z_max = z


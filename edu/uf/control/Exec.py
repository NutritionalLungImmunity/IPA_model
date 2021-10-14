from edu.uf.interactable.Hepcidin import Hepcidin
from edu.uf.interactable.Hepcidin import Iron, Transferrin
from edu.uf.interactable.Afumigatus import Afumigatus
from edu.uf.interactable.Agent import Cell
from edu.uf.interactable.Agent import Phagocyte
import timeit
from edu.uf.geometry.Voxel import Voxel
#from edu.uf.geometry.Interact import Interact



class Exec():

    @staticmethod
    def exec(grid, xbin, ybin, zbin):
#        Interact.interact(grid, xbin, ybin, zbin) ### CAUTION!!!!!
        for x in range(xbin):
            for y in range(ybin):
                for z in range(zbin):
                    grid[x][y][z].interact()
                    Exec._gc(grid[x][y][z])
                    grid[x][y][z].loop(xbin, ybin, zbin, grid)
                    grid[x][y][z].quadrant.update_chemokines(grid[x][y][z].molecules)
                    #grid[x][y][z].update()
                    #grid[x][y][z].move()
                    #grid[x][y][z].hyphal_growth(xbin, ybin, zbin, grid)
                    grid[x][y][z].degrade()
    @staticmethod
    def recruit(recruiters, grid, quadrants):
        for recruiter in recruiters:
            recruiter.recruit(grid, quadrants)

    @staticmethod
    def diffusion(diffusion, grid):

        for _, molecules in grid[0][0][0].molecules.items():
            for i in range(molecules.NUM_STATES):
                if molecules.name is not Iron.name and molecules.name is not Hepcidin.name and molecules.name is not Transferrin.name:
                    diffusion[0].solver(grid, molecules.name, i)
                elif molecules.name is Transferrin.name:
                    diffusion[1].solver(grid, molecules.name, i)

    @staticmethod
    def reset_count(LIST_MOLECULE_CLASSES):
        for MOL in LIST_MOLECULE_CLASSES:
            for i in range(len(MOL.total_molecule)):
                MOL.total_molecule[i] = 0

    @staticmethod
    def _gc(voxel):
        #from edu.uf.interactable.Cells import Cell, Phagocyte

        tmp_agents = voxel.get_cells().copy()
        for k, v in tmp_agents.items():
            if isinstance(v, Cell):
                if v.is_dead():
                    if isinstance(v, Phagocyte):
                        phagosome = v.get_phagosome()
                        Exec._release_phagosome(phagosome, voxel)
                    voxel.remove_cell(k)
                #elif isinstance(v, Phagocyte):
                #    tmp_phagosome = v.get_phagosome().agents.copy()
                #    for f, a in tmp_phagosome.items():
                #        if isinstance(a, Cell):
                #            if a.is_dead():
                #                v.get_phagosome().remove_agent(f)
                elif isinstance(v, Afumigatus):
                    if v.is_internalizing():
                        voxel.remove_cell(k)

    @staticmethod
    def _release_phagosome(phagosome, voxel):
        #from edu.uf.interactable.Cells import Afumigatus
        from edu.uf.interactable.Afumigatus import Afumigatus
        for _, agent in phagosome.agents.items():
            if isinstance(agent, Afumigatus) and agent.status != Afumigatus.DEAD:
                voxel.set_agent(agent)
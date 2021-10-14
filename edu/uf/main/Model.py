#from edu.uf.interactable.Cells import *
#from edu.uf.interactable.Molecules import *
#from edu.uf.interactable.Liver import Liver
#from edu.uf.interactable.AntiTNFa import AntiTNFa
from edu.uf.interactable.Hemopexin import *
from edu.uf.geometry.Recruit import *
from edu.uf.geometry.Diffusion import *
from edu.uf.control.Exec import *
from edu.uf.interactable.MIP1B import *
from random import random, randint
import timeit
import math
import cProfile

class Model():

    SA = False
    aspergillus_count = 0
    max_count = 0
    iter = 0

    lip_on = 0

    INJ_TIME = None

    INJ_HPX = -1

    @staticmethod
    def create_grid(xbin, ybin, zbin):
        from edu.uf.geometry.Voxel import Voxel
        grid = [[[Voxel(x, y, z) for z in range(zbin)] for y in range(ybin)] for x in range(xbin)]
        for x in range(xbin):
            for y in range(ybin):
                for z in range(zbin):
                    if x + 1 < xbin: grid[x][y][z].neighbors.append(grid[x + 1][y][z])
                    if y + 1 < ybin: grid[x][y][z].neighbors.append(grid[x][y + 1][z])
                    if x - 1 >= 0: grid[x][y][z].neighbors.append(grid[x - 1][y][z])
                    if y - 1 >= 0: grid[x][y][z].neighbors.append(grid[x][y - 1][z])
                    if z + 1 < zbin: grid[x][y][z].neighbors.append(grid[x][y][z + 1])
                    if z - 1 >= 0: grid[x][y][z].neighbors.append(grid[x][y][z - 1])
        return (grid)

    @staticmethod
    def create_periodic_grid(xbin, ybin, zbin):
        from edu.uf.geometry.Voxel import Voxel
        grid = [[[Voxel(x, y, z) for z in range(zbin)] for y in range(ybin)] for x in range(xbin)]
        for x in range(xbin):
            for y in range(ybin):
                for z in range(zbin):
                    if x + 1 < xbin:
                        grid[x][y][z].neighbors.append(grid[x + 1][y][z])
                    else:
                        grid[x][y][z].neighbors.append(grid[0][y][z])
                    if y + 1 < ybin:
                        grid[x][y][z].neighbors.append(grid[x][y + 1][z])
                    else:
                        grid[x][y][z].neighbors.append(grid[x][0][z])
                    if x - 1 >= 0:
                        grid[x][y][z].neighbors.append(grid[x - 1][y][z])
                    else:
                        grid[x][y][z].neighbors.append(grid[xbin-1][y][z])
                    if y - 1 >= 0:
                        grid[x][y][z].neighbors.append(grid[x][y - 1][z])
                    else:
                        grid[x][y][z].neighbors.append(grid[x][ybin-1][z])
                    if z + 1 < zbin:
                        grid[x][y][z].neighbors.append(grid[x][y][z + 1])
                    else:
                        grid[x][y][z].neighbors.append(grid[x][y][0])
                    if z - 1 >= 0:
                        grid[x][y][z].neighbors.append(grid[x][y][z - 1])
                    else:
                        grid[x][y][z].neighbors.append(grid[x][y][zbin-1])
        return (grid)

    @staticmethod
    def create_quadrant(grid, xbin, ybin, zbin, x_size, y_size, z_size):
        from edu.uf.geometry.Voxel import Quadrant
        quadrants = []
        for i in range(math.ceil(xbin/x_size)):
            for j in range(math.ceil(ybin/y_size)):
                for k in range(math.ceil(zbin/z_size)):
                    x_min = i * x_size
                    x_max = x_min + x_size
                    y_min = j * y_size
                    y_max = y_min + y_size
                    z_min = k * z_size
                    z_max = z_min + z_size
                    q = Quadrant()
                    quadrants.append(q)
                    for x in range(x_min, x_max):
                        for y in range(y_min, y_max):
                            for z in range(z_min, z_max):
                                if x < xbin and y < ybin and z < zbin:
                                    grid[x][y][z].set_quadrant(q)

        return(quadrants)

    @staticmethod
    def initialize_molecules(grid, xbin, ybin, zbin):
        resting = 0
        for x in range(xbin):
            for y in range(ybin):
                for z in range(zbin):

                    """ Iron related """
                    grid[x][y][z].set_molecule(Iron.name, Iron([resting*random()]))
                    grid[x][y][z].set_molecule(TAFC.name, TAFC([0,0]))
                    grid[x][y][z].set_molecule(Lactoferrin.name, Lactoferrin([resting*random(),resting*random(),resting*random()]))
                    grid[x][y][z].set_molecule(Transferrin.name, Transferrin([
                        Constants.DEFAULT_APOTF_CONCENTRATION,
                        Constants.DEFAULT_TFFE_CONCENTRATION,
                        Constants.DEFAULT_TFFE2_CONCENTRATION
                    ]))
                    grid[x][y][z].set_molecule(Hepcidin.name, Hepcidin([Constants.THRESHOLD_HEP * Constants.VOXEL_VOL]))

                    """ Cytokines (all ZERO) """
                    grid[x][y][z].set_molecule(IL6.name, IL6([resting*random()]))
                    grid[x][y][z].set_molecule(TNFa.name, TNFa([resting*random()]))
                    grid[x][y][z].set_molecule(IL10.name, IL10([resting*random()]))
                    grid[x][y][z].set_molecule(TGFb.name, TGFb([resting*random()]))

                    """ Chemokines (all ZERO) """
                    grid[x][y][z].set_molecule(MIP2.name, MIP2([resting*random()]))
                    grid[x][y][z].set_molecule(MIP1B.name, MIP1B([resting*random()]))

                    """ anti-TNF-a """
                    #grid[x][y][z].set_molecule(AntiTNFa.name, AntiTNFa([Constants.ANTI_TNFA_SYSTEM_CONCENTRATION]))

                    """ hemorrhage """
                    #grid[x][y][z].set_molecule(Hemolysin.name, Hemolysin([0]))
                    #grid[x][y][z].set_molecule(Hemoglobin.name, Hemoglobin([0]))

                    """ EstB """
                    #grid[x][y][z].set_molecule(EstB.name, EstB([0]))

                    """Hemopexin"""
                    #grid[x][y][z].set_molecule(Hemopexin.name, Hemopexin([0]))


    @staticmethod
    def initialize_liver(grid, xbin, ybin, zbin):
        for x in range(xbin):
            for y in range(ybin):
                for z in range(zbin):
                    grid[x][y][z].set_agent(Liver())

    @staticmethod
    def initialize_pneumocytes(grid, xbin, ybin, zbin, num_cells=640):
        k = 0
        while k < num_cells:
            x = randint(0, xbin-1)
            y = randint(0, ybin-1)
            z = randint(0, zbin-1)
            if len(grid[x][y][z].get_cells()) == 0:
                grid[x][y][z].set_agent(Pneumocytes())
                k = k + 1

    @staticmethod
    def initialize_macrophage(grid, xmin, xmax, ymin, ymax, zmin, zmax, num_macrophages = 100):
        for _ in range(num_macrophages):
            x = randint(xmin, xmax-1)
            y = randint(ymin, ymax-1)
            z = randint(zmin, zmax-1)
            m = Macrophage(iron_pool=Constants.MA_INTERNAL_IRON)
            #m.status = Macrophage.ACTIVE
            grid[x][y][z].set_agent(m)
            #grid[x][y][z].set_agent(Neutrophil())

    def initialize_neutrophils(grid, xmin, xmax, ymin, ymax, zmin, zmax, num_macrophages = 100):
        for _ in range(num_macrophages):
            x = randint(xmin, xmax-1)
            y = randint(ymin, ymax-1)
            z = randint(zmin, zmax-1)
            grid[x][y][z].set_agent(Neutrophil())

    @staticmethod
    def initialize_erytrocyte(grid, xbin, ybin, zbin, num_erytrocytes = 320):
        for x in range(xbin):
            for y in range(ybin):
                for z in range(zbin):
                    grid[x][y][z].set_agent(Erythrocyte(num_erytrocytes))

    @staticmethod
    def infect(num_aspergillus, grid, xmin, xmax, ymin, ymax, zmin, zmax):
        l = []
        for _ in range(num_aspergillus):
            x = randint(xmin, xmax-1)
            y = randint(ymin, ymax-1)
            z = randint(zmin, zmax-1)
            a = Afumigatus(x_root=x, y_root=y, z_root=z, x_tip=x, y_tip=y, z_tip=z, ironPool = Constants.CONIDIA_INIT_IRON)
            #a.status = Afumigatus.SWELLING_CONIDIA
            #a.status = Afumigatus.GERM_TUBE
            #a.status = Afumigatus.HYPHAE
            l.append(a)
            grid[x][y][z].set_agent(a)
        return(l)

    @staticmethod
    def inject_hemopexin(grid, xbin, ybin, zbin, hemopexin_concentration):
        for x in range(xbin):
            for y in range(ybin):
                for z in range(zbin):
                    grid[x][y][z].get_molecule(Hemopexin.name).values[0] = hemopexin_concentration


    @staticmethod
    def run(iterations, xbin, ybin, zbin, grid, quadrants, recruiters, diffusion, list, input, print_lattice=False):
        l = [0,1,2]
        for k in range(iterations):
            Model.iter = k
            if k != 0:
                #pass
                shuffle(l)
            #tic = timeit.default_timer()
            if Afumigatus.total_cells > 1.1e5:
                return
            #tic = timeit.default_timer()
            for ii in l:

                if ii == 0:
                    #tic0 = timeit.default_timer()
                    Exec.exec(grid, xbin, ybin, zbin)
                    #toc0 = timeit.default_timer()
                    #print("exec ", (toc - tic))
                elif ii == 1:
                    #tic1 = timeit.default_timer()
                    Exec.recruit(recruiters, grid, quadrants)
                    #toc1 = timeit.default_timer()
                    #print("recruit ", (toc - tic))
                elif ii == 2:
                    #tic2 = timeit.default_timer()
                    Exec.diffusion(diffusion, grid)
                    #toc2 = timeit.default_timer()
            #toc = timeit.default_timer()
            #print((toc0 - tic0), (toc1 - tic1), (toc2 - tic2))
            #print(Afumigatus.total_cells, (toc - tic))
            if not Model.SA:
                pass
                Model.print_statistics(k, input, list)
            else:
                Model.count_afumigatus()
                #Model.print_statistics(k, input) #### DELET
            Exec.reset_count([Iron, Transferrin, Lactoferrin, TAFC, Hepcidin, IL6, IL10, TNFa, TGFb, MIP2, MIP1B, AntiTNFa, Hemoglobin, Hemolysin, EstB, Hemopexin])

            if False:#k == 720:
                Constants.ANTI_TNFA_SYSTEM_CONCENTRATION = Constants.ANTI_TNFA_SYSTEM_CONCENTRATION_REF

            if k == Model.INJ_HPX:
                Model.inject_hemopexin(grid, xbin, ybin, zbin, Constants.HEMOPEXIN_SYSTEM_CONCENTRATION_REF)
                Constants.HEMOPEXIN_SYSTEM_CONCENTRATION = Constants.HEMOPEXIN_SYSTEM_CONCENTRATION_REF

            if False: #k == Model.INJ_TIME:
                Constants.ESTB_SYSTEM_CONCENTRATION = 1e-6*6.4e-11
                for x in range(xbin):
                    for y in range(ybin):
                        for z in range(zbin):
                            grid[x][y][z].get_molecule(EstB.name).values[0] = 1e-6*6.4e-11

            if False:#k > 400:
                for x in range(xbin):
                    for y in range(ybin):
                        for z in range(zbin):
                            grid[x][y][z].get_molecule(Hemoglobin.name).dec(grid[x][y][z].get_molecule(Hemoglobin.name).values[0])

            if False and k%30==0:
                resting = 0
                active = 0
                restingc = 0
                inactive = 0
                anergic = 0
                other = 0
                swelling = 0
                hyphae = 0
                aother = 0
                sterile = 0
                Model.lip_on = 0
                for x in range(xbin):
                    for y in range(ybin):
                        for z in range(zbin):
                            for _, a in grid[x][y][z].cells.items():
                                if type(a) is Macrophage:
                                    if a.status == Macrophage.RESTING or a.status == Macrophage.ACTIVATING or a.status == Macrophage.INTERACTING:
                                        resting = resting + 1
                                    elif a.status == Macrophage.ACTIVE:
                                        active = active + 1
                                    elif a.status == Macrophage.INACTIVE:
                                        inactive = inactive + 1
                                    elif a.status == Macrophage.ANERGIC:
                                        anergic = anergic + 1
                                    else:
                                        other = other + 1
                                """
                                if type(a) is Afumigatus:
                                    #Model.lip_on += a.boolean_network[Afumigatus.LIP]

                                    if a.status == Afumigatus.RESTING_CONIDIA:
                                        restingc = restingc + 1
                                    elif a.status == Afumigatus.SWELLING_CONIDIA or a.status == Afumigatus.GERM_TUBE:
                                        swelling = swelling + 1
                                    elif a.status == Afumigatus.HYPHAE:
                                        hyphae = hyphae + 1
                                    else:
                                        aother = aother + 1
                                """
                print((k, resting, active, inactive, anergic, Afumigatus.total_cells))
                                #if a.status == Afumigatus.RESTING_CONIDIA:
                                #    resting = resting + 1
                                #elif a.status == Afumigatus.SWELLING_CONIDIA or a.status == Afumigatus.GERM_TUBE:
                                #    swelling = swelling + 1
                                #elif a.status == Afumigatus.HYPHAE:
                                #    hyphae = hyphae + 1
                                #elif a.status == Afumigatus.DYING or a.status == Afumigatus.DEAD:
                                #    dead = dead + 1
                                #elif a.status == Afumigatus.STERILE_CONIDIA:
                                #    sterile = sterile + 1
                #print(swelling, hyphae, aother)
                #print((resting, swelling, hyphae, dead, sterile))

                                #if a.branchable or a.growable:
                                #    print(a.iron_pool)
                                #if type(a) is Neutrophil:
                                #    if a.state == Macrophage.FREE:
                                #        free = free + 1
                                #    elif a.status == Macrophage.INTERACTING:
                                #        internalizing = internalizing + 1
                                #    else :
                                #        other = other + 1
                #print((free, internalizing, other))

            if print_lattice:
                import random
                n = open("/Users/henriquedeassis/Documents/aspergillus_tmp_files/neutrophils" + str(k) + ".csv", "w")
                f = open("/Users/henriquedeassis/Documents/aspergillus_tmp_files/aspergillus" + str(k) + ".csv", "w")
                w = open("/Users/henriquedeassis/Documents/aspergillus_tmp_files/macrophage" + str(k) + ".csv", "w")
                for x in range(xbin):
                    for y in range(ybin):
                        for z in range(zbin):
                            for _, a in grid[x][y][z].cells.items():
                                if type(a) is Neutrophil:
                                    n.write(str(x + random.random()) + "," + str(y + random.random()) + "," + str(z + random.random()) + "\n")
                                if type(a) is Afumigatus:
                                    f.write(str(a.x_tip) + "," + str(a.y_tip) + "," + str(a.z_tip) + "\n")
                                    #f.write(str(x) + "," + str(y) + "," + str(z) + "\n")
                                if type(a) is Macrophage:
                                    w.write(str(x + random.random()) + "," + str(y + random.random()) + "," + str(z + random.random()) + "\n")

                n.close()
                f.close()
                w.close()


            """
            f Afumigatus.total_cells < 4:
                for a in list:
                    print(a, a.status, a.iron_pool, a.state)
            if Afumigatus.total_cells < 4:
                for x in range(xbin):
                    for y in range(xbin):
                        for z in range(zbin):
                            for _, m in grid[x][y][z].aspergillus.items():
                                print(m, m.status, m.iron_pool)

            """

    @staticmethod
    def has_afumigatus(voxel):
        for p in voxel.interactables:
            if type(p) is Afumigatus:
                return True
        return False



    @staticmethod
    def print_statistics(k, input, cfu=None):
        from edu.uf.geometry.Voxel import Voxel
        """
        print(str(k) + "\t" +  # 1
              str(Afumigatus.total_cells) + "\t" +  # 2
              str(Afumigatus.total_iron) + "\t" +
              str(TAFC.total_molecule[0]) + "\t" +
              str(TAFC.total_molecule[1]) + "\t" +
              str((TAFC.total_molecule[0] + TAFC.total_molecule[1]))
              )
        """

        #engulfed = 0
        #for l in cfu:
        #    if l.engulfed:
        #        engulfed += 1

        print(str(k) + "\t" +   # 1
              #str(Afumigatus.total_iron) + "\t" +
              str(Afumigatus.total_cells) + "\t" +  # 2
              #str(Model.lip_on) + "\t" +  # 2
              #str(Afumigatus.total_sterile_conidia) + "\t" +  # 2
              #str(Afumigatus.total_iron) + "\t" +  # 2
              #str(1-math.exp(-Afumigatus.total_iron/(Afumigatus.total_cells*Constants.Kd_LIP*Constants.HYPHAE_VOL))) + "\t" +  # 2
              #input[1] + "\t" +  # 2
              #input[2] + "\t" +  # 2
              #input[3] + "\t" +  # 2
              #input[4] + "\t" +  # 2
              #input[5] + "\t" +  # 2
              #input[6] + "\t" +  # 2
              #input[7] + "\t" +  # 2
              #str(Erythrocyte.hemo) + "\t" +
              #str(Voxel.count_refusal) + "\t" +  # 2

              str((TAFC.total_molecule[0] + TAFC.total_molecule[1])) + "\t" +
              str(TAFC.total_molecule[0]) + "\t" +
              str(TAFC.total_molecule[1]) + "\t" +
              str(Lactoferrin.total_molecule[0]) + "\t" +
              str(Lactoferrin.total_molecule[1]) + "\t" +
              str(Lactoferrin.total_molecule[2]) + "\t" +
              #str((Transferrin.total_molecule[1] + 2*Transferrin.total_molecule[2]) + 4*Hemoglobin.total_molecule[0]) + "\t" +
              str((Transferrin.total_molecule[0] + Transferrin.total_molecule[1] + Transferrin.total_molecule[2])) + "\t" +
              str(Transferrin.total_molecule[0]) + "\t" +
              str(Transferrin.total_molecule[1]) + "\t" +
              str(Transferrin.total_molecule[2]) + "\t" +
              str(Hepcidin.total_molecule[0]) + "\t" +
              str(TGFb.total_molecule[0]) + "\t" +
              str(IL6.total_molecule[0]) + "\t" +
              str(IL10.total_molecule[0]) + "\t" +
              #str(AntiTNFa.total_molecule[0]) + "\t" +
              #str(EstB.total_molecule[0]) + "\t" +
              #str(Constants.NEUTROPHIL_HALF_LIFE) + "\t" +
              str(TNFa.total_molecule[0]) + "\t" +
              str(MIP1B.total_molecule[0]) + "\t" +
              str(MIP2.total_molecule[0]) + "\t" +
              #str(Hemolysin.total_molecule[0]) + "\t" +  # 7
              #str(Hemoglobin.total_molecule[0]) + "\t" +  # 7
              #str(Hemopexin.total_molecule[0]) + "\t" +  # 7
              #str(Erythrocyte.total_cells) + "\t" +  # 7
              #str(engulfed) + "\t" +
              #str(Constants.ITER_TO_GROW) +  "\t" +
              str(Macrophage.total_cells) + "\t" +  # 7
              str(Neutrophil.total_cells))  # 8
        """
        print(str(k) + "\t" + # 1
              str(Afumigatus.total_cells) + "\t" +  # 2
              str(Iron.total_molecule[0]) + "\t" +  # 3
              str(Afumigatus.total_iron) + "\t" +  # 4
              str(Macrophage.total_iron) + "\t" +  # 5
              str((Transferrin.total_molecule[1] + 2 * Transferrin.total_molecule[2])) + "\t" +  # 6
              str(TAFC.total_molecule[1]) + "\t" +  # 7
              str((Afumigatus.total_iron + Transferrin.total_molecule[1] + 2 * Transferrin.total_molecule[2] +
                   TAFC.total_molecule[1] + Macrophage.total_iron + Neutrophil.total_iron + Iron.total_molecule[0] +
                   Lactoferrin.total_molecule[1] + 2 * Lactoferrin.total_molecule[2])) + "\t" +  # 8
              str(Macrophage.total_cells) + "\t" +  # 9
              str(Neutrophil.total_cells))  # 10


        print(str(k) + "\t" +
              str(Afumigatus.total_cells) + "\t" +  # 2
              str(TNFa.total_molecule[0]) + "\t" +  #3
              str(MIP2.total_molecule[0]) + "\t" +  # 4
              str(MIP1B.total_molecule[0]) + "\t" +  #5
              #str((Transferrin.total_molecule[0] + Transferrin.total_molecule[1] + Transferrin.total_molecule[2])) + "\t" +  # 4
              str(Macrophage.total_cells) + "\t" +  # 6
              str(Neutrophil.total_cells))  # 7


        print(str(Afumigatus.total_cells) + "\t" + #1
              str(Iron.total_iron) + "\t" + #2
              str(Afumigatus.total_iron) + "\t" +  #3
              str(Macrophage.total_iron) + "\t" +  #4
              str(Neutrophil.total_iron) + "\t" +  #5
              str((Transferrin.total_transferrin[0] + Transferrin.total_transferrin[1] + Transferrin.total_transferrin[2])) + "\t" +  #6
              str((Transferrin.total_transferrin[1] + 2 * Transferrin.total_transferrin[2])) + "\t" + #7
              str(TAFC.total_tafc[1]) + "\t" + #8
              str((Lactoferrin.total_lactoferrin[1] + 2 * Lactoferrin.total_lactoferrin[2])) + "\t" + #9
              str((Afumigatus.total_iron + Transferrin.total_transferrin[1] + 2 * Transferrin.total_transferrin[2] +
                  TAFC.total_tafc[1] + Macrophage.total_iron + Neutrophil.total_iron + Iron.total_iron +
                   Lactoferrin.total_lactoferrin[1] + 2 * Lactoferrin.total_lactoferrin[2])) + "\t" + #10
              str(Macrophage.total_cells) + "\t" + #11
              str(Neutrophil.total_cells)) #12

        print(str(Afumigatus.total_cells) + "\t" + str(Iron.total_iron) + "\t" + str(TAFC.total_tafc[0]) + "\t" + str(TAFC.total_tafc[1]) + "\t" + \
              str(Transferrin.total_transferrin[0]) + "\t" + str(Transferrin.total_transferrin[1]) + "\t" + \
              str(Lactoferrin.total_lactoferrin[0]) + "\t" + str(Lactoferrin.total_lactoferrin[1]) + "\t" + str(ROS.total_ros) + "\t" + \
              str((Transferrin.total_transferrin[0] + Transferrin.total_transferrin[1])) + "\t" + \
              str(Afumigatus.total_iron) + "\t" + str(Macrophage.total_iron) + "\t" + str(Neutrophil.total_iron) + "\t" + \
              str((Iron.total_iron + TAFC.total_tafc[1] + Transferrin.total_transferrin[1] + \
                 Afumigatus.total_iron + Macrophage.total_iron + Neutrophil.total_iron + Lactoferrin.total_lactoferrin[1])) + "\t" + \
              str(IL6.total_il6[0]) + "\t" + str(IL8.total_il8[0]) + "\t" + str(IL10.total_il10[0]) + "\t" + str(TNFa.total_tnfa[0]) + "\t" + str(Hepcidin.total_hep[0]) + "\t" + \
              str(Macrophage.total_cells) + "\t" + str(Neutrophil.total_cells) + "\t" + str(Pneumocytes.total_cells))
        """

    @staticmethod
    def count_afumigatus():
        Model.aspergillus_count = Model.aspergillus_count + Afumigatus.total_cells
        if Afumigatus.total_cells > Model.max_count:
            Model.max_count = Afumigatus.total_cells

    @staticmethod
    def read_parameters(parameters):
        Constants.D = float(parameters[1])
        # pr swell
        Constants.PR_ASPERGILLUS_CHANGE = float(parameters[2])  ##
        # iteration to change state
        Constants.ITER_TO_CHANGE_STATE = int(parameters[3])
        # iterations to rest
        Constants.ITER_TO_REST = int(parameters[4])
        # iterations to grow (growth rate^-1)
        Constants.ITER_TO_GROW = int(parameters[5])
        # pr branch
        Constants.PR_BRANCH = float(parameters[6])  ##
        # turnover rate (1-x)
        Constants.TURNOVER_RATE = float(parameters[7]) ######
        # LACTOFERRIN_QTTY
        Constants.LAC_QTTY = float(parameters[8])
        # TAFC uptake rate
        Constants.TAFC_UP = float(parameters[9])
        # TAFC_QTTY
        # cytokines half life
        Constants.IL10_HALF_LIFE = float(parameters[10]) #####
        Constants.IL6_HALF_LIFE = float(parameters[10])
        Constants.TNF_HALF_LIFE = float(parameters[10])
        Constants.TGF_HALF_LIFE = float(parameters[10])
        Constants.MIP2_HALF_LIFE = float(parameters[10])
        Constants.MIP1B_HALF_LIFE = float(parameters[10])
        # Kd
        Constants.Kd_HEMO = float(parameters[11])
        Constants.Kd_IL6 = float(parameters[12])
        Constants.Kd_IL10 = float(parameters[13])
        Constants.Kd_MIP1B = float(parameters[14])
        Constants.Kd_MIP2 = float(parameters[15])
        Constants.Kd_TNF = float(parameters[16])
        Constants.Kd_TGF = float(parameters[17])
        Constants.Kd_Hep = float(parameters[18])
        # Kd LIP
        Constants.Kd_LIP = float(parameters[19])
        # Macrophage iron export rate
        Constants.MA_IRON_EXPORT_RATE = float(parameters[20])
        # move rate
        Constants.MA_MOVE_RATE_REST = float(parameters[21])
        Constants.MA_MOVE_RATE_ACT = float(parameters[21]),
        # Pr phagocytosis/kill&interaction
        Constants.PR_MA_PHAG = float(parameters[22]) ######
        Constants.PR_N_PHAG = float(parameters[23])
        Constants.PR_P_INT = float(parameters[24])
        Constants.PR_N_HYPHAE = float(parameters[25])
        Constants.PR_MA_HYPHAE = float(parameters[26])
        Constants.PR_KILL = float(parameters[27])
        # K_M
        Constants.K_M_TF_TAFC = float(parameters[28])
        Constants.K_M_TF_LAC = float(parameters[29])
        # cells half life
        Constants.NEUTROPHIL_HALF_LIFE = float(parameters[30])
        Constants.MA_HALF_LIFE = float(parameters[31])
        # HEP_CURVE
        #Constants.HEP_INTERCEPT = float(parameters[37])
        Constants.HEP_SLOPE = float(parameters[32])
        # Tf CURVE
        #Constants.TF_INTERCEPT = float(parameters[39])
        Constants.TF_SLOPE = float(parameters[33])

        Constants.IL6_THRESHOLD = math.pow(10, (Constants.THRESHOLD_LOG_HEP - Constants.HEP_INTERCEPT) / Constants.HEP_SLOPE)
        Constants.DEFAULT_TF_CONCENTRATION = (Constants.TF_INTERCEPT + Constants.TF_SLOPE* Constants.THRESHOLD_LOG_HEP) * Constants.VOXEL_VOL
        Constants.DEFAULT_APOTF_CONCENTRATION = Constants.DEFAULT_APOTF_REL_CONCENTRATION * Constants.DEFAULT_TF_CONCENTRATION
        Constants.DEFAULT_TFFE_CONCENTRATION = Constants.DEFAULT_TFFE_REL_CONCENTRATION * Constants.DEFAULT_TF_CONCENTRATION
        Constants.DEFAULT_TFFE2_CONCENTRATION = Constants.DEFAULT_TFFE2_REL_CONCENTRATION * Constants.DEFAULT_TF_CONCENTRATION

        Constants.CONIDIA_INIT_IRON = Constants.Kd_LIP * Constants.CONIDIA_VOL
        #pr phagocyte erythrocyte
        Constants.PR_MA_PHAG_ERYT = float(parameters[34])
        Constants.HEMOGLOBIN_UPTAKE_RATE = float(parameters[35])

    @staticmethod
    def read_ma_hyphae_parameters(parameters):
        Constants.ITER_TO_GROW = int(parameters[1])
        Constants.PR_MA_HYPHAE = float(parameters[2])

    @staticmethod
    def read_ma_tafcup_parameters(parameters):
        Constants.TAFC_UP = float(parameters[1])
        Constants.PR_MA_HYPHAE = float(parameters[2])

    @staticmethod
    def read_EstB_parameters(parameters):
        Constants.ESTB_HALF_LIFE = float(parameters[1])
        Model.INJ_TIME = int(parameters[2]) * 30

        Constants.PR_MA_PHAG_ERYT = 0
        Constants.HEMOGLOBIN_UPTAKE_RATE = 2 * 0.01 * 0.0010416666666666667

    @staticmethod
    def read_hemo_parameter(parameters):
        Constants.PR_MA_PHAG_ERYT = 0.00702550791573453 * float(parameters[1])
        Constants.HEMOGLOBIN_UPTAKE_RATE = 0
        Constants.HEMOLYSIN_QTTY = 0

    @staticmethod
    def read_hemoH_parameter(parameters):
        Constants.PR_MA_PHAG_ERYT = 0
        Constants.HEMOGLOBIN_UPTAKE_RATE = 0.01 * 0.0010416666666666667 * float(parameters[1])

    @staticmethod
    def print_sa():
        print(
            str(Constants.D) + "," +  # = float(parameters[1])
            str(Constants.PR_ASPERGILLUS_CHANGE) + "," +  # = float(parameters[2])  ##
            str(Constants.ITER_TO_CHANGE_STATE) + "," +  # = int(parameters[3])
            str(Constants.ITER_TO_REST) + "," +  # = int(parameters[4])
            str(Constants.ITER_TO_GROW) + "," +  # = int(parameters[5])
            str(Constants.PR_BRANCH) + "," +  # = float(parameters[6])  ##
            str(Constants.TURNOVER_RATE) + "," +  # = float(parameters[7]) ######
            str(Constants.LAC_QTTY) + "," +  # = float(parameters[8])
            str(Constants.TAFC_UP) + "," +  # = float(parameters[9])
            str(Constants.IL10_HALF_LIFE) + "," +  # = float(parameters[12]) #####
            str(Constants.Kd_HEMO) + "," +  # = float(parameters[10])
            str(Constants.Kd_IL6) + "," +  # = float(parameters[14])
            str(Constants.Kd_IL10) + "," +  # = float(parameters[15])
            str(Constants.Kd_MIP1B) + "," +  # = float(parameters[16])
            str(Constants.Kd_MIP2) + "," +  # = float(parameters[17])
            str(Constants.Kd_TNF) + "," +  # = float(parameters[18])
            str(Constants.Kd_TGF) + "," +  # = float(parameters[19])
            str(Constants.Kd_Hep) + "," +  # = float(parameters[20])
            str(Constants.Kd_LIP) + "," +  # = float(parameters[21])
            str(Constants.MA_IRON_EXPORT_RATE) + "," +  # = float(parameters[24]
            str(Constants.MA_MOVE_RATE_REST) + "," +  #T = float(parameters[26])
            str(Constants.PR_MA_PHAG) + "," +  # = float(parameters[27]) ######
            str(Constants.PR_N_PHAG) + "," +  # = float(parameters[28])
            str(Constants.PR_P_INT) + "," +  # = float(parameters[29])
            str(Constants.PR_N_HYPHAE) + "," +  # = float(parameters[30])
            str(Constants.PR_MA_HYPHAE) + "," +  # = float(parameters[31])
            str(Constants.PR_KILL) + "," +  # = float(parameters[32])
            str(Constants.K_M_TF_TAFC) + "," +  # = float(parameters[33])
            str(Constants.K_M_TF_LAC) + "," +  # = float(parameters[34])
            str(Constants.NEUTROPHIL_HALF_LIFE) + "," +  # = float(parameters[35])
            str(Constants.MA_HALF_LIFE) + "," +  # = float(parameters[36])
            str(Constants.HEP_SLOPE) + "," +  # = float(parameters[37])
            str(Constants.TF_SLOPE) + "," +  # = float(parameters[38])
            str(Constants.PR_MA_PHAG_ERYT) + "," +  # =float(parameters[41])
            str(Constants.HEMOGLOBIN_UPTAKE_RATE) + "," +  # =float(parameters[41])

            str(Model.aspergillus_count) + "," +
            str(Afumigatus.total_cells) + "," +
            str(Model.max_count) + "," +
            str(Model.iter)
        )


    @staticmethod
    def main(input, output):
        Model.SA = False
        #bla()
        if Model.SA:
            Model.read_parameters(input)

        input = ["0", "1920", "15", "640"]#, input[40]]#, input[1], input[2], input[3]]

        xbin = 10
        ybin = 10
        zbin = 10
        xquadrant = 3
        yquadrant = 3
        zquadrant = 3

        grid = Model.create_periodic_grid(xbin, ybin, zbin)
        quadrants = Model.create_quadrant(grid, xbin, ybin, zbin, xquadrant, yquadrant, zquadrant)
        Model.initialize_molecules(grid, xbin, ybin, zbin)
        Model.initialize_pneumocytes(grid, xbin, ybin, zbin, int(input[3]))
        Model.initialize_liver(grid, xbin, ybin, zbin) ###
        #Model.initialize_erytrocyte(grid, xbin, ybin, zbin, Constants.MAX_ERYTHROCYTE_VOXEL)

        Model.initialize_macrophage(grid, 0, xbin, 0, ybin, 0, zbin, int(input[2]))
        list = Model.infect(int(input[1]), grid, 0, xbin, 0, ybin, 0, zbin)

        #Constants.MA_HALF_LIFE = 0
        recruiters = []
        recruiters.append(MacrophageRecruiter())
        #recruiters.append(NeutrophilRecruiter())

        #print(Constants.PR_P_INT)
        D = Constants.D/(4*(30/Constants.TIME_STEP_SIZE))
        diffusion_closed = ADI(D, 4)
        diffusion_closed.create_periodic_matrixes(xbin, ybin, zbin)
        diffusion_open = ADI(D, 4)
        diffusion_open.create_periodic_matrixes(xbin, ybin, zbin)
        #diffusion_open.create_open_matrixes(xbin, ybin, zbin)
        #int(9*(30/Constants.TIME_STEP_SIZE))  int(250*(30/Constants.TIME_STEP_SIZE))
        tic = timeit.default_timer();
        #int(72*2*(30/Constants.TIME_STEP_SIZE))
        Model.run(1440, xbin, ybin, zbin, grid, quadrants, recruiters, [diffusion_open, diffusion_closed], list, input, print_lattice=False)
        toc = timeit.default_timer();
        print((toc - tic))
        #Model.print_sa()

if __name__ == "__main__":
    Model.main(sys.argv, "")

    #for i in range(1):
    #    Model.main(str(i), str())

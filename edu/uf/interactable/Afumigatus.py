from edu.uf.interactable.Transferrin import *
import math

class Afumigatus(Cell):
    name = "Afumigatus"
    InitAfumigatusBooleanState = [1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

    RESTING_CONIDIA = 0
    SWELLING_CONIDIA = 1
    GERM_TUBE = 2
    HYPHAE = 3
    DYING = 4
    DEAD = 5
    STERILE_CONIDIA = 6

    FREE = 0
    INTERNALIZING = 1
    RELEASING = 2
    ENGAGED = 3

    hapX = 0
    sreA = 1
    HapX = 2
    SreA = 3
    RIA = 4
    EstB = 5
    MirB = 6
    SidA = 7
    TAFC = 8
    ICP = 9
    LIP = 10
    CccA = 11
    FC0fe = 12
    FC1fe = 13
    VAC = 14
    ROS = 15
    Yap1 = 16
    SOD2_3 = 17
    Cat1_2 = 18
    ThP = 19
    Fe = 20
    O = 21
    #TAFCBI = 22
    SPECIES_NUM = 22

    total_iron = 0
    total_cells = 0
    total_sterile_conidia = 0
    total_cfu = set()

    tree_septa = set()

    iter = 0


    def __init__(self, x_root=0, y_root=0, z_root=0, x_tip=0, y_tip=0, z_tip=0,
                 #dx=random(), dy=random(), dz=random(), \
                 dx=None, dy=None, dz=None,\
                 growth_iteration = 0, ironPool = 0, status = 0, state = 0, isRoot = True):

        ################################################
        if dx is None or dy is None or dz is None:
            dx, dy, dz = Util.random_sphere_point()
        ################################################

        self.id = Id.get_id()
        self.iron_pool = ironPool
        self.state = state
        self.status = status
        self.is_root = isRoot
        self.x_root = x_root
        self.y_root = y_root
        self.z_root = z_root
        self.x_tip = x_tip
        self.y_tip = y_tip
        self.z_tip = z_tip
        self.dx = dx
        self.dy = dy
        self.dz = dz
        self.engulfed = False

        self.iter_to_die = 0 # CHANGE HERE!!!!

        self.cfu = None

        self.growable = True
        self.branchable = False
        self.activation_iteration = 0
        self.growth_iteration = growth_iteration
        self.boolean_network = Afumigatus.InitAfumigatusBooleanState.copy()

        self.next_septa = None
        self.next_branch = None
        self.previous_septa = None
        #self.Fe = False
        self.bn_iteration = 0

        self._normalize()

        self.iter = 0

        Afumigatus.total_iron = Afumigatus.total_iron + ironPool
        Afumigatus.total_cells = Afumigatus.total_cells + 1

    def _normalize(self):
        S = math.sqrt(self.dx**2  + self.dy**2 + self.dz**2)
        self.dx = self.dx/S
        self.dy = self.dy/S
        self.dz = self.dz/S
        #self.x_tip = self.x_root + self.dx
        #self.y_tip = self.y_root + self.dy
        #self.z_tip = self.z_root + self.dz

    def move(self, oldVoxel, newVoxel):
        pass

    def is_internalizing(self):
        return self.state == Afumigatus.INTERNALIZING

    def elongate(self):

        septa = None
        if self.growable and self.boolean_network[Afumigatus.LIP] == 1:
            if self.status == Afumigatus.HYPHAE:
                if self.growth_iteration >= Constants.ITER_TO_GROW:
                    self.growth_iteration = 0
                    self.growable = False
                    self.branchable = True
                    self.iron_pool = self.iron_pool / 2.0;

                    self._normalize()
                    self.next_septa = Afumigatus(x_root = self.x_tip, y_root = self.y_tip, z_root = self.z_tip, \
                                                 x_tip=self.x_tip + self.dx, y_tip=self.y_tip + self.dy, z_tip=self.z_tip + self.dz,\
                                                 dx=self.dx, dy=self.dy, dz=self.dz, growth_iteration=0, \
                                                 ironPool=0, status=Afumigatus.HYPHAE, state=self.state, isRoot=False)

                    self.next_septa.previous_septa = self
                    self.next_septa.iron_pool = self.iron_pool
                    septa = self.next_septa
                else:
                    self.growth_iteration = self.growth_iteration + 1
            elif self.status == Afumigatus.GERM_TUBE:
                if self.growth_iteration >= Constants.ITER_TO_GROW:
                    self.status = Afumigatus.HYPHAE
                    self._normalize()
                    self.x_tip = self.x_root + self.dx
                    self.y_tip = self.y_root + self.dy
                    self.z_tip = self.z_root + self.dz

                else:
                    self.growth_iteration = self.growth_iteration + 1
        return septa

    def branch(self, phi=None, pr_branch=Constants.PR_BRANCH):
        branch = None
        if self.branchable and self.status == Afumigatus.HYPHAE and self.boolean_network[Afumigatus.LIP] == 1:
            if random() < 1:#pr_branch:
                if phi == None:
                    phi = 2*random()*math.pi

                self.iron_pool = self.iron_pool / 2.0
                growth_vector = np.array([self.dx, self.dy, self.dz])
                B = np.array([growth_vector, [random(), random(), random()], [random(), random(), random()]]).transpose()
                B = Util.gram_schmidt(B)
                B_inv = np.transpose(B)
                R = Util.rotatation_matrix(phi)
                growth_vector = B.dot(R.dot(B_inv.dot(growth_vector.transpose())))

                self._normalize() #UNNECESSARY?
                self.next_branch = Afumigatus(x_root=self.x_tip, y_root=self.y_tip, z_root=self.z_tip, \
                                             x_tip=self.x_tip + growth_vector[0], y_tip=self.y_tip + growth_vector[1], z_tip=self.z_tip + growth_vector[2], \
                                              dx=growth_vector[0], dy=growth_vector[1], dz=growth_vector[2], growth_iteration = -1, \
                                              ironPool=0, status=Afumigatus.HYPHAE, state=self.state, isRoot=False)

                self.next_branch.iron_pool = self.iron_pool
                self.next_branch.previous_septa = self
                branch = self.next_branch
            self.branchable = False
        return branch


    def interact(self, interactable):
        itype = type(interactable)
        if itype is Afumigatus:
            return False
        elif itype is Iron:
            if self.status == Afumigatus.DYING or self.status == Afumigatus.DEAD:
                interactable.inc(self.iron_pool)
                self.inc_iron_pool(-self.iron_pool)
            return True
        elif itype is Transferrin:
            return False
        elif itype is Macrophage:
            if interactable.engaged:
                return True
            if interactable.status != Macrophage.APOPTOTIC and interactable.status != Macrophage.NECROTIC and interactable.status != Macrophage.DEAD:
                if self.status != Afumigatus.RESTING_CONIDIA:
                    pr_interact = Constants.PR_MA_HYPHAE if self.status == Afumigatus.HYPHAE else Constants.PR_MA_PHAG
                    if random() < pr_interact:
                        Phagocyte.int_aspergillus(interactable, self, self.status != Afumigatus.HYPHAE)
                        if self.status == Afumigatus.HYPHAE and interactable.status == Macrophage.ACTIVE:
                            self.status = Afumigatus.DYING
                            if self.next_septa is not None:
                                self.next_septa.is_root = True
                            if self.next_branch is not None:
                                self.next_branch.is_root = True
                        else:
                            if self.status == Afumigatus.HYPHAE and interactable.status == Macrophage.ACTIVE:
                                interactable.engaged = True

            return True
        return interactable.interact(self)

    """
    def interact(self, interactable):
        itype = type(interactable)
        if itype is Afumigatus:
            return False
        elif itype is Iron:
            if self.status == Afumigatus.DYING or self.status == Afumigatus.DEAD:
                interactable.inc(self.iron_pool)
                self.inc_iron_pool(-self.iron_pool)
            return True
        elif itype is Transferrin:
            return False
        elif itype is Macrophage:
            if interactable.engaged is not None:
                #print((interactable, interactable.engaged, self))
                if interactable.engaged is self:
                    if self.iter_to_die < Constants.ITERATION_TO_KILL:
                        self.iter_to_die += 1
                    else:
                        #print(self.iter_to_die)
                        self.iter_to_die = 0
                        interactable.engaged.status = Afumigatus.DYING

                        if interactable.engaged.next_septa is not None:
                            interactable.engaged.next_septa.is_root = True
                        if interactable.engaged.next_branch is not None:
                            interactable.engaged.next_branch.is_root = True
                        interactable.engaged = None
                return True

            if interactable.status != Macrophage.APOPTOTIC and interactable.status != Macrophage.NECROTIC and interactable.status != Macrophage.DEAD:
                if self.status != Afumigatus.RESTING_CONIDIA:
                    pr_interact = Constants.PR_MA_HYPHAE if self.status == Afumigatus.HYPHAE else Constants.PR_MA_PHAG
                    if random() < pr_interact:
                        Phagocyte.int_aspergillus(interactable, self, self.status != Afumigatus.HYPHAE)
                        if self.status == Afumigatus.HYPHAE and interactable.status == Macrophage.ACTIVE:
                            interactable.engaged = self
                            self.iter_to_die = 1
                            #print(self.iter_to_die)

            return True
        return interactable.interact(self)
    """

    def process_boolean_network(self):
        if self.bn_iteration == 15:
            #print(self.iron_pool)
            temp = [0 for _ in range(Afumigatus.SPECIES_NUM)]

            temp[Afumigatus.hapX] = -self.boolean_network[Afumigatus.SreA] + 1
            temp[Afumigatus.sreA] = -self.boolean_network[Afumigatus.HapX] + 1
            temp[Afumigatus.HapX] = self.boolean_network[Afumigatus.hapX] & (-self.boolean_network[Afumigatus.LIP] + 1)
            temp[Afumigatus.SreA] = self.boolean_network[Afumigatus.sreA] & self.boolean_network[Afumigatus.LIP]
            temp[Afumigatus.RIA] = -self.boolean_network[Afumigatus.SreA] + 1
            temp[Afumigatus.EstB] = -self.boolean_network[Afumigatus.SreA] + 1
            temp[Afumigatus.MirB] = self.boolean_network[Afumigatus.HapX] & (-self.boolean_network[Afumigatus.SreA] + 1)
            temp[Afumigatus.SidA] = self.boolean_network[Afumigatus.HapX] & (-self.boolean_network[Afumigatus.SreA] + 1)
            temp[Afumigatus.TAFC] = self.boolean_network[Afumigatus.SidA]
            temp[Afumigatus.ICP] = (-self.boolean_network[Afumigatus.HapX] + 1) & (self.boolean_network[Afumigatus.VAC] | self.boolean_network[Afumigatus.FC1fe])
            temp[Afumigatus.LIP] = (self.boolean_network[Afumigatus.Fe] & self.boolean_network[Afumigatus.RIA]) | self._lip_activation()
            temp[Afumigatus.CccA] = -self.boolean_network[Afumigatus.HapX] + 1
            temp[Afumigatus.FC0fe] = self.boolean_network[Afumigatus.SidA]
            temp[Afumigatus.FC1fe] = self.boolean_network[Afumigatus.LIP] & self.boolean_network[Afumigatus.FC0fe]
            temp[Afumigatus.VAC] = self.boolean_network[Afumigatus.LIP] & self.boolean_network[Afumigatus.CccA]
            # temp[Afumigatus.ROS] = self.boolean_network[Afumigatus.LIP] | \
            #                        (self.boolean_network[Afumigatus.O] & (- (self.boolean_network[Afumigatus.SOD2_3] & self.boolean_network[Afumigatus.ThP] \
            #                         & self.boolean_network[Afumigatus.Cat1_2]) + 1)) \
            #                        | (self.boolean_network[Afumigatus.ROS] & (- (self.boolean_network[Afumigatus.SOD2_3] \
            #                         & (self.boolean_network[Afumigatus.ThP] | self.boolean_network[Afumigatus.Cat1_2])) + 1))
            temp[Afumigatus.ROS] = (self.boolean_network[Afumigatus.O] & (- (self.boolean_network[Afumigatus.SOD2_3] & self.boolean_network[Afumigatus.ThP] \
                                   & self.boolean_network[Afumigatus.Cat1_2]) + 1)) \
                                  | (self.boolean_network[Afumigatus.ROS] & (- (self.boolean_network[Afumigatus.SOD2_3] \
                                   & (self.boolean_network[Afumigatus.ThP] | self.boolean_network[Afumigatus.Cat1_2])) + 1))
            temp[Afumigatus.Yap1] = self.boolean_network[Afumigatus.ROS]
            temp[Afumigatus.SOD2_3] = self.boolean_network[Afumigatus.Yap1]
            temp[Afumigatus.Cat1_2] = self.boolean_network[Afumigatus.Yap1] & (-self.boolean_network[Afumigatus.HapX] + 1)
            temp[Afumigatus.ThP] = self.boolean_network[Afumigatus.Yap1]

            temp[Afumigatus.Fe] = 0 # might change according to iron environment?
            temp[Afumigatus.O] = 0
            #temp[Afumigatus.TAFCBI] = 0

            #print(self.boolean_network)
            for i in range(Afumigatus.SPECIES_NUM):
                self.boolean_network[i] = temp[i]
            self.bn_iteration = 0
        else:
            self.bn_iteration = self.bn_iteration + 1

    def update_status(self):
        self.activation_iteration = self.activation_iteration + 1
        if self.status == Afumigatus.RESTING_CONIDIA and \
                self.activation_iteration >= Constants.ITER_TO_SWELLING and \
                random() < Constants.PR_ASPERGILLUS_CHANGE:
            self.status = Afumigatus.SWELLING_CONIDIA
            self.activation_iteration = 0
        elif self.status == Afumigatus.SWELLING_CONIDIA and \
                self.activation_iteration >= Constants.ITER_TO_GERMINATE:# and \
                #random() < Constants.PR_ASPERGILLUS_CHANGE:
            self.status = Afumigatus.GERM_TUBE
            self.activation_iteration = 0
        elif self.status == Afumigatus.DYING:
            self.die()

        if self.next_septa == None:
            self.growable = True

        if self.state == Afumigatus.INTERNALIZING or self.state == Afumigatus.RELEASING:
            self.state = Afumigatus.FREE

        self.diffuse_iron()
        if self.next_branch is None:
            pass #REVIEW!!!
            #self.growable = True

        if self.status == Afumigatus.HYPHAE:
            self.iter += 1
            Afumigatus.iter = self.iter if self.iter > Afumigatus.iter else Afumigatus.iter

    def is_dead(self):
        return self.status == Afumigatus.DEAD

    def diffuse_iron(self, afumigatus=None):
        if afumigatus == None:
            if self.is_root:
                Afumigatus.tree_septa = set()
                Afumigatus.tree_septa.add(self)
                self.diffuse_iron(self)
                total_iron = 0
                for a in Afumigatus.tree_septa:
                    total_iron = total_iron + a.iron_pool
                total_iron = total_iron/len(Afumigatus.tree_septa)
                for a in Afumigatus.tree_septa:
                    a.iron_pool = total_iron
        else:
            if afumigatus.next_septa is not None and afumigatus.next_branch is None:
                Afumigatus.tree_septa.add(afumigatus.next_septa)
                self.diffuse_iron(afumigatus.next_septa)
            elif afumigatus.next_septa is not None and afumigatus.next_branch is not None:
                Afumigatus.tree_septa.add(afumigatus.next_septa)
                Afumigatus.tree_septa.add(afumigatus.next_branch)
                self.diffuse_iron(afumigatus.next_branch)
                self.diffuse_iron(afumigatus.next_septa)

    def diffuse_iron_legacy(self, afumigatus=None):
        if afumigatus == None:
            if self.is_root:
                self.diffuse_iron(self)
        else:
            if afumigatus.next_septa is not None and afumigatus.next_branch is None:
                current_iron_pool = afumigatus.iron_pool
                next_iron_pool = afumigatus.next_septa.iron_pool
                iron_pool = (current_iron_pool + next_iron_pool) / 2.0
                afumigatus.iron_pool = iron_pool
                afumigatus.next_septa.iron_pool = iron_pool
                self.diffuse_iron(afumigatus.next_septa)
            elif afumigatus.next_septa is not None and afumigatus.next_branch is not None:
                current_iron_pool = afumigatus.iron_pool
                next_iron_pool = afumigatus.next_septa.iron_pool
                branch_iron_pool = afumigatus.next_branch.iron_pool
                iron_pool = (current_iron_pool + next_iron_pool + branch_iron_pool) / 3.0
                afumigatus.iron__pool = iron_pool
                afumigatus.next_septa.iron_pool = iron_pool
                afumigatus.next_branch.iron_pool = iron_pool
                self.diffuse_iron(afumigatus.next_branch)
                self.diffuse_iron(afumigatus.next_septa)

    def count_cfu(self, afumigatus, cfu=-1):
        if afumigatus.cfu is not None:
            return None
        else:
            if cfu == -1:
                cfu = Id.get_id()
                afumigatus.cfu = cfu
                Afumigatus.total_cfu.add(str(cfu))
                #print(cfu)
                #print(len(Afumigatus.total_cfu))
            else:
                afumigatus.cfu = cfu
                Afumigatus.total_cfu.add(str(cfu))
                #print(len(Afumigatus.total_cfu))

            if afumigatus.next_septa is not None and afumigatus.next_branch is None:
                self.count_cfu(afumigatus.next_septa, cfu)
            elif afumigatus.next_septa is not None and afumigatus.next_branch is not None:
                self.count_cfu(afumigatus.next_branch, cfu)
                self.count_cfu(afumigatus.next_septa, cfu)



    def inc_iron_pool(self, qtty):
        self.iron_pool = self.iron_pool + qtty
        Afumigatus.total_iron = Afumigatus.total_iron + qtty

    def die(self):
        if self.status != Afumigatus.DEAD:
            if self.status == Afumigatus.STERILE_CONIDIA:
                Afumigatus.total_sterile_conidia = Afumigatus.total_sterile_conidia - 1
            self.status = Afumigatus.DEAD
            Afumigatus.total_cells = Afumigatus.total_cells - 1

            #Afumigatus.total_iron = Afumigatus.total_iron - self.iron_pool

    def _lip_activation(self):
        return random() < Util.activation_function(self.iron_pool, Constants.Kd_LIP, h=1, v=Constants.HYPHAE_VOL)

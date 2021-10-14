from edu.uf.interactable.Agent import *

class Macrophage(Phagocyte):

    name = "Macrophage"

    chemokine = None

    InitMacrophageBooleanState = [1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1]

    Dectin1 = 0
    TNFa = 1
    IL6 = 2
    Ft = 3
    DMT1 = 4
    LIP = 5
    TFR = 6
    Fe2 = 7
    IRP = 8
    Hep = 9
    Fpn = 10
    TFBI = 11
    Bglucan = 12

    SPECIES_NUM = 13

    total_cells = 0
    total_iron = 0

    def _init(self, iron_pool):
        Macrophage.total_cells = Macrophage.total_cells + 1
        self.boolean_network = Macrophage.InitMacrophageBooleanState
        self.status = Macrophage.RESTING
        self.state = Macrophage.FREE
        self.fpn = True
        self.fpn_iteration = 0
        self.tf = False
        self.max_move_step = None
        self.tnfa = False
        self.engaged = None ### CHANGE HERE!!!
        Macrophage.total_iron = Macrophage.total_iron + iron_pool

    def _get_max_conidia(self):
        return Constants.MA_MAX_CONIDIA

    def get_max_move_steps(self):
        if self.max_move_step is None:
            if self.status == Macrophage.ACTIVE:
                self.max_move_step = np.random.poisson(Constants.MA_MOVE_RATE_ACT)
            else:
                self.max_move_step = np.random.poisson(Constants.MA_MOVE_RATE_REST)
        return self.max_move_step

    def interact(self, interactable):
        if(type(interactable) is Macrophage):
            return False
        return(interactable.interact(self))

    def process_boolean_network(self):
        pass
        # temp = [0 for _ in range(Macrophage.SPECIES_NUM)]
        #
        # temp[Macrophage.Dectin1] = self.boolean_network[Macrophage.Bglucan]
        # temp[Macrophage.TNFa] = self.boolean_network[Macrophage.Dectin1]
        # temp[Macrophage.IL6] = self.boolean_network[Macrophage.Dectin1]
        # temp[Macrophage.Ft] = self.boolean_network[Macrophage.TNFa] | (-self.boolean_network[Macrophage.IRP] + 1)
        # temp[Macrophage.DMT1] = self.boolean_network[Macrophage.TNFa]
        # temp[Macrophage.LIP] = ((-self.boolean_network[Macrophage.Ft] + 1) | \
        #                         (self.boolean_network[Macrophage.DMT1] & self.boolean_network[Macrophage.Fe2]) | \
        #                         (self.boolean_network[Macrophage.TFBI] & self.boolean_network[Macrophage.TFR]) | (-self.boolean_network[Macrophage.Fpn] + 1))
        # temp[Macrophage.TFR] = self.boolean_network[Macrophage.IRP]
        # temp[Macrophage.Fe2] = 0 #self.hasIron ? 1: 0
        # temp[Macrophage.IRP] = -self.boolean_network[Macrophage.LIP] + 1
        # temp[Macrophage.Hep] = self.boolean_network[Macrophage.IL6]
        # temp[Macrophage.Fpn] = (-self.boolean_network[Macrophage.IRP] + 1) | (-self.boolean_network[Macrophage.Hep] + 1)
        # temp[Macrophage.TFBI] = 0 #self.hasTfbi ? 1: 0
        # temp[Macrophage.Bglucan] = 0 #self.attachedFungus ? 1: 0
        #
        # for i in range(Macrophage.SPECIES_NUM):
        #     self.boolean_network[i] = temp[i]


    def update_status(self):
        from edu.uf.interactable.Afumigatus import Afumigatus
        if self.status == Macrophage.DEAD:
            return
        if self.status == Macrophage.NECROTIC:
            self.die()
            for _, a in self.phagosome.agents.items():
                a.state = Afumigatus.RELEASING
        elif len(self.phagosome.agents) > Constants.MA_MAX_CONIDIA:
            self.status = Phagocyte.NECROTIC
        elif self.status == Macrophage.ACTIVE:
            if self.status_iteration >= Constants.ITER_TO_REST:#30*6:#Constants.ITER_TO_CHANGE_STATE:
                self.status_iteration = 0
                self.tnfa = False
                self.status = Macrophage.RESTING
                #self.state = Macrophage.FREE
            else:
                self.status_iteration = self.status_iteration + 1
        elif self.status == Macrophage.INACTIVE:
            if self.status_iteration >= Constants.ITER_TO_CHANGE_STATE:
                self.status_iteration = 0
                self.status = Macrophage.RESTING
            else:
                pass
                self.status_iteration = self.status_iteration + 1
        elif self.status == Macrophage.ACTIVATING:
            if self.status_iteration >= Constants.ITER_TO_CHANGE_STATE:
                self.status_iteration = 0
                self.status = Macrophage.ACTIVE
            else:
                self.status_iteration = self.status_iteration + 1
        elif self.status == Macrophage.INACTIVATING:
            if self.status_iteration >= Constants.ITER_TO_CHANGE_STATE:
                self.status_iteration = 0
                self.status = Macrophage.INACTIVE
            else:
                self.status_iteration = self.status_iteration + 1
        elif self.status == Macrophage.ANERGIC:
            if self.status_iteration >= Constants.ITER_TO_CHANGE_STATE:
                self.status_iteration = 0
                self.status = Macrophage.RESTING
            else:
                self.status_iteration = self.status_iteration + 1
        #if not (self.status == Macrophage.DEAD or self.status == Macrophage.NECROTIC or self.status == Macrophage.APOPTOTIC):
        #    if random() < Util.activation_function((self.iron_pool - Constants.MA_INTERNAL_IRON), Constants.Kd_MA_IRON, Constants.STD_UNIT_T, Constants.MA_VOL):
        #        self.status = Macrophage.ANERGIC
        #        self.status_iteration = 0

        #self.engaged = False  CHANGE HERE!!!

        if random() < Constants.MA_HALF_LIFE and len(self.phagosome.agents) == 0 and Macrophage.total_cells > Constants.MIN_MA:
            self.die()

        if not self.fpn:
            if self.fpn_iteration >= Constants.ITER_TO_CHANGE_STATE:
                self.fpn_iteration = 0
                self.fpn = True
            else:
                self.fpn_iteration = self.fpn_iteration + 1
        self.move_step = 0
        self.max_move_step = None

    def is_dead(self):
        return self.status == Macrophage.DEAD

    def inc_iron_pool(self, qtty):
        self.iron_pool = self.iron_pool + qtty
        Macrophage.total_iron = Macrophage.total_iron + qtty

    def die(self):
        if self.status != Macrophage.DEAD:
            self.status = Macrophage.DEAD
            Macrophage.total_cells = Macrophage.total_cells - 1

    def attracted_by(self):
        return Macrophage.chemokine.name
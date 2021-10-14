from edu.uf.interactable.ROS import *

class Neutrophil(Phagocyte):

    name = "Neutrophils"


    chemokine = None
    total_cells = 0
    total_iron = 0

    def _init(self, iron_pool):
        Neutrophil.total_cells = Neutrophil.total_cells + 1
        #init boolean network
        self.status = Neutrophil.RESTING
        self.state = Macrophage.FREE
        Neutrophil.total_iron = Neutrophil.total_iron + iron_pool
        self.max_move_step = None
        self.tnfa = False
        self.engaged = False

    def _get_max_conidia(self):
        return Constants.N_MAX_CONIDIA

    def _get_p_phagocytose(self):
        return Constants.PR_N_PHAG

    def get_max_move_steps(self): ##REVIEW
        if self.max_move_step is None:
            if self.status == Macrophage.ACTIVE:
                self.max_move_step = np.random.poisson(Constants.MA_MOVE_RATE_REST)
            else:
                self.max_move_step = np.random.poisson(Constants.MA_MOVE_RATE_REST)
        return self.max_move_step

    def process_boolean_network(self):
        pass

    def update_status(self):
        if self.status == Neutrophil.DEAD:
            return
        if self.status == Neutrophil.NECROTIC or self.status == Neutrophil.APOPTOTIC:
            self.die()
            for _, a in self.phagosome.agents.items():
                a.state = Afumigatus.RELEASING
        elif random() < Constants.NEUTROPHIL_HALF_LIFE:
            self.status = Neutrophil.APOPTOTIC
        elif self.status == Neutrophil.ACTIVE:
            if self.status_iteration > Constants.ITER_TO_CHANGE_STATE:
                self.status_iteration = 0
                self.tnfa = False
                self.status = Neutrophil.RESTING
                self.state = Neutrophil.FREE
            else:
                self.status_iteration = self.status_iteration + 1
        elif self.status == Neutrophil.ACTIVATING:
            if self.status_iteration > Constants.ITER_TO_CHANGE_STATE:
                self.status_iteration = 0
                self.status = Neutrophil.ACTIVE
            else:
                self.status_iteration = self.status_iteration + 1
        self.move_step = 0
        self.max_move_step = None
        self.engaged = False

    def is_dead(self):
        return self.status == Neutrophil.DEAD

    def inc_iron_pool(self, qtty):
        self.iron_pool = self.iron_pool + qtty
        Neutrophil.total_iron = Neutrophil.total_iron + qtty

    def interact(self, interactable):
        itype = type(interactable)
        if itype is Neutrophil:
            return False
        if itype is Afumigatus:
            if self.engaged:
                return True
            if self.status != Neutrophil.APOPTOTIC and self.status != Neutrophil.NECROTIC and self.status != Neutrophil.DEAD:
                if interactable.status == Afumigatus.HYPHAE or interactable.status == Afumigatus.GERM_TUBE:
                    #pr = Constants.PR_N_HYPHAE if interactable.status == Afumigatus.HYPHAE else Constants.PR_N_HYPHAE*0.1
                    pr = Constants.PR_N_HYPHAE
                    if random() < pr:
                        Phagocyte.int_aspergillus(self, interactable)
                        interactable.status = Afumigatus.DYING
                    else:
                        self.engaged = True
                elif interactable.status == Afumigatus.SWELLING_CONIDIA:
                    if random() < Constants.PR_N_PHAG:
                        Phagocyte.int_aspergillus(self, interactable)
                    else:
                        pass
                        #interactable.status = Afumigatus.STERILE_CONIDIA
                        #Afumigatus.total_sterile_conidia = Afumigatus.total_sterile_conidia + 1
            return True
        if itype is Macrophage:
            if self.status == Neutrophil.APOPTOTIC and len(interactable.phagosome.agents) == 0:
                interactable.phagosome.agents = self.phagosome.agents
                interactable.inc_iron_pool(self.iron_pool)
                self.inc_iron_pool(self.iron_pool)
                self.die()
                interactable.status = Macrophage.INACTIVE
            return True
        if itype is Transferrin:
            return False
        if itype is TAFC:
            return False
        if itype is Iron:
            if self.status == Neutrophil.NECROTIC:# or self.status == Neutrophil.APOPTOTIC or self.status == Neutrophil.DEAD:
                interactable.inc(self.iron_pool)
                self.inc_iron_pool(-self.iron_pool)
            return False
        if itype is ROS and interactable.state == Neutrophil.INTERACTING:
            if self.status == Neutrophil.ACTIVE:
                interactable.inc(0)
            return True
        return interactable.interact(self)

    def die(self):
        if self.status != Neutrophil.DEAD:
            self.status = Neutrophil.DEAD
            Neutrophil.total_cells = Neutrophil.total_cells - 1

    def attracted_by(self):
        return Neutrophil.chemokine.name
from cc3d.core.PySteppables import *
import numpy as np
import os
import Parameters

IFNWash = False #Whether the plate is prestimulated with IFNe before infection

min_to_mcs = 10.0  # min/mcs
hours_to_mcs = min_to_mcs / 60.0 # hours/mcs
days_to_mcs = min_to_mcs / 1440.0  # day/mcs
hours_to_simulate = 80.0  # 10 in the original model

virus_diffusion_coefficient = 1.0/10.0 #vl^2 / min
IFNe_diffusion_coefficient = 1.0/10.0 #vl^2 / min

Replicate = Parameters.R

folder_path = '/Users/joaponte/Desktop/IFNModel/'
if not os.path.exists(folder_path):
    os.makedirs(folder_path)

#Cell Transition Model
FluModel_string = '''        
        model FluModel()

        //State Variables and Transitions
        V1: -> T  ; -beta * V * T ;                             // Susceptible Cells
        V2: -> I1 ;  beta * V * T - k * I1 ;                    // Early Infected Cells
        V3: -> I2 ;  k * I1 - delta_d * I2 / (K_delta + I2) ;   // Late Infected Cells
        V4: -> V  ;  p * I2 - c * V ;                           // Extracellular Virus
        V5: -> D  ;  delta_d * I2 / (K_delta + I2) ;            // Cleared Infected Cells (for Bookkeeping)

        //Parameters
        beta = 2.4* 10^(-4) ;                                   // Virus Infective
        p = 1.6 ;                                               // Virus Production
        c = 13.0 ;                                              // Virus Clearance
        k = 4.0 ;                                               // Eclipse phase
        delta_d = 1.6 * 10^6 ;                                  // Infected Cell Clearance
        K_delta = 4.5 * 10^5 ;                                  // Half Saturation Constant         

        // Initial Conditions ;
        T0 = 1.0*10^7;
        T = T0  ;                                               // Initial Number of Uninfected Cells
        I1 = 75.0 ;                                             // Initial Number of Infected Cells
end'''

# Viral Replication Model
viral_model_string = '''
    E7a: H ->           ; H*k61*V                     ;
    E8a: -> V           ; H*k71*V/(1.0+k72*IFNe*7E-5) ;
    E8b: V ->           ; k73*V                       ;

    //Parameters
    k61 = 0.635     ;
    k71 = 1.537     ;
    k72 = 47.883    ;
    k73 = 0.197     ;

    //Initial Conditions
    V =  0.0      ; 
    H = 1.0          ;
    
    //Inputs
    IFNe  =  0.0     ;
'''

IFN_model_string = '''
    //Equations
    E2a: -> IFN         ; H*(k11*RIGI*V+k12*(V^n)/(k13+(V^n))+k14*IRF7P)    ;
    E2b: IFN ->         ; k21*IFN                                           ;
    E4a: -> STATP       ; H*k31*IFNe/(k32+k33*IFNe)                         ;
    E4b: STATP ->       ; t3*STATP                                          ;
    E5a: -> IRF7        ; H*(k41*STATP+k42*IRF7P)                           ;
    E5b: IRF7 ->        ; t4*IRF7                                           ;
    E6a: -> IRF7P       ; H*k51*IRF7                                        ;
    E6b: IRF7P ->       ; t5*IRF7P                                          ;

    //Parameters
    //k11 = 10.0^(5)  ; 
    k11 = 0.0       ; 
    k12 = 9.746     ; 
    k13 = 12.511    ; 
    k14 = 13.562    ;
    k21 = 10.385    ;
    k31 = 45.922    ;
    k32 = 5.464     ;
    k33 = 0.068     ;
    t3  = 0.3       ;
    k41 = 0.115     ;
    k42 = 1.053     ;
    t4  = 0.75      ;
    k51 = 0.202     ;
    t5  = 0.3       ;
    n   = 3.0       ;
    RIGI = 1.0      ;
    
    // Inputs
    H    = 0.0      ;
    IFNe = 0.0      ;
    V = 0.0         ;
'''

## Global Parameters
# IFN Decay Rate
t2 = 3.481

class ODEModelSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def start(self):
        # Store Initial Number of Cells
        self.shared_steppable_vars['InitialNumberCells'] = len(self.cell_list)

        # Set Max Simulation Steps
        self.get_xml_element('simulation_steps').cdata = hours_to_simulate / hours_to_mcs

        # Load Original FLU ODE Model
        self.add_free_floating_antimony(model_string=FluModel_string, model_name='FluModel',
                                        step_size=days_to_mcs)

        # Load Viral Model inside Cells
        self.add_antimony_to_cell_types(model_string=viral_model_string, model_name='VModel',
                                        cell_types=[self.U], step_size=hours_to_mcs)

        # Load IFN Model inside Cells
        self.add_antimony_to_cell_types(model_string=IFN_model_string, model_name='IModel',
                                        cell_types=[self.U], step_size=hours_to_mcs)

        # Parameter Scan
        # self.sbml.FluModel['beta'] *= Multiplier
        # for cell in self.cell_list:
            # cell.sbml.IModel['k11'] *= Multiplier
            # cell.sbml.VModel['k73'] *= Multiplier

        # Initial conditions: infected cell in the center
        cell = self.cell_field[self.dim.x // 2, self.dim.y // 2, 0]
        cell.type = self.I1
        cell.sbml.VModel['V'] = 6.9e-8
        self.sbml.FluModel['I1'] = 1.0 / self.shared_steppable_vars['InitialNumberCells']
        self.sbml.FluModel['V'] = 0.0
        
        #Set prestimulated internal protein values
        if IFNWash:
            for cell in self.cell_list_by_type(self.U,self.I1):
                cell.sbml.IModel['IFN'] = 0.035
                cell.sbml.IModel['IRF7'] = 0.097
                cell.sbml.IModel['IRF7P'] = 0.028
                cell.sbml.IModel['STATP'] = 0.714   

class CellularModelSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def start(self):
        # Set IFNe diffusion parameters
        self.get_xml_element('IFNe_dc').cdata = IFNe_diffusion_coefficient * min_to_mcs
        self.get_xml_element('IFNe_decay').cdata = t2 * hours_to_mcs

        # Set Virus diffusion parameters
        self.get_xml_element('virus_dc').cdata = virus_diffusion_coefficient * min_to_mcs
        self.get_xml_element('virus_decay').cdata = self.sbml.FluModel['c'] * days_to_mcs

        # Set secretors
        self.secretorIFN = self.get_field_secretor("IFNe")
        self.secretorV = self.get_field_secretor("Virus")

    def step(self, mcs):
        ## Measure amount of IFNe in the Field
        self.shared_steppable_vars['ExtracellularIFN_Field'] = 0
        for cell in self.cell_list_by_type(self.U,self.I1,self.I2):
            self.shared_steppable_vars['ExtracellularIFN_Field'] += self.secretorIFN.amountSeenByCell(cell)

        ## Production of IFNe
        # E2b: IFN -> IFNe; k21 * IFN ;
        for cell in self.cell_list_by_type(self.U,self.I1,self.I2):
            intracellularIFN = cell.sbml.IModel['IFN']
            k21 = cell.sbml.IModel['k21'] * hours_to_mcs
            p = k21 * intracellularIFN
            self.secretorIFN.secreteInsideCellTotalCount(cell, p / cell.volume)

        ## Measure amount of extracellular virus field
        self.shared_steppable_vars['ExtracellularVirus_Field'] = 0
        for cell in self.cell_list:
            V = self.secretorV.amountSeenByCell(cell)
            self.shared_steppable_vars['ExtracellularVirus_Field'] += V

        ## Production of extracellular virus
        # E8b: V -> ; k73 * V
        for cell in self.cell_list_by_type(self.I2):
            k73 = cell.sbml.VModel['k73'] * hours_to_mcs
            Virus = cell.sbml.VModel['V']
            p = k73 * Virus * 1094460.28
            self.secretorV.secreteInsideCellTotalCount(cell, p / cell.volume)

        ## P to D transition
        # E7a: P -> ; P * k61 * V;
        for cell in self.cell_list_by_type(self.I2):
            k61 = cell.sbml.VModel['k61'] * hours_to_mcs
            H = cell.sbml.VModel['H']
            V = cell.sbml.VModel['V']
            r = k61 * V * (1-H)
            p_I2toD = 1.0 - np.exp(-r)
            if np.random.random() < p_I2toD:
                cell.type = self.DEAD

        ## I1 to I2 transition
        # E2: I1 -> I2 ; k * I1
        for cell in self.cell_list_by_type(self.I1):
            k = self.sbml.FluModel['k'] * days_to_mcs
            r = k
            p_T1oI2 = 1.0 - np.exp(-r)
            if np.random.random() < p_T1oI2:
                cell.type = self.I2

        ## U to I1 transition
        # E1: T -> I1 ; beta * V * T
        for cell in self.cell_list_by_type(self.U):
            b = self.sbml.FluModel['beta'] * self.shared_steppable_vars['InitialNumberCells'] * days_to_mcs
            V = self.secretorV.amountSeenByCell(cell)
            r = b * V
            p_UtoI1 = 1.0 - np.exp(-r)
            if np.random.random() < p_UtoI1:
                cell.type = self.I1
                cell.sbml.VModel['V'] = 6.9e-8

        ## Updating Cellular Models
        for cell in self.cell_list:
            cell.sbml.VModel['IFNe'] = self.secretorIFN.amountSeenByCell(cell)
            cell.sbml.IModel['IFNe'] = self.secretorIFN.amountSeenByCell(cell)
            cell.sbml.IModel['H'] = cell.sbml.VModel['H']
            cell.sbml.IModel['V'] = cell.sbml.VModel['V']

        self.timestep_sbml()

class OutputSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def start(self):
        # Output Cellular Data
        file_name1 = 'FullModelCellular_%i.txt' % (Replicate)
        self.output1 = open(folder_path + file_name1, 'w')
        self.output1.write("%s,%s,%s,%s,%s,%s,%s\n" %('Time','U','I1','I2','D','Ve','IFNe'))
        self.output1.flush()

        # Output Intracellular Data
        file_name2 = 'FullModelIntracellular_%i.txt' % (Replicate)
        self.output2 = open(folder_path + file_name2, 'w')
        self.output2.write("%s,%s,%s,%s,%s,%s,%s,%s,%s\n" %
                          ('Time','V','H','P','IFNe','STATP','IRF7','IRF7P','IFN'))
        self.output2.flush()

        #IFNe secretor
        self.secretorIFNe = self.get_field_secretor("IFNe")
            
    def step(self, mcs):
        Time = mcs * hours_to_mcs
        U = len(self.cell_list_by_type(self.U)) / self.shared_steppable_vars['InitialNumberCells']
        I1 =  len(self.cell_list_by_type(self.I1)) / self.shared_steppable_vars['InitialNumberCells']
        I2 = len(self.cell_list_by_type(self.I2)) / self.shared_steppable_vars['InitialNumberCells']
        D =  len(self.cell_list_by_type(self.DEAD)) / self.shared_steppable_vars['InitialNumberCells']
        Ve = self.shared_steppable_vars['ExtracellularVirus_Field']
        IFNe = 0.0
        for cell in self.cell_list:
            IFNe += self.secretorIFNe.amountSeenByCell(cell)

        self.output1.write("%e,%e,%e,%e,%e,%e,%e\n" % (Time,U,I1,I2,D,Ve,IFNe))
        self.output1.flush()

        L = len(self.cell_list_by_type(self.U,self.I1,self.I2))
        P = L / self.shared_steppable_vars['InitialNumberCells']
        V = 0.0
        H = 0.0
        STATP = 0.0
        IRF7 = 0.0
        IRF7P = 0.0
        IFN = 0.0
        for cell in self.cell_list_by_type(self.U,self.I1,self.I2):
            V += cell.sbml.VModel['V'] / L
            H += cell.sbml.VModel['H'] / L
            STATP += cell.sbml.IModel['STATP'] / L
            IRF7 += cell.sbml.IModel['IRF7'] / L
            IRF7P += cell.sbml.IModel['IRF7P'] / L
            IFN += cell.sbml.IModel['IFN'] / L
        IFNe = self.shared_steppable_vars['ExtracellularIFN_Field'] \
                          / self.shared_steppable_vars['InitialNumberCells']
        self.output2.write("%e,%e,%e,%e,%e,%e,%e,%e,%e\n" %
                           (Time,V,H,P,IFNe,STATP,IRF7,IRF7P,IFN))
        self.output2.flush()

class PlaqueAssaySteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def start(self):
        file_name3 = 'PlaqueAssay_%i.txt' % (Replicate)
        self.output3 = open(folder_path + file_name3, 'w')
        self.output3.write("%s,%s,%s,%s\n" % ('Time', 'avgI1rd', 'avgI2rd', 'avgDrd'))
        self.output3.flush()

    def step(self, mcs):
        # Measure area occupied by D cells and assume its a circle
        volume_D = 0.0
        for cell in self.cell_list_by_type(self.DEAD):
            volume_D += cell.volume
        avgDrd = np.sqrt(volume_D/np.pi)

        volume_I2 = 0.0
        for cell in self.cell_list_by_type(self.I2):
            volume_I2 += cell.volume
        avgI2rd = np.sqrt((volume_D+volume_I2)/np.pi)

        volume_I1 = 0.0
        for cell in self.cell_list_by_type(self.I1):
            volume_I1 += cell.volume
        avgI1rd = np.sqrt((volume_D+volume_I2+volume_I1)/np.pi)

        Time = mcs * hours_to_mcs
        self.output3.write("%e,%e,%e,%e\n" % (Time, avgI1rd, avgI2rd, avgDrd))
        self.output3.flush()
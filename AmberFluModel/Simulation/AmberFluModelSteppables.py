from cc3d.core.PySteppables import *
import numpy as np
import os

plot_ODEModel = True
plot_CellModel = True
plot_Residuals = False
Data_writeout = False

## How to determine V
# 1 pulls from the scalar virus from the ODE original model (no feedback in the cellular model)
# 2 pulls from the scalar virus from the cellular model (feedback in the cellular model but no field)
# 3 pulls from the virus field
how_to_determine_V = 3

min_to_mcs = 10.0  # min/mcs
days_to_mcs = min_to_mcs / 1440.0  # day/mcs
days_to_simulate = 4.0 #10 in the original model

'''Smith AP, Moquin DJ, Bernhauerova V, Smith AM. Influenza virus infection model with density dependence 
supports biphasic viral decay. Frontiers in microbiology. 2018 Jul 10;9:1554.'''

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

class AmberFluModelSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def start(self):
        # Uptading max simulation steps using scaling factor to simulate 10 days
        self.get_xml_element('simulation_steps').cdata = days_to_simulate / days_to_mcs

        # Adding free floating antimony model
        self.add_free_floating_antimony(model_string=FluModel_string, model_name='FluModel',
                                        step_size=days_to_mcs)

        # Changing initial values according to discussions with Amber Smith
        self.sbml.FluModel['I1'] = 0.0
        self.sbml.FluModel['V'] = 75.0

    def step(self, mcs):
        self.timestep_sbml()

class CellularModelSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def start(self):
        # set initial model parameters
        self.initial_uninfected = len(self.cell_list_by_type(self.U))
        self.ExtracellularVirus = self.sbml.FluModel['V']
        self.get_xml_element('virus_decay').cdata = self.sbml.FluModel['c'] * days_to_mcs

    def step(self, mcs):
        # Transition rule from U to I1
        secretor = self.get_field_secretor("Virus")
        for cell in self.cell_list_by_type(self.U):
            # Determine V from scalar virus from the ODE
            if how_to_determine_V == 1:
                b = self.sbml.FluModel['beta'] * self.sbml.FluModel['T0'] * days_to_mcs
                V = self.sbml.FluModel['V'] / self.sbml.FluModel['T0']

            # Determine V from scalar virus from the cellular model
            if how_to_determine_V == 2:
                b = self.sbml.FluModel['beta'] * self.initial_uninfected * days_to_mcs
                V = self.ExtracellularVirus / self.initial_uninfected

            # Determine V from the virus field
            if how_to_determine_V == 3:
                b = self.sbml.FluModel['beta'] * self.initial_uninfected * days_to_mcs
                uptake_probability = 0.0000001
                uptake = secretor.uptakeInsideCellTotalCount(cell, 1E6, uptake_probability)
                V = abs(uptake.tot_amount) / uptake_probability
                secretor.secreteInsideCellTotalCount(cell, abs(uptake.tot_amount) / cell.volume)

            # Calculate the probability of infection of individual cells based on the amount of virus PER cell
            p_UtoI1 = b * V
            if np.random.random() < p_UtoI1:
                cell.type = self.I1

        # Transition rule from I1 to I2
        k = self.sbml.FluModel['k'] * days_to_mcs
        p_T1oI2 = k
        for cell in self.cell_list_by_type(self.I1):
            if np.random.random() < p_T1oI2:
                cell.type = self.I2

        # Transition rule from I2 to D
        K_delta = self.sbml.FluModel['K_delta'] / self.sbml.FluModel['T0'] * self.initial_uninfected
        delta_d = self.sbml.FluModel['delta_d'] / self.sbml.FluModel['T0'] * self.initial_uninfected
        I2 = len(self.cell_list_by_type(self.I2))
        p_T2toD = delta_d / (K_delta + I2) * days_to_mcs
        for cell in self.cell_list_by_type(self.I2):
            if np.random.random() < p_T2toD:
                cell.type = self.DEAD

        # Production of extracellular virus
        secretor = self.get_field_secretor("Virus")
        V = self.ExtracellularVirus
        p = self.sbml.FluModel['p'] / self.initial_uninfected * self.sbml.FluModel['T0'] * days_to_mcs
        c = self.sbml.FluModel['c'] * days_to_mcs
        for cell in self.cell_list_by_type(self.I2):
            release = secretor.secreteInsideCellTotalCount(cell, p / cell.volume)
            self.ExtracellularVirus += release.tot_amount
        self.ExtracellularVirus -= c * V

        # Measure amount of extracellular virus field
        self.ExtracellularVirus_Field = 0
        for cell in self.cell_list:
            uptake_probability = 0.0000001
            uptake = secretor.uptakeInsideCellTotalCount(cell, 1E6, uptake_probability)
            V = abs(uptake.tot_amount) / uptake_probability
            self.ExtracellularVirus_Field += V
            secretor.secreteInsideCellTotalCount(cell, abs(uptake.tot_amount) / cell.volume)

        # Dictonary to pass information between steppables
        self.shared_steppable_vars['ExtracellularVirus'] = self.ExtracellularVirus
        self.shared_steppable_vars['ExtracellularVirus_Field'] = self.ExtracellularVirus_Field

class PlotSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def start(self):
        self.initial_uninfected = len(self.cell_list_by_type(self.U))
        if (plot_ODEModel == True) or (plot_CellModel == True):
            self.plot_win = self.add_new_plot_window(title='Flu Model Cells',
                                                     x_axis_title='Days',
                                                     y_axis_title='Variables', x_scale_type='linear',
                                                     y_scale_type='linear',
                                                     grid=False, config_options={'legend': True})

            self.plot_win2 = self.add_new_plot_window(title='Flu Model Virus',
                                                      x_axis_title='Days',
                                                      y_axis_title='Virus', x_scale_type='linear',
                                                      y_scale_type='linear',
                                                      grid=False, config_options={'legend': True})

            if plot_ODEModel == True:
                self.plot_win.add_plot("ODET", style='Dots', color='blue', size=5)
                self.plot_win.add_plot("ODEI1", style='Dots', color='orange', size=5)
                self.plot_win.add_plot("ODEI2", style='Dots', color='red', size=5)
                self.plot_win.add_plot("ODED", style='Dots', color='purple', size=5)
                self.plot_win2.add_plot("ODEV", style='Dots', color='blue', size=5)

            if plot_CellModel == True:
                self.plot_win.add_plot("CC3DT", style='Lines', color='blue', size=5)
                self.plot_win.add_plot("CC3DI1", style='Lines', color='orange', size=5)
                self.plot_win.add_plot("CC3DI2", style='Lines', color='red', size=5)
                self.plot_win.add_plot("CC3DD", style='Lines', color='purple', size=5)
                self.plot_win2.add_plot("CC3DV", style='Lines', color='blue', size=5)

    def step(self, mcs):
        if (plot_ODEModel == True) or (plot_CellModel == True):
            if plot_ODEModel == True:
                self.plot_win.add_data_point("ODET", mcs * days_to_mcs,
                                             self.sbml.FluModel['T'] / self.sbml.FluModel['T0'])
                self.plot_win.add_data_point("ODEI1", mcs * days_to_mcs,
                                             self.sbml.FluModel['I1'] / self.sbml.FluModel['T0'])
                self.plot_win.add_data_point("ODEI2", mcs * days_to_mcs,
                                             self.sbml.FluModel['I2'] / self.sbml.FluModel['T0'])
                self.plot_win.add_data_point("ODED", mcs * days_to_mcs,
                                             self.sbml.FluModel['D'] / self.sbml.FluModel['T0'])
                self.plot_win2.add_data_point("ODEV", mcs * days_to_mcs, np.log10(self.sbml.FluModel['V']))

            if plot_CellModel == True:
                self.plot_win.add_data_point("CC3DT", mcs * days_to_mcs,
                                              len(self.cell_list_by_type(self.U)) / self.initial_uninfected)
                self.plot_win.add_data_point("CC3DI1", mcs * days_to_mcs,
                                              len(self.cell_list_by_type(self.I1)) / self.initial_uninfected)
                self.plot_win.add_data_point("CC3DI2", mcs * days_to_mcs,
                                              len(self.cell_list_by_type(self.I2)) / self.initial_uninfected)
                self.plot_win.add_data_point("CC3DD", mcs * days_to_mcs,
                                              len(self.cell_list_by_type(self.DEAD)) / self.initial_uninfected)
                if how_to_determine_V == 3:
                    self.plot_win2.add_data_point("CC3DV", mcs * days_to_mcs,
                                              np.log10(self.shared_steppable_vars['ExtracellularVirus_Field']))
                else:
                    self.plot_win2.add_data_point("CC3DV", mcs * days_to_mcs,
                                              np.log10(self.shared_steppable_vars['ExtracellularVirus']))

class StatisticsSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def start(self):
        self.initial_uninfected = len(self.cell_list)

        self.cellular_infection = False
        self.cellular_infection_time = 0.0
        self.Ambersmodel_infection = False
        self.Ambersmodel_infection_time = 0.0
        self.infection_threshold = 0.1

        if plot_Residuals:
            self.plot_win5 = self.add_new_plot_window(title='Residuals',
                                                      x_axis_title='days',
                                                      y_axis_title='Variables', x_scale_type='linear',
                                                      y_scale_type='linear',
                                                      grid=False, config_options={'legend': True})
            self.plot_win5.add_plot("dU", style='Lines', color='blue', size=5)
            self.plot_win5.add_plot("dI1", style='Lines', color='yellow', size=5)
            self.plot_win5.add_plot("dI2", style='Lines', color='red', size=5)
            self.plot_win5.add_plot("dD", style='Lines', color='purple', size=5)

    def step(self, mcs):
        if self.cellular_infection == False:
            if len(self.cell_list_by_type(self.I1)) / self.initial_uninfected >= self.infection_threshold:
                self.cellular_infection_time = mcs
                self.cellular_infection = True

        if self.Ambersmodel_infection == False:
            if self.sbml.FluModel['I1'] / self.sbml.FluModel['T0'] >= self.infection_threshold:
                self.Ambersmodel_infection_time = mcs
                self.Ambersmodel_infection = True

        # print("Cellular Infection = ", self.cellular_infection_time * days_to_mcs)
        # print("ODE Infection = ", self.Ambersmodel_infection_time * days_to_mcs)

        dU = (len(self.cell_list_by_type(self.U)) / self.initial_uninfected) - (
                self.sbml.FluModel['T'] / self.sbml.FluModel['T0'])
        dI1 = (len(self.cell_list_by_type(self.I1)) / self.initial_uninfected) - (
                self.sbml.FluModel['I1'] / self.sbml.FluModel['T0'])
        dI2 = (len(self.cell_list_by_type(self.I2)) / self.initial_uninfected) - (
                self.sbml.FluModel['I2'] / self.sbml.FluModel['T0'])
        dD = (len(self.cell_list_by_type(self.DEAD)) / self.initial_uninfected) - (
                self.sbml.FluModel['D'] / self.sbml.FluModel['T0'])

        if plot_Residuals:
            self.plot_win5.add_data_point("dU", mcs * days_to_mcs, dU)
            self.plot_win5.add_data_point("dI1", mcs * days_to_mcs, dI1)
            self.plot_win5.add_data_point("dI2", mcs * days_to_mcs, dI2)
            self.plot_win5.add_data_point("dD", mcs * days_to_mcs, dD)

class Data_OutputSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def start(self):
        if Data_writeout:
            folder_path = '/Users/Josua/AmberFluModelData/'
            if not os.path.exists(folder_path):
                os.makedirs(folder_path)

            file_name = 'ODEModel.txt'
            self.output1 = open(folder_path + file_name, 'w')
            self.output1.write(
                "%s,%s,%s,%s,%s,%s\n" % ('Time','T', 'I1', 'I2', 'D', 'V'))
            self.output1.flush()

            file_name = 'CC3DModel.txt'
            self.output2 = open(folder_path + file_name, 'w')
            self.output2.write(
                "%s,%s,%s,%s,%s,%s\n" % ('Time', 'T', 'I1', 'I2', 'D', 'V'))
            self.output2.flush()

    def step(self, mcs):
        if Data_writeout:
            # Record variables from ODE model
            d = mcs * days_to_mcs
            AT = self.sbml.FluModel['T']
            AI1 = self.sbml.FluModel['I1']
            AI2 = self.sbml.FluModel['I2']
            AD = self.sbml.FluModel['D']
            AV = self.sbml.FluModel['V']

            self.output1.write("%.2f,%.2f,%.2f,%.2f,%.2f,%.2f\n" % (d, AT, AI1, AI2, AD, AV))
            self.output1.flush()

            # Record variables from Cellularized Model
            U = len(self.cell_list_by_type(self.U))
            I1 = len(self.cell_list_by_type(self.I1))
            I2 = len(self.cell_list_by_type(self.I2))
            D = len(self.cell_list_by_type(self.DEAD))
            if how_to_determine_V == 3:
                V = self.shared_steppable_vars['ExtracellularVirus_Field']
            else:
                V = self.shared_steppable_vars['ExtracellularVirus']

            self.output2.write("%.2f,%.2f,%.2f,%.2f,%.2f,%.2f\n" % (d, U, I1, I2, D, V))
            self.output2.flush()

    def finish(self):
        if Data_writeout:
            self.output1.close()
            self.output2.close()

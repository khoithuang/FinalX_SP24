import sys
import numpy as np
from PyQt5 import QtWidgets as qtw
from matplotlib import pyplot as plt
from Air import air, units, StateDataForPlotting  # Assuming these are properly defined in a separate module


class DualCycleModel:
    """
    Represents an air standard dual cycle, consisting of five main thermodynamic processes:
    1. Isentropic compression (1 to 2)
    2. Constant volume heat addition (2 to 3)
    3. Constant pressure heat addition (3 to 4)
    4. Isentropic expansion (4 to 5)
    5. Constant volume heat rejection (5 to 1)
    """

    def __init__(self, p_initial=1E5, t_initial=300, v_cylinder=0.1, r=18, rc=1.2, p_ratio=1.5):
        """
        Initializes the dual cycle model with given parameters.
        :param p_initial: Initial pressure in Pa.
        :param t_initial: Initial temperature in K.
        :param v_cylinder: Total cylinder volume in m^3.
        :param r: Compression ratio (V1/V2).
        :param rc: Cutoff ratio (V4/V3).
        :param p_ratio: Pressure ratio (P3/P2).
        """
        self.air = air()
        self.units = units()
        self.units.SI = True  # Using SI units

        self.r = r
        self.rc = rc
        self.p_ratio = p_ratio
        self.p_initial = p_initial
        self.t_initial = t_initial
        self.v_cylinder = v_cylinder

        self.calculate_cycle()

    def calculate_cycle(self):
        """
        Calculates all states and processes for the dual cycle based on initial conditions and cycle parameters.
        """
        # State 1: Initial state
        self.State1 = self.air.set(P=self.p_initial, T=self.t_initial)
        v1 = self.v_cylinder

        # State 2: After isentropic compression
        v2 = v1 / self.r
        self.State2 = self.air.set(v=v2, s=self.State1.s)

        # State 3: Constant volume heat addition
        p3 = self.State2.P * self.p_ratio
        self.State3 = self.air.set(P=p3, v=v2)

        # State 4: Constant pressure heat addition
        v4 = v2 * self.rc
        self.State4 = self.air.set(P=p3, v=v4)

        # State 5: Isentropic expansion to initial volume
        self.State5 = self.air.set(v=v1, s=self.State4.s)

        # Update thermodynamic properties for the cycle
        self.update_thermodynamic_properties()

    def update_thermodynamic_properties(self):
        """
        Updates the work and heat transfer calculations for the cycle based on state properties.
        """
        # Work done during compression and expansion
        self.work_compression = self.State2.u - self.State1.u
        self.work_expansion = self.State5.u - self.State4.u

        # Heat added and rejected
        self.heat_added_cv = self.State3.u - self.State2.u
        self.heat_added_cp = self.State4.h - self.State3.h
        self.heat_rejected = self.State5.u - self.State1.u

        # Thermal efficiency
        self.efficiency = 1 - (self.heat_rejected / (self.heat_added_cv + self.heat_added_cp))


class DualCycleController:
    """
    Controls interactions between the dual cycle model and the view.
    """

    def __init__(self, model, view):
        self.model = model
        self.view = view
        self.view.set_controller(self)

    def update_parameters(self, p_initial, t_initial, v_cylinder, r, rc, p_ratio):
        """
        Updates the cycle parameters in the model and recalculates the cycle.
        """
        self.model.p_initial = p_initial
        self.model.t_initial = t_initial
        self.model.v_cylinder = v_cylinder
        self.model.r = r
        self.model.rc = rc
        self.model.p_ratio = p_ratio
        self.model.calculate_cycle()
        self.view.update_view()


class DualCycleView(qtw.QWidget):
    """
    Provides the GUI components for the dual cycle simulation.
    """

    def __init__(self):
        super().__init__()
        self.init_ui()

    def set_controller(self, controller):
        self.controller = controller

    def init_ui(self):
        # Create layout and widgets
        self.layout = qtw.QVBoxLayout()
        self.parameter_inputs = qtw.QGroupBox("Cycle Parameters")
        self.plt_layout = qtw.QHBoxLayout()
        self.figure, self.ax = plt.subplots()
        self.canvas = FigureCanvasQTAgg(self.figure)

        # Set up parameter inputs
        self.p_initial_input = qtw.QLineEdit()
        self.t_initial_input = qtw.QLineEdit()
        self.v_cylinder_input = qtw.QLineEdit()
        self.r_input = qtw.QLineEdit()
        self.rc_input = qtw.QLineEdit()
        self.p_ratio_input = qtw.QLineEdit()
        self.update_button = qtw.QPushButton("Update Cycle")
        self.update_button.clicked.connect(self.on_update_clicked)

        # Assemble layout
        param_layout = qtw.QFormLayout()
        param_layout.addRow("Initial Pressure (Pa)", self.p_initial_input)
        param_layout.addRow("Initial Temperature (K)", self.t_initial_input)
        param_layout.addRow("Cylinder Volume (m^3)", self.v_cylinder_input)
        param_layout.addRow("Compression Ratio (r)", self.r_input)
        param_layout.addRow("Cutoff Ratio (rc)", self.rc_input)
        param_layout.addRow("Pressure Ratio (P3/P2)", self.p_ratio_input)
        self.parameter_inputs.setLayout(param_layout)

        self.layout.addWidget(self.parameter_inputs)
        self.layout.addWidget(self.update_button)
        self.layout.addWidget(self.canvas)
        self.setLayout(self.layout)

    def on_update_clicked(self):
        # Read input values and update model through the controller
        p_initial = float(self.p_initial_input.text())
        t_initial = float(self.t_initial_input.text())
        v_cylinder = float(self.v_cylinder_input.text())
        r = float(self.r_input.text())
        rc = float(self.rc_input.text())
        p_ratio = float(self.p_ratio_input.text())

        self.controller.update_parameters(p_initial, t_initial, v_cylinder, r, rc, p_ratio)

    def update_view(self):
        # Update the plot with new cycle data
        self.ax.clear()
        self.plot_cycle()
        self.canvas.draw()

    def plot_cycle(self):
        # Plotting the PV diagram or other thermodynamic properties
        states = [self.controller.model.State1, self.controller.model.State2, self.controller.model.State3,
                  self.controller.model.State4, self.controller.model.State5]
        p = [state.P for state in states]
        v = [state.v for state in states]
        self.ax.plot(v, p, 'o-')
        self.ax.set_xlabel("Volume (m^3)")
        self.ax.set_ylabel("Pressure (Pa)")
        self.ax.set_title("Pressure-Volume Diagram for Dual Cycle")


def main():
    app = qtw.QApplication(sys.argv)
    model = DualCycleModel()
    view = DualCycleView()
    controller = DualCycleController(model, view)
    view.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()

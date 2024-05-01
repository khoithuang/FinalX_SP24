# I got help from ChatGPT
from PyQt5 import QtWidgets
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import numpy as np
from mpmath import quad


class TakeoffView(QtWidgets.QWidget):
    """
   The GUI view component of the application that displays input fields for weight and thrust,
   a button to trigger calculation, and a canvas to plot the graph.
   """
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Airplane Take-off Distance Calculator")
        self.layout = QtWidgets.QVBoxLayout(self)

        # Input fields
        self.weight_input = QtWidgets.QLineEdit(self)
        self.thrust_input = QtWidgets.QLineEdit(self)
        self.calculate_button = QtWidgets.QPushButton("Calculate", self)

        self.layout.addWidget(self.weight_input)
        self.layout.addWidget(self.thrust_input)
        self.layout.addWidget(self.calculate_button)

        # Matplotlib Graph
        self.figure = Figure()
        self.canvas = FigureCanvas(self.figure)
        self.layout.addWidget(self.canvas)

    def plot_graph(self, thrust, weights, stos):
        """
        Plots the take-off distance for different weights at a given thrust.

        Parameters:
        thrust (list): List of thrust values for plotting.
        weights (list): List of aircraft weights in lbs.
        stos (list): Calculated STO distances corresponding to the weights.
        """
        ax = self.figure.add_subplot(111)
        ax.clear()
        for weight, sto in zip(weights, stos):
            ax.plot(thrust, sto, marker='o', label=f'Weight: {weight} lbs')
        ax.set_xlabel('Thrust (lbs)')
        ax.set_ylabel('Take-off Distance (ft)')
        ax.legend()
        self.canvas.draw()


class TakeoffController:
    """
    The controller component that handles interaction between the view and the model.
    """
    def __init__(self, view, model):
        self.view = view
        self.model = model
        self.view.calculate_button.clicked.connect(self.on_calculate_clicked)

    def on_calculate_clicked(self):
        """
        Handles the button click event to compute and plot STO based on input values.
        """
        weight = float(self.view.weight_input.text())
        thrust = float(self.view.thrust_input.text())

        weights = [weight - 10000, weight, weight + 10000]
        stos = [self.model.calculate_sto(w, thrust) for w in weights]

        thrust_values = [thrust] * len(weights)
        self.view.plot_graph(thrust_values, weights, stos)


class TakeoffModel:
    """
    The model component that performs calculations for take-off distances.
    """
    def __init__(self):
        self.gc = 32.174  # gravitational constant in lbm*ft/(lbf*s^2)

    def calculate_sto(thrust, weight, cl_max, cd, rho, s, gc):
        """
        Calculates the short take-off distance using given parameters.

        Parameters:
        thrust (float): Engine thrust in lbs.
        weight (float): Aircraft weight in lbs.
        cl_max (float): Maximum lift coefficient.
        cd (float): Drag coefficient.
        rho (float): Air density in slugs/ft^3.
        s (float): Wing surface area in ft^2.
        gc (float): Gravitational constant in lbm*ft/(lbf*s^2).

        Returns:
        float: Computed STO distance in feet.
        """
        def v_stall(weight, rho, s, cl_max):
            return np.sqrt(weight / (0.5 * rho * s * cl_max))

        def sto_integral(v, a, b):
            return v / (a - b * v ** 2)

        v_stall_value = v_stall(weight, rho, s, cl_max)
        v_to = 1.2 * v_stall_value
        a = (gc * thrust) / weight
        b = (gc / weight) * (0.5 * rho * s * cd)

        sto_distance, _ = quad(sto_integral, 0, v_to, args=(a, b))
        return sto_distance

def main():
    """
    Entry point for the application, sets up the MVC components and starts the GUI.
    """
    app = QtWidgets.QApplication([])
    view = TakeoffView()
    model = TakeoffModel()
    controller = TakeoffController(view, model)
    view.show()
    app.exec()

if __name__ == "__main__":
    main()

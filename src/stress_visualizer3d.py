import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

from cube_renderer import CubeRenderer
from stress_plotter import StressPlotter
from stress_tensor import StressTensor


class StressVisualizer:
    def __init__(self):
        # Create figure and axes
        self.fig = plt.figure(figsize=(14, 8))
        self.ax_3d = self.fig.add_subplot(121, projection='3d')
        self.ax_plot = self.fig.add_subplot(122)

        # Initialize components
        self.stress_tensor = StressTensor()
        self.cube_renderer = CubeRenderer(self.ax_3d)
        self.stress_plotter = StressPlotter(self.ax_plot)

        # Add sliders for rotation angles
        self.create_sliders()

    def create_sliders(self):
        """Create sliders for controlling rotation angles."""

        # Add sliders for rotation angles
        self.ax_theta_x = plt.axes([0.2, 0.02, 0.6, 0.03])
        self.ax_theta_y = plt.axes([0.2, 0.06, 0.6, 0.03])
        self.ax_theta_z = plt.axes([0.2, 0.10, 0.6, 0.03])

        self.theta_x_slider = Slider(
            ax=self.ax_theta_x, label='θx (degrees)', valmin=0, valmax=360,
            valinit=0)
        self.theta_y_slider = Slider(
            ax=self.ax_theta_y, label='θy (degrees)', valmin=0, valmax=360,
            valinit=0)
        self.theta_z_slider = Slider(
            ax=self.ax_theta_z, label='θz (degrees)', valmin=0, valmax=360,
            valinit=0)

        self.theta_x_slider.on_changed(self.update)
        self.theta_y_slider.on_changed(self.update)
        self.theta_z_slider.on_changed(self.update)

    def update(self, val):
        theta_x = self.theta_x_slider.val
        theta_y = self.theta_y_slider.val
        theta_z = self.theta_z_slider.val

        # Get rotation matrix and transformed stress tensor
        R = self.stress_tensor.get_rotation_matrix(theta_x, theta_y, theta_z)
        transformed_stress = self.stress_tensor.transform_stress_tensor(
            theta_x, theta_y, theta_z)

        # Update 3D visualization
        self.cube_renderer.draw_cube(R, transformed_stress)

        # Calculate and plot stress variation
        stress_data = self.stress_tensor.calculate_stress_variation(
            theta_x, theta_z)
        self.stress_plotter.plot_stress_variation(stress_data, theta_y)

        # Add text with current stress values
        stress_text = self.stress_tensor.get_stress_text(
            transformed_stress, theta_x, theta_y, theta_z)
        self.stress_plotter.add_stress_text(stress_text)

        plt.draw()

    def show(self):
        """Display the visualization"""
        self.fig.subplots_adjust(bottom=0.15, top=0.95, wspace=0.3)
        self.update(None)  # Initial update
        plt.show()


if __name__ == "__main__":
    visualizer = StressVisualizer()
    visualizer.show()

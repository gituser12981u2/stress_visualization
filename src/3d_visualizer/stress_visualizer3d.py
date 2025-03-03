from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, CheckButtons

from cube_renderer import CubeRenderer
from stress_plotter import StressPlotter
from stress_tensor import StressTensor


class StressVisualizer:
    def __init__(self):
        # Create figure and axes
        self.fig = plt.figure(figsize=(14, 8))

        # Flag for graph visibility
        self.show_graph = True

        # Use the full figure initially with two panels
        self.setup_dual_panel_layout()

        # Initialize components
        self.stress_tensor = StressTensor()
        self.cube_renderer = CubeRenderer(self.ax_3d)
        self.stress_plotter = StressPlotter(self.ax_plot)

        # Add sliders for rotation angles
        self.create_sliders()

        # Add checkbox for toggling graph
        self.create_checkbox()

    def setup_dual_panel_layout(self):
        """Set up the dual panel layout with 3D view and plot."""
        self.fig.clear()
        self.gs = GridSpec(2, 2, height_ratios=[5, 1], figure=self.fig)
        self.ax_3d = self.fig.add_subplot(self.gs[0, 0], projection='3d')
        self.ax_plot = self.fig.add_subplot(self.gs[0, 1])

    def setup_single_panel_layout(self):
        """Set up the single panel layout with just 3D view."""
        self.fig.clear()
        self.gs_full = GridSpec(2, 1, height_ratios=[5, 1], figure=self.fig)
        # Use a larger subplot that takes up more of the figure
        self.ax_3d = self.fig.add_subplot(self.gs_full[0, 0], projection='3d')
        self.ax_plot = None

    def create_checkbox(self):
        """Create checkbox to toggle the graph visibility."""
        self.ax_checkbox = plt.axes([0.01, 0.15, 0.1, 0.05])
        self.check = CheckButtons(
            self.ax_checkbox, ['Show Graph'], [self.show_graph])
        self.check.on_clicked(self.toggle_graph)

    def toggle_graph(self, label):
        """Toggle visibility of the stress vs. angle graph and adjust layout"""
        self.show_graph = not self.show_graph

        # Clear the figure and recreate the layout
        self.fig.clear()

        if self.show_graph:
            # Two-panel layout when graph is visible
            self.gs = GridSpec(2, 2, height_ratios=[5, 1], figure=self.fig)
            self.ax_3d = self.fig.add_subplot(self.gs[0, 0], projection='3d')
            self.ax_plot = self.fig.add_subplot(self.gs[0, 1])
            self.cube_renderer.ax = self.ax_3d
            self.stress_plotter.ax = self.ax_plot
        else:
            # Single-panel layout when graph is hidden
            self.gs_full = GridSpec(2, 1, height_ratios=[
                                    5, 1], figure=self.fig)
            self.ax_3d = self.fig.add_subplot(
                self.gs_full[0, 0], projection='3d')
            self.ax_plot = None
            self.cube_renderer.ax = self.ax_3d

        # Recreate the sliders and checkbox
        self.create_sliders()
        self.create_checkbox()

        self.fig.tight_layout(rect=[0, 0.15, 1, 0.95])

        # Update the visualization
        self.update(None)

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

        # Update 2D plot only if visible
        if self.show_graph and self.ax_plot is not None:
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
        self.fig.tight_layout(rect=[0, 0.15, 1, 0.95])
        self.update(None)  # Initial update
        plt.show()


if __name__ == "__main__":
    visualizer = StressVisualizer()
    visualizer.show()

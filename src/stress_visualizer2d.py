import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import matplotlib.gridspec as gridspec


class StressVisualizer2D:
    def __init__(self):
        self.sigma_x = 1.0
        self.tau_xy = 0.5

        self.fig = plt.figure(figsize=(15, 8))
        gs = gridspec.GridSpec(2, 2, height_ratios=[4, 1])

        # Create subplots
        self.ax_element = plt.subplot(gs[0, 0])  # Element visualization
        self.ax_graph = plt.subplot(gs[0, 1])    # Stress vs angle plot

        self.ax_slider = plt.subplot(gs[1, :])
        self.theta_slider = Slider(
            ax=self.ax_slider,
            label='Rotation Angle θ (degrees)',
            valmin=0,
            valmax=180,
            valinit=0,
            slidermin=None,
            slidermax=None
        )
        self.theta_slider.on_changed(self.update)

    def calculate_stresses(self, theta_deg):
        """Calculate transformed stresses for a given angle"""
        theta = np.radians(theta_deg)

        # Normal stress
        sigma_x_prime = (self.sigma_x/2) + \
            (self.sigma_x/2) * np.cos(2*theta) + \
            self.tau_xy * np.sin(2*theta)

        # Shear stress
        tau_xy_prime = -self.sigma_x/2 * np.sin(2*theta) + \
            self.tau_xy * np.cos(2*theta)

        return sigma_x_prime, tau_xy_prime

    def draw_element(self, angle_deg):
        """Draw the rotated element with stress arrows"""
        self.ax_element.clear()
        angle = np.radians(angle_deg)

        # Create square element
        size = 1
        corners = np.array([
            [-size/2, -size/2],
            [size/2, -size/2],
            [size/2, size/2],
            [-size/2, size/2]
        ])

        # Rotate corners
        rotation_matrix = np.array([
            [np.cos(angle), -np.sin(angle)],
            [np.sin(angle), np.cos(angle)]
        ])
        rotated_corners = corners @ rotation_matrix.T

        # Draw element
        self.ax_element.fill(rotated_corners[:, 0], rotated_corners[:, 1],
                             facecolor='none', edgecolor='black')

        # Calculate stresses at current angle
        sigma_prime, tau_prime = self.calculate_stresses(angle_deg)

        # Draw stress arrows
        arrow_scale = 0.3
        center = np.array([0, 0])

        # Get rotated directions
        # Direction for horizontal faces
        normal_x = rotation_matrix @ np.array([1, 0])
        # Direction for vertical faces
        normal_y = rotation_matrix @ np.array([0, 1])

        # Normal stress arrows (red) for horizontal faces
        if abs(sigma_prime) > 0.01:
            # Left and right faces
            self.ax_element.arrow(
                center[0] - normal_x[0]*size/2, center[1] - normal_x[1]*size/2,
                -normal_x[0]*arrow_scale*sigma_prime, -
                normal_x[1]*arrow_scale*sigma_prime,
                color='red', width=0.02, head_width=0.1
            )
            self.ax_element.arrow(
                center[0] + normal_x[0]*size/2, center[1] + normal_x[1]*size/2,
                normal_x[0]*arrow_scale*sigma_prime, normal_x[1] *
                arrow_scale*sigma_prime,
                color='red', width=0.02, head_width=0.1
            )

            # Top and bottom faces
            self.ax_element.arrow(
                center[0] - normal_y[0]*size/2, center[1] - normal_y[1]*size/2,
                -normal_y[0]*arrow_scale*sigma_prime, -
                normal_y[1]*arrow_scale*sigma_prime,
                color='red', width=0.02, head_width=0.1
            )
            self.ax_element.arrow(
                center[0] + normal_y[0]*size/2, center[1] + normal_y[1]*size/2,
                normal_y[0]*arrow_scale*sigma_prime, normal_y[1] *
                arrow_scale*sigma_prime,
                color='red', width=0.02, head_width=0.1
            )

        # Shear stress arrows (blue) following right-hand rule
        if abs(tau_prime) > 0.01:
            # Left face (clockwise)
            self.ax_element.arrow(
                center[0] - normal_x[0]*size/2, center[1] - normal_x[1]*size/2,
                -normal_y[0]*arrow_scale*tau_prime, -
                normal_y[1]*arrow_scale*tau_prime,
                color='blue', width=0.02, head_width=0.1
            )

            # Right face (counterclockwise)
            self.ax_element.arrow(
                center[0] + normal_x[0]*size/2, center[1] + normal_x[1]*size/2,
                normal_y[0]*arrow_scale *
                tau_prime, normal_y[1]*arrow_scale*tau_prime,
                color='blue', width=0.02, head_width=0.1
            )

            # Bottom face (clockwise)
            self.ax_element.arrow(
                center[0] - normal_y[0]*size/2, center[1] - normal_y[1]*size/2,
                normal_x[0]*arrow_scale *
                tau_prime, normal_x[1]*arrow_scale*tau_prime,
                color='blue', width=0.02, head_width=0.1
            )

            # Top face (counterclockwise)
            self.ax_element.arrow(
                center[0] + normal_y[0]*size/2, center[1] + normal_y[1]*size/2,
                -normal_x[0]*arrow_scale*tau_prime, -
                normal_x[1]*arrow_scale*tau_prime,
                color='blue', width=0.02, head_width=0.1
            )

        self.ax_element.set_aspect('equal')
        self.ax_element.set_xlim(-1.5, 1.5)
        self.ax_element.set_ylim(-1.5, 1.5)

        # Add detailed grid
        self.ax_element.grid(True, which='major',
                             linestyle='-', linewidth='0.5', color='gray')
        self.ax_element.grid(True, which='minor', linestyle=':',
                             linewidth='0.5', color='gray', alpha=0.5)
        self.ax_element.minorticks_on()

        # Add legend
        if abs(sigma_prime) > 0.01 or abs(tau_prime) > 0.01:
            self.ax_element.plot([], [], 'r-', label='Normal stress (σ)')
            self.ax_element.plot([], [], 'b-', label='Shear stress (τ)')
            self.ax_element.legend(loc='upper right')

        self.ax_element.set_title(f'Element at θ = {int(angle_deg)}°\n'
                                  f'σx\' = {sigma_prime:.2f}, τxy\' = {tau_prime:.2f}')

    def draw_graph(self, current_angle):
        """Draw stress vs angle plot"""
        self.ax_graph.clear()

        # Calculate stresses for all angles
        angles = np.linspace(0, 180, 181)
        normal_stresses = []
        shear_stresses = []

        for angle in angles:
            sigma, tau = self.calculate_stresses(angle)
            normal_stresses.append(sigma)
            shear_stresses.append(tau)

        # Plot stress variations
        self.ax_graph.plot(angles, normal_stresses, 'r-', label='σx\'')
        self.ax_graph.plot(angles, shear_stresses, 'b-', label='τxy\'')

        # Mark current angle
        current_sigma, current_tau = self.calculate_stresses(current_angle)
        self.ax_graph.plot(current_angle, current_sigma, 'ro')
        self.ax_graph.plot(current_angle, current_tau, 'bo')
        self.ax_graph.axvline(x=current_angle, color='k',
                              linestyle='--', alpha=0.5)

        # Add detailed grid
        self.ax_graph.grid(True, which='major', linestyle='-',
                           linewidth='0.5', color='gray')
        self.ax_graph.grid(True, which='minor', linestyle=':',
                           linewidth='0.5', color='gray', alpha=0.5)
        self.ax_graph.minorticks_on()

        # Customize plot
        self.ax_graph.set_xlabel('Angle θ (degrees)')
        self.ax_graph.set_ylabel('Stress')
        self.ax_graph.legend()
        self.ax_graph.set_title('Stress Transformation')

    def update(self, val):
        """Update visualization when slider changes"""
        angle = self.theta_slider.val
        self.draw_element(angle)
        self.draw_graph(angle)
        plt.draw()

    def show(self):
        """Display the visualization"""
        self.update(0)
        plt.show()


if __name__ == "__main__":
    viz = StressVisualizer2D()
    viz.show()

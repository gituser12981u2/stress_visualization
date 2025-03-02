import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import matplotlib.gridspec as gridspec


class StressVisualizer2D:
    def __init__(self):
        self.sigma_x = -80.0
        self.sigma_y = 50.0
        self.tau_xy = -25.0

        # Current angle
        self.current_angle = 0

        self.fig = plt.figure(figsize=(15, 8))
        gs = gridspec.GridSpec(3, 2, height_ratios=[4, 0.5, 0.5])

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
            valstep=1.0,
            slidermin=None,
            slidermax=None
        )
        self.theta_slider.on_changed(self.update)

    def calculate_stresses(self, theta_deg):
        """Calculate transformed stresses for a given angle"""
        theta = np.radians(theta_deg)

        # Calculate average and half-difference of normal stresses
        avg_normal = (self.sigma_x + self.sigma_y) / 2
        half_diff = (self.sigma_x - self.sigma_y) / 2

        # Normal stress's
        sigma_x_prime = avg_normal + half_diff * \
            np.cos(2*theta) + self.tau_xy * np.sin(2*theta)
        sigma_y_prime = avg_normal - half_diff * \
            np.cos(2*theta) - self.tau_xy * np.sin(2*theta)

        # Shear stress
        tau_xy_prime = -half_diff * \
            np.sin(2*theta) + self.tau_xy * np.cos(2*theta)

        return sigma_x_prime, sigma_y_prime, tau_xy_prime

    def calculate_max_stress_magnitude(self):
        """Calculate the maximum stress magnitude across all angles"""
        angles = np.linspace(0, 180, 181)
        max_magnitude = 0

        for angle in angles:
            sigma_x, sigma_y, tau = self.calculate_stresses(angle)
            max_magnitude = max(max_magnitude, abs(
                sigma_x), abs(sigma_y), abs(tau))

        return max_magnitude

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
        sigma_x_prime, sigma_y_prime, tau_prime = self.calculate_stresses(
            angle_deg)

        # Dynamically scale arrows
        max_stress = self.calculate_max_stress_magnitude()

        # Draw stress arrows
        base_arrow_scale = 0.8
        element_size = 1.0

        # Compute adaptive scale
        adaptive_scale = min(
            base_arrow_scale,
            element_size / (2 * max(abs(max_stress), 1))
        )

        center = np.array([0, 0])

        # Get rotated directions
        # Direction for horizontal faces
        normal_x = rotation_matrix @ np.array([1, 0])
        # Direction for vertical faces
        normal_y = rotation_matrix @ np.array([0, 1])

        # Normal stress arrows (red) for horizontal faces
        if abs(sigma_x_prime) > 0.01:
            # Left and right faces
            self.ax_element.arrow(
                center[0] - normal_x[0]*size/2, center[1] - normal_x[1]*size/2,
                -normal_x[0]*adaptive_scale*sigma_x_prime, -
                normal_x[1]*adaptive_scale*sigma_x_prime,
                color='red', width=0.02, head_width=0.1
            )
            self.ax_element.arrow(
                center[0] + normal_x[0]*size/2, center[1] + normal_x[1]*size/2,
                normal_x[0]*adaptive_scale*sigma_x_prime, normal_x[1] *
                adaptive_scale*sigma_x_prime,
                color='red', width=0.02, head_width=0.1
            )

        # Normal stress arrows for y-faces (green)
        if abs(sigma_y_prime) > 0.01:
            # Top and bottom faces (y-normal)
            self.ax_element.arrow(
                center[0] - normal_y[0]*size/2, center[1] - normal_y[1]*size/2,
                -normal_y[0]*adaptive_scale*sigma_y_prime, -
                normal_y[1]*adaptive_scale*sigma_y_prime,
                color='green', width=0.02, head_width=0.1
            )
            self.ax_element.arrow(
                center[0] + normal_y[0]*size/2, center[1] + normal_y[1]*size/2,
                normal_y[0]*adaptive_scale*sigma_y_prime, normal_y[1] *
                adaptive_scale*sigma_y_prime,
                color='green', width=0.02, head_width=0.1
            )

        # Shear stress arrows (blue) following right-hand rule
        if abs(tau_prime) > 0.01:
            # Left face (clockwise)
            self.ax_element.arrow(
                center[0] - normal_x[0]*size/2, center[1] - normal_x[1]*size/2,
                -normal_y[0]*adaptive_scale*tau_prime, -
                normal_y[1]*adaptive_scale*tau_prime,
                color='blue', width=0.02, head_width=0.1
            )

            # Right face (counterclockwise)
            self.ax_element.arrow(
                center[0] + normal_x[0]*size/2, center[1] + normal_x[1]*size/2,
                normal_y[0]*adaptive_scale *
                tau_prime, normal_y[1]*adaptive_scale*tau_prime,
                color='blue', width=0.02, head_width=0.1
            )

            # Bottom face (clockwise)
            self.ax_element.arrow(
                center[0] - normal_y[0]*size/2, center[1] - normal_y[1]*size/2,
                normal_x[0]*adaptive_scale *
                tau_prime, normal_x[1]*adaptive_scale*tau_prime,
                color='blue', width=0.02, head_width=0.1
            )

            # Top face (counterclockwise)
            self.ax_element.arrow(
                center[0] + normal_y[0]*size/2, center[1] + normal_y[1]*size/2,
                -normal_x[0]*adaptive_scale*tau_prime, -
                normal_x[1]*adaptive_scale*tau_prime,
                color='blue', width=0.02, head_width=0.1
            )

        view_margin = 1.5 + max(1.0, abs(max_stress) * adaptive_scale * 1.2)
        self.ax_element.set_aspect('equal')
        self.ax_element.set_xlim(-view_margin, view_margin)
        self.ax_element.set_ylim(-view_margin, view_margin)

        # Add detailed grid
        self.ax_element.grid(True, which='major',
                             linestyle='-', linewidth='0.5', color='gray')
        self.ax_element.grid(True, which='minor', linestyle=':',
                             linewidth='0.5', color='gray', alpha=0.5)
        self.ax_element.minorticks_on()

        # Add legend
        if abs(sigma_x_prime) > 0.01 or abs(tau_prime) > 0.01:
            self.ax_element.plot([], [], 'r-', label='Normal stress (σ)')
            self.ax_element.plot([], [], 'g-', label='Normal stress y (σy)')
            self.ax_element.plot([], [], 'b-', label='Shear stress (τ)')
            self.ax_element.legend(loc='upper right')

        self.ax_element.set_title(f'Element at θ = {int(angle_deg)}°\n'
                                  f'σx\' = {sigma_x_prime:.2f}, σy\' = {sigma_y_prime:.2f}, τxy\' = {tau_prime:.2f}')

    def draw_graph(self, current_angle):
        """Draw stress vs angle plot"""
        self.ax_graph.clear()

        # Calculate stresses for all angles
        angles = np.linspace(0, 180, 181)
        sigma_x_values = []
        sigma_y_values = []
        tau_values = []

        for angle in angles:
            sigma_x, sigma_y, tau = self.calculate_stresses(angle)
            sigma_x_values.append(sigma_x)
            sigma_y_values.append(sigma_y)
            tau_values.append(tau)

        # Plot stress variations
        self.ax_graph.plot(angles, sigma_x_values, 'r-', label='σx\'')
        self.ax_graph.plot(angles, sigma_y_values, 'g-', label='σy\'')
        self.ax_graph.plot(angles, tau_values, 'b-', label='τxy\'')

        # Mark current angle
        current_sigma_x, current_sigma_y, current_tau = self.calculate_stresses(
            current_angle)
        self.ax_graph.plot(current_angle, current_sigma_x, 'ro')
        self.ax_graph.plot(current_angle, current_sigma_y, 'go')
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

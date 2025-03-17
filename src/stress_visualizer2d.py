"""
Stress Tensor 2D Visualization Module

This module provides an interactive visualization of 2D stress tensors using
a cube representation with color-coded stress arrows. It allows one to rotate
the element around three axes and observe how stress components transform
under rotation.

The visualization includes:
1. A 2D element with color-coded stress arrows
2. Stress variation plots showing how components change with rotation angle
3. Principal stress calculation and visualization
3. Interactive sliders to control rotation angles

Classes:
    StressVisualizer2D: Main class for the 2D stress visualization application

Usage:
    Run this script directly to launch the visualization

    $ python stress_visualizer2d.py
"""


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import matplotlib.gridspec as gridspec


class StressVisualizer2D:
    """
    Interactive 2D stress transformation visualization application.

    This class creates a matplotlib-based GUI that visualizes stress
    transformation in 2D--plane stress. It displays a square element with
    stress arrows that update as the element is rotated, and plots tress
    variation with angle.

    Attributes:
        sigma_x (float): Normal stress in x-direction
        sigma_y (float): Normal stress in y-direction
        tau_xy (float): Shear stress
        current_angle (float): Current rotation angle in degrees
        sigma_1 (float): Maximum principal stress
        sigma_2 (float): Minimum principal stress
        principal_angle (float): Principal angle in degrees
        fig (matplotlib.figure.Figure): The main figure containing all plots
        ax_element (matplotlib.axes.Axes): Axes for the rotated element
        ax_graph (matplotlib.axes.Axes): Axes for the stress variation plot
    """

    def __init__(self):
        """Initialize the 2D stress visualization with default values."""
        self.sigma_x = 100.0
        self.sigma_y = 10.0
        self.tau_xy = 60.0

        # Current angle
        self.current_angle = 0

        # Calculate the principal stresses and angle
        self.calculate_principal_values()

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
            slidermin=None,
            slidermax=None
        )
        self.theta_slider.on_changed(self.update)

    def calculate_stresses(self, theta_deg):
        """
        Calculate transformed stresses for a given angle.

        Args:
            theta_deg (float): Rotation angle in degrees

        Returns:
            tuple: (sigma_x_prime, sigma_y_prime, tau_xy_prime) transformed
            stresses
        """
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
        """
        Calculate the maximum stress magnitude across all angles.

        This is used for scaling the stress arrows appropriately.

        Returns:
            float: Maximum absolute stress value
        """
        angles = np.linspace(0, 180, 181)
        max_magnitude = 0

        for angle in angles:
            sigma_x, sigma_y, tau = self.calculate_stresses(angle)
            max_magnitude = max(max_magnitude, abs(
                sigma_x), abs(sigma_y), abs(tau))

        return max_magnitude

    def calculate_principal_values(self):
        """
        Calculate principal stresses and principal angle.

        Sets:
            self.sigma_1: Maximum principal stress
            self.sigma2: Minimum principal stress
            self.principal_angle: Principal angle in degrees
        """
        avg_normal = (self.sigma_x + self.sigma_y) / 2
        half_diff = (self.sigma_x - self.sigma_y) / 2

        # Principal stresses
        r = np.sqrt(half_diff**2 + self.tau_xy**2)
        self.sigma_1 = avg_normal + r
        self.sigma_2 = avg_normal - r

        # Principal angle
        if abs(half_diff) < 1e-10:
            self.principal_angle = 45.0 if self.tau_xy > 0 else -45.0
        else:
            self.principal_angle = np.degrees(
                0.5 * np.arctan2(self.tau_xy, half_diff))

        # Make angle between 0 and 90 degrees
        if self.principal_angle < 0:
            self.principal_angle += 90

        self.principal_angle = self.principal_angle % 90

    def draw_element(self, angle_deg):
        """
        Draw the rotated element with stress arrows.

        Args:
            angle_deg (float): Rotation angle in degrees
        """
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
                -normal_y[0]*adaptive_scale*tau_prime,
                -normal_y[1]*adaptive_scale*tau_prime,
                color='blue', width=0.02, head_width=0.1
            )

            # Right face (counterclockwise)
            self.ax_element.arrow(
                center[0] + normal_x[0]*size/2, center[1] + normal_x[1]*size/2,
                normal_y[0]*adaptive_scale*tau_prime,
                normal_y[1]*adaptive_scale*tau_prime,
                color='blue', width=0.02, head_width=0.1
            )

            # Bottom face (clockwise)
            self.ax_element.arrow(
                center[0] - normal_y[0]*size/2, center[1] - normal_y[1]*size/2,
                -normal_x[0]*adaptive_scale*tau_prime,
                -normal_x[1]*adaptive_scale*tau_prime,
                color='blue', width=0.02, head_width=0.1
            )

            # Top face (counterclockwise)
            self.ax_element.arrow(
                center[0] + normal_y[0]*size/2, center[1] + normal_y[1]*size/2,
                normal_x[0]*adaptive_scale*tau_prime,
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
            self.ax_element.plot([], [], 'r-', label='Normal stress x (σx)')
            self.ax_element.plot([], [], 'g-', label='Normal stress y (σy)')
            self.ax_element.plot([], [], 'b-', label='Shear stress (τxy)')
            self.ax_element.legend(loc='upper right')

        self.ax_element.set_title(f'Element at θ = {int(angle_deg)}°\n'
                                  f'σx\' = {sigma_x_prime:.2f}, σy\' = {sigma_y_prime:.2f}, τxy\' = {tau_prime:.2f}')

    def draw_graph(self, current_angle):
        """
        Draw stress vs angle plot.

        Args:
            current_angle (float): Current rotation angle to highlight on the
            plot
        """
        self.ax_graph.clear()

        # Calculate stresses for all angles
        angles = np.linspace(0, 180, 361)
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

        # Add horizontal lines for principal stresses
        self.ax_graph.axhline(y=self.sigma_1, color='m', linestyle='--',
                              alpha=0.7, label=f'σ₁ = {self.sigma_1:.2f}')
        self.ax_graph.axhline(y=self.sigma_2, color='c', linestyle='--',
                              alpha=0.7, label=f'σ₂ = {self.sigma_2:.2f}')

        # Add vertical lines for principal angle
        self.ax_graph.axvline(x=self.principal_angle, color='m', linestyle=':', alpha=0.7,
                              label=f'θp = {self.principal_angle:.2f}°')
        self.ax_graph.axvline(x=self.principal_angle+90,
                              color='c', linestyle=':', alpha=0.7)

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

        # Display principal values in text box
        principal_text = (f"Principal Stresses:\n"
                          f"σ₁ = {self.sigma_1:.2f}\n"
                          f"σ₂ = {self.sigma_2:.2f}\n"
                          f"Principal Angle (θp) = {self.principal_angle:.2f}°")

        self.ax_graph.text(0.02, 0.02, principal_text,
                           transform=self.ax_graph.transAxes,
                           bbox=dict(facecolor='white', alpha=0.7), fontsize=9)

    def update(self, val):
        """
        Update visualization when slider changes.

        Args:
            val: The slider value (unused but required by matplotlib)
        """
        angle = self.theta_slider.val
        self.draw_element(angle)
        self.draw_graph(angle)
        plt.draw()

    def show(self):
        """Display the visualization."""
        self.update(0)
        plt.show()


if __name__ == "__main__":
    viz = StressVisualizer2D()
    viz.show()

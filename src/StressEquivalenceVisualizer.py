import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from matplotlib.patches import Polygon
import matplotlib.gridspec as gridspec


class StressEquivalenceVisualizer:
    """
    Visualizer to demonstrate why stress transformations work by showing
    force equilibrium on different cutting planes through a stress element.
    """

    def __init__(self):
        self.sigma_x = -0.3
        self.sigma_y = 1.0
        self.tau_xy = 0.2

        # Current angle and visualization settings
        self.current_angle = 30.0  # Start with a non-zero angle
        self.show_force_components = True
        self.show_force_balance = True

        # Calculate principal stresses and angle
        self.calculate_principal_values()

        # Create separate figures for each visualization
        self.fig_original = plt.figure(figsize=(8, 7), num="Original Element")
        self.ax_original = self.fig_original.add_subplot(111)

        self.fig_rotated = plt.figure(
            figsize=(8, 7), num="Transformed Element")
        self.ax_rotated = self.fig_rotated.add_subplot(111)

        self.fig_fbd = plt.figure(figsize=(10, 8), num="Force Balance")
        gs_fbd = gridspec.GridSpec(3, 1, height_ratios=[4, 0.5, 0.5])
        self.ax_fbd = plt.subplot(gs_fbd[0])

        # Controls on the force balance figure
        self.ax_angle_slider = plt.subplot(gs_fbd[1])
        self.theta_slider = Slider(
            ax=self.ax_angle_slider,
            label='Rotation Angle θ (degrees)',
            valmin=0,
            valmax=180,
            valinit=self.current_angle
        )
        self.theta_slider.on_changed(self.update)

        # Create a panel for buttons and explanation
        button_panel = plt.subplot(gs_fbd[2])
        button_panel.axis('off')

        # Add buttons
        # button_width = 0.2
        # button_height = 0.6

        # self.components_button = Button(
        #     plt.axes([0.2, 0.15, button_width, button_height]),
        #     'Toggle Components',
        #     color='lightgoldenrodyellow',
        #     hovercolor='0.975'
        # )
        # self.components_button.on_clicked(self.toggle_components)

        # self.balance_button = Button(
        #     plt.axes([0.5, 0.15, button_width, button_height]),
        #     'Toggle Force Balance',
        #     color='lightgoldenrodyellow',
        #     hovercolor='0.975'
        # )
        # self.balance_button.on_clicked(self.toggle_balance)

        # Add explanation text box at the bottom of the FBD figure
        explanation_text = (
            # "WHY STRESS TRANSFORMATIONS WORK: At any point in a body, stresses must balance forces in all directions. "
            # "When cutting through the stress element at any angle, the resulting forces must still balance. "
            # "The transformed stresses represent the force distribution needed to maintain equilibrium."
        )

        self.fig_fbd.text(0.5, 0.02, explanation_text, ha='center', fontsize=9,
                          bbox=dict(facecolor='lightyellow', alpha=0.9, boxstyle='round'))
        # self.ax_explanation.axis('off')

    def calculate_principal_values(self):
        """Calculate principal stresses and principal angle"""
        avg_normal = (self.sigma_x + self.sigma_y) / 2
        half_diff = (self.sigma_x - self.sigma_y) / 2

        # Principal stresses
        r = np.sqrt(half_diff**2 + self.tau_xy**2)
        self.sigma_1 = avg_normal + r
        self.sigma_2 = avg_normal - r

        # Principal angle (in degrees)
        if abs(half_diff) < 1e-10:  # Prevent division by zero
            self.principal_angle = 45.0 if self.tau_xy > 0 else -45.0
        else:
            self.principal_angle = np.degrees(
                0.5 * np.arctan2(self.tau_xy, half_diff))

        # Make angle between 0 and 90 degrees for conventional representation
        if self.principal_angle < 0:
            self.principal_angle += 90

        # Ensure principal_angle is within [0, 90] range
        self.principal_angle = self.principal_angle % 90

    def calculate_stresses(self, theta_deg):
        """Calculate transformed stresses for a given angle"""
        theta = np.radians(theta_deg)

        # Calculate average and half-difference of normal stresses
        avg_normal = (self.sigma_x + self.sigma_y) / 2
        half_diff = (self.sigma_x - self.sigma_y) / 2

        # Normal stresses
        sigma_x_prime = avg_normal + half_diff * \
            np.cos(2*theta) + self.tau_xy * np.sin(2*theta)
        sigma_y_prime = avg_normal - half_diff * \
            np.cos(2*theta) - self.tau_xy * np.sin(2*theta)

        # Shear stress
        tau_xy_prime = -half_diff * \
            np.sin(2*theta) + self.tau_xy * np.cos(2*theta)

        return sigma_x_prime, sigma_y_prime, tau_xy_prime

    def draw_original_element(self):
        """Draw the element in its original orientation with stresses"""
        self.ax_original.clear()

        # Create square element
        size = 1
        corners = np.array([
            [-size/2, -size/2],
            [size/2, -size/2],
            [size/2, size/2],
            [-size/2, size/2]
        ])

        # Draw element
        self.ax_original.add_patch(
            Polygon(corners, fill=False, edgecolor='black'))

        # Determine scale factor for arrows - improved scaling
        max_stress = max(abs(self.sigma_x), abs(
            self.sigma_y), abs(self.tau_xy))
        arrow_scale = 0.2 / max(1.0, max_stress/100)

        center = np.array([0, 0])

        # Draw normal stress arrows on x-faces
        if abs(self.sigma_x) > 0.01:
            # Left face
            self.ax_original.arrow(
                -size/2, 0,
                -arrow_scale*self.sigma_x, 0,
                head_width=0.1, head_length=0.1, fc='red', ec='red', width=0.02
            )
            # Right face
            self.ax_original.arrow(
                size/2, 0,
                arrow_scale*self.sigma_x, 0,
                head_width=0.1, head_length=0.1, fc='red', ec='red', width=0.02
            )

        # Draw normal stress arrows on y-faces
        if abs(self.sigma_y) > 0.01:
            # Bottom face
            self.ax_original.arrow(
                0, -size/2,
                0, -arrow_scale*self.sigma_y,
                head_width=0.1, head_length=0.1, fc='green', ec='green', width=0.02
            )
            # Top face
            self.ax_original.arrow(
                0, size/2,
                0, arrow_scale*self.sigma_y,
                head_width=0.1, head_length=0.1, fc='green', ec='green', width=0.02
            )

        # Draw shear stress arrows
        if abs(self.tau_xy) > 0.01:
            # Left face
            self.ax_original.arrow(
                -size/2, size/4,
                0, arrow_scale*self.tau_xy,
                head_width=0.1, head_length=0.1, fc='blue', ec='blue', width=0.02
            )
            self.ax_original.arrow(
                -size/2, -size/4,
                0, -arrow_scale*self.tau_xy,
                head_width=0.1, head_length=0.1, fc='blue', ec='blue', width=0.02
            )

            # Right face
            self.ax_original.arrow(
                size/2, size/4,
                0, -arrow_scale*self.tau_xy,
                head_width=0.1, head_length=0.1, fc='blue', ec='blue', width=0.02
            )
            self.ax_original.arrow(
                size/2, -size/4,
                0, arrow_scale*self.tau_xy,
                head_width=0.1, head_length=0.1, fc='blue', ec='blue', width=0.02
            )

            # Bottom face
            self.ax_original.arrow(
                -size/4, -size/2,
                -arrow_scale*self.tau_xy, 0,
                head_width=0.1, head_length=0.1, fc='blue', ec='blue', width=0.02
            )
            self.ax_original.arrow(
                size/4, -size/2,
                arrow_scale*self.tau_xy, 0,
                head_width=0.1, head_length=0.1, fc='blue', ec='blue', width=0.02
            )

            # Top face
            self.ax_original.arrow(
                -size/4, size/2,
                arrow_scale*self.tau_xy, 0,
                head_width=0.1, head_length=0.1, fc='blue', ec='blue', width=0.02
            )
            self.ax_original.arrow(
                size/4, size/2,
                -arrow_scale*self.tau_xy, 0,
                head_width=0.1, head_length=0.1, fc='blue', ec='blue', width=0.02
            )

        # Add a line representing the cutting plane at the given angle
        angle_rad = np.radians(self.current_angle)
        dx = np.cos(angle_rad)
        dy = np.sin(angle_rad)

        # Calculate intersection points with the square
        if abs(dx) > abs(dy):
            # Intersect with left/right sides
            t1 = (-size/2 - center[0]) / dx
            t2 = (size/2 - center[0]) / dx

            p1 = center + t1 * np.array([dx, dy])
            p2 = center + t2 * np.array([dx, dy])
        else:
            # Intersect with top/bottom sides
            t1 = (-size/2 - center[1]) / dy
            t2 = (size/2 - center[1]) / dy

            p1 = center + t1 * np.array([dx, dy])
            p2 = center + t2 * np.array([dx, dy])

        # Draw the cutting plane
        self.ax_original.plot([p1[0], p2[0]], [p1[1], p2[1]], 'k--', lw=2)

        # Draw a small angle arc to show the rotation angle
        angle_radius = 0.3
        self.ax_original.plot(
            [0, angle_radius], [0, 0], 'k-', lw=1, alpha=0.7
        )
        arc_angles = np.linspace(0, self.current_angle, 50)
        arc_x = angle_radius * np.cos(np.radians(arc_angles))
        arc_y = angle_radius * np.sin(np.radians(arc_angles))
        self.ax_original.plot(arc_x, arc_y, 'k-', lw=1, alpha=0.7)
        # Add the angle label
        arc_label_angle = self.current_angle / 2
        arc_label_x = angle_radius * 1.2 * np.cos(np.radians(arc_label_angle))
        arc_label_y = angle_radius * 1.2 * np.sin(np.radians(arc_label_angle))
        self.ax_original.text(
            arc_label_x, arc_label_y, f'θ = {self.current_angle:.1f}°',
            ha='center', va='center', fontsize=10
        )

        # Set plot limits and labels
        view_margin = 2.0
        self.ax_original.set_aspect('equal')
        self.ax_original.set_xlim(-view_margin, view_margin)
        self.ax_original.set_ylim(-view_margin, view_margin)
        self.ax_original.grid(True)

        # Add legend and title - with smaller, more compact legend
        self.ax_original.plot([], [], 'r-', label='σₓ')
        self.ax_original.plot([], [], 'g-', label='σᵧ')
        self.ax_original.plot([], [], 'b-', label='τₓᵧ')
        self.ax_original.plot([], [], 'k--', label='Cutting plane')
        self.ax_original.legend(
            loc='upper right', fontsize=9, framealpha=0.7, handlelength=1.5)

        self.ax_original.set_title(
            'Original Stress Element with Cutting Plane')

    def draw_rotated_element(self):
        """Draw the element rotated to align with the cutting plane and display transformed stresses"""
        self.ax_rotated.clear()

        # Get transformed stresses
        sigma_n, sigma_t, tau_nt = self.calculate_stresses(self.current_angle)

        # Create square element
        size = 1
        corners = np.array([
            [-size/2, -size/2],
            [size/2, -size/2],
            [size/2, size/2],
            [-size/2, size/2]
        ])

        # Draw element
        self.ax_rotated.add_patch(
            Polygon(corners, fill=False, edgecolor='black'))

        # Determine scale factor for arrows
        max_stress = max(abs(sigma_n), abs(sigma_t), abs(tau_nt))
        arrow_scale = 0.3 / max(1.0, max_stress/50)

        center = np.array([0, 0])

        # Draw normal stress arrows on x-faces (now representing sigma_n)
        if abs(sigma_n) > 0.01:
            # Left face
            self.ax_rotated.arrow(
                -size/2, 0,
                -arrow_scale*sigma_n, 0,
                head_width=0.1, head_length=0.1, fc='red', ec='red', width=0.02
            )
            # Right face
            self.ax_rotated.arrow(
                size/2, 0,
                arrow_scale*sigma_n, 0,
                head_width=0.1, head_length=0.1, fc='red', ec='red', width=0.02
            )

        # Draw normal stress arrows on y-faces (now representing sigma_t)
        if abs(sigma_t) > 0.01:
            # Bottom face
            self.ax_rotated.arrow(
                0, -size/2,
                0, -arrow_scale*sigma_t,
                head_width=0.1, head_length=0.1, fc='green', ec='green', width=0.02
            )
            # Top face
            self.ax_rotated.arrow(
                0, size/2,
                0, arrow_scale*sigma_t,
                head_width=0.1, head_length=0.1, fc='green', ec='green', width=0.02
            )

        # Draw shear stress arrows (now representing tau_nt)
        if abs(tau_nt) > 0.01:
            # Left face
            self.ax_rotated.arrow(
                -size/2, size/4,
                0, arrow_scale*tau_nt,
                head_width=0.1, head_length=0.1, fc='blue', ec='blue', width=0.02
            )
            self.ax_rotated.arrow(
                -size/2, -size/4,
                0, -arrow_scale*tau_nt,
                head_width=0.1, head_length=0.1, fc='blue', ec='blue', width=0.02
            )

            # Right face
            self.ax_rotated.arrow(
                size/2, size/4,
                0, -arrow_scale*tau_nt,
                head_width=0.1, head_length=0.1, fc='blue', ec='blue', width=0.02
            )
            self.ax_rotated.arrow(
                size/2, -size/4,
                0, arrow_scale*tau_nt,
                head_width=0.1, head_length=0.1, fc='blue', ec='blue', width=0.02
            )

            # Bottom face
            self.ax_rotated.arrow(
                -size/4, -size/2,
                -arrow_scale*tau_nt, 0,
                head_width=0.1, head_length=0.1, fc='blue', ec='blue', width=0.02
            )
            self.ax_rotated.arrow(
                size/4, -size/2,
                arrow_scale*tau_nt, 0,
                head_width=0.1, head_length=0.1, fc='blue', ec='blue', width=0.02
            )

            # Top face
            self.ax_rotated.arrow(
                -size/4, size/2,
                arrow_scale*tau_nt, 0,
                head_width=0.1, head_length=0.1, fc='blue', ec='blue', width=0.02
            )
            self.ax_rotated.arrow(
                size/4, size/2,
                -arrow_scale*tau_nt, 0,
                head_width=0.1, head_length=0.1, fc='blue', ec='blue', width=0.02
            )

        # Add a horizontal line representing the original cutting plane (now horizontal)
        self.ax_rotated.plot([-size/2, size/2], [0, 0], 'k--', lw=2)

        # Set plot limits and labels
        view_margin = 2.0
        self.ax_rotated.set_aspect('equal')
        self.ax_rotated.set_xlim(-view_margin, view_margin)
        self.ax_rotated.set_ylim(-view_margin, view_margin)
        self.ax_rotated.grid(True)

        # Add legend and title - with smaller, more compact legend
        self.ax_rotated.plot([], [], 'r-', label='σₙ')
        self.ax_rotated.plot([], [], 'g-', label='σₜ')
        self.ax_rotated.plot([], [], 'b-', label='τₙₜ')
        self.ax_rotated.plot([], [], 'k--', label='Cutting plane')
        self.ax_rotated.legend(loc='upper right', fontsize=9,
                               framealpha=0.7, handlelength=1.5)

        stress_text = (f"Transformed stresses:\n"
                       f"σₙ = {sigma_n:.2f}\n"
                       f"σₜ = {sigma_t:.2f}\n"
                       f"τₙₜ = {tau_nt:.2f}")

        self.ax_rotated.text(
            0.02, 0.02, stress_text, transform=self.ax_rotated.transAxes,
            bbox=dict(facecolor='white', alpha=0.7), fontsize=10
        )

        self.ax_rotated.set_title('Transformed Stress Element')

    def draw_force_balance(self):
        """Draw free body diagram showing force balance on the cutting plane"""
        self.ax_fbd.clear()

        # Get angle and transformed stresses
        angle_rad = np.radians(self.current_angle)
        sigma_n, sigma_t, tau_nt = self.calculate_stresses(self.current_angle)

        # Define the cutting plane segment length
        plane_length = 2.0

        # Calculate the area factor (for visualization purposes)
        area_factor = 0.5  # Scale the forces for better visualization

        # Draw the cutting plane
        self.ax_fbd.plot([-plane_length/2, plane_length/2], [0, 0], 'k-', lw=2)

        # Calculate force components from original stresses
        normal_from_sigma_x = self.sigma_x * \
            np.cos(angle_rad) * np.cos(angle_rad) * area_factor
        shear_from_sigma_x = self.sigma_x * \
            np.cos(angle_rad) * np.sin(angle_rad) * area_factor

        normal_from_sigma_y = self.sigma_y * \
            np.sin(angle_rad) * np.sin(angle_rad) * area_factor
        shear_from_sigma_y = self.sigma_y * \
            np.sin(angle_rad) * np.cos(angle_rad) * area_factor

        normal_from_tau = 2 * self.tau_xy * \
            np.sin(angle_rad) * np.cos(angle_rad) * area_factor
        shear_from_tau = self.tau_xy * \
            (np.cos(angle_rad)**2 - np.sin(angle_rad)**2) * area_factor

        # Total forces on the plane
        total_normal_force = normal_from_sigma_x + normal_from_sigma_y + normal_from_tau
        total_shear_force = shear_from_sigma_x + shear_from_sigma_y + shear_from_tau

        # Calculate the total transformed stresses (should match those from calculate_stresses)
        transformed_normal = total_normal_force / area_factor
        transformed_shear = total_shear_force / area_factor

        # Scale forces for display
        force_scale = 0.5

        # Draw individual force components if enabled
        if self.show_force_components:
            # Normal force contributions (pointing away from plane is positive)
            components_y = 0.5  # offset for stacking components visually

            # From sigma_x
            self.ax_fbd.arrow(
                0, components_y,
                0, normal_from_sigma_x * force_scale,
                head_width=0.1, head_length=0.1, fc='red', ec='red', width=0.02,
                alpha=0.5, label='From σₓ (normal)'
            )

            # From sigma_y
            self.ax_fbd.arrow(
                0, components_y,
                0, normal_from_sigma_y * force_scale,
                head_width=0.1, head_length=0.1, fc='green', ec='green', width=0.02,
                alpha=0.5, label='From σᵧ (normal)'
            )

            # From tau_xy
            self.ax_fbd.arrow(
                0, components_y,
                0, normal_from_tau * force_scale,
                head_width=0.1, head_length=0.1, fc='blue', ec='blue', width=0.02,
                alpha=0.5, label='From τₓᵧ (normal)'
            )

            # Shear force contributions (right is positive)
            components_x = 0.5  # offset for clarity

            # From sigma_x
            self.ax_fbd.arrow(
                components_x, 0,
                shear_from_sigma_x * force_scale, 0,
                head_width=0.1, head_length=0.1, fc='red', ec='red', width=0.02,
                alpha=0.5, label='From σₓ (shear)'
            )

            # From sigma_y
            self.ax_fbd.arrow(
                components_x, 0,
                shear_from_sigma_y * force_scale, 0,
                head_width=0.1, head_length=0.1, fc='green', ec='green', width=0.02,
                alpha=0.5, label='From σᵧ (shear)'
            )

            # From tau_xy
            self.ax_fbd.arrow(
                components_x, 0,
                shear_from_tau * force_scale, 0,
                head_width=0.1, head_length=0.1, fc='blue', ec='blue', width=0.02,
                alpha=0.5, label='From τₓᵧ (shear)'
            )

        # Draw the resultant forces on the plane
        if self.show_force_balance:
            # Normal force (perpendicular to cutting plane)
            self.ax_fbd.arrow(
                0, 0,
                0, total_normal_force * force_scale,
                head_width=0.1, head_length=0.1, fc='purple', ec='purple', width=0.03,
                label=f'Total normal force (σₙ = {transformed_normal:.2f})'
            )

            # Shear force (along cutting plane)
            self.ax_fbd.arrow(
                0, 0,
                total_shear_force * force_scale, 0,
                head_width=0.1, head_length=0.1, fc='orange', ec='orange', width=0.03,
                label=f'Total shear force (τₙₜ = {transformed_shear:.2f})'
            )

        # Set plot limits and labels
        view_margin = 2.0
        self.ax_fbd.set_aspect('equal')
        self.ax_fbd.set_xlim(-view_margin, view_margin)
        self.ax_fbd.set_ylim(-view_margin, view_margin)
        self.ax_fbd.grid(True)

        # Draw the angle indicator
        self.ax_fbd.plot(
            [-plane_length/2, -plane_length/2], [0, -0.5], 'k-', lw=1, alpha=0.7
        )
        arc_x = -plane_length/2 + 0.3 * np.cos(np.linspace(0, -angle_rad, 50))
        arc_y = 0.3 * np.sin(np.linspace(0, -angle_rad, 50))
        self.ax_fbd.plot(arc_x, arc_y, 'k-', lw=1, alpha=0.7)

        # Draw reference coordinate system
        self.ax_fbd.arrow(
            -plane_length/2, -0.8,
            0.3, 0,
            head_width=0.06, head_length=0.06, fc='black', ec='black'
        )
        self.ax_fbd.arrow(
            -plane_length/2, -0.8,
            0, 0.3,
            head_width=0.06, head_length=0.06, fc='black', ec='black'
        )
        self.ax_fbd.text(-plane_length/2 + 0.35, -0.85, 'x', fontsize=10)
        self.ax_fbd.text(-plane_length/2 - 0.1, -0.5, 'y', fontsize=10)

        # Add explanation text
        stress_balance_text = (
            f"Force balance on cutting plane at θ = {self.current_angle:.1f}°:\n\n"
            f"Normal stress contributions:\n"
            f"  From σₓ = {normal_from_sigma_x/area_factor:.2f}\n"
            f"  From σᵧ = {normal_from_sigma_y/area_factor:.2f}\n"
            f"  From τₓᵧ = {normal_from_tau/area_factor:.2f}\n"
            f"  Total σₙ = {transformed_normal:.2f}\n\n"
            f"Shear stress contributions:\n"
            f"  From σₓ = {shear_from_sigma_x/area_factor:.2f}\n"
            f"  From σᵧ = {shear_from_sigma_y/area_factor:.2f}\n"
            f"  From τₓᵧ = {shear_from_tau/area_factor:.2f}\n"
            f"  Total τₙₜ = {transformed_shear:.2f}"
        )

        # Place text on right side
        self.ax_fbd.text(
            plane_length/2 + 0.1, -view_margin + 0.3, stress_balance_text,
            fontsize=9, bbox=dict(facecolor='white', alpha=0.7)
        )

        # Add compact legend with smaller font
        handles, labels = self.ax_fbd.get_legend_handles_labels()
        # If there are many items, split them into columns
        if len(handles) > 4:
            ncol = 2
        else:
            ncol = 1

        self.ax_fbd.legend(
            handles, labels,
            loc='upper left',
            fontsize=8,
            framealpha=0.7,
            ncol=ncol,
            handlelength=1.2
        )

        # Add title
        self.ax_fbd.set_title('Force Balance on Cutting Plane')

    def draw_explanation(self):
        """Not needed in the new layout as explanation is added directly to the figure"""
        pass

    def toggle_components(self, event):
        """Toggle display of force components"""
        self.show_force_components = not self.show_force_components
        self.update(None)

    def toggle_balance(self, event):
        """Toggle display of force balance"""
        self.show_force_balance = not self.show_force_balance
        self.update(None)

    def update(self, val):
        """Update visualization when slider changes"""
        if val is not None:
            self.current_angle = self.theta_slider.val

        self.draw_original_element()
        self.draw_rotated_element()
        self.draw_force_balance()

        # Update all figures
        self.fig_original.canvas.draw_idle()
        self.fig_rotated.canvas.draw_idle()
        self.fig_fbd.canvas.draw_idle()

    def show(self):
        """Display the visualization with all figures"""
        # Initial draw
        self.update(None)

        # Adjust layout for each figure
        self.fig_original.tight_layout()
        self.fig_rotated.tight_layout()
        self.fig_fbd.tight_layout()

        # Show all figures
        plt.show()


if __name__ == "__main__":
    viz = StressEquivalenceVisualizer()
    viz.show()

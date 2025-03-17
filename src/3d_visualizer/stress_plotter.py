"""
Stress Plotter Module

This module provides functionality for plotting stress variations with respect
to rotation angle. It is used in conjunction with 3D stress visualization to
show how stress components change as an element is rotated.

Classes:
    StressPlotter: Handles 3D plotting of stress variation with rotation angle
"""


class StressPlotter:
    """
    Class to handle the 3D stress plotting and analysis.

    This class creates and updates 3D plots showing how stress components
    (normal and shear stresses) vary as a function of rotation angle.

    Attributes:
        ax (matplotlib.axes.Axes): The matplotlib axes where the plots will be
        drawn
    """

    def __init__(self, ax):
        """
        Initialize with a matplotlib axis.

        Args:
            ax (matplotlib.axes.Axes): the matplotlib axes for plotting
        """
        self.ax = ax

    def plot_stress_variation(self, stress_data, current_angle):
        """
        Plot stress variation with rotation angle

        Args:
            stress_data (dict): Dictionary containing stress data with keys:
                - 'angles': Array of rotation angles
                - 'σx', 'σy', 'σz': Normal stress components
                - 'τxy' 'τxz', 'τyz': Shear stress components
            current_angle (float): The current rotation angle to highlight
        """
        self.ax.clear()

        # Plot normal stresses
        self.ax.plot(stress_data['angles'],
                     stress_data['σx'], 'r-', label='σx')

        # Plot shear stresses
        self.ax.plot(stress_data['angles'],
                     stress_data['τxy'], 'b-', label='τxy')
        self.ax.plot(stress_data['angles'],
                     stress_data['τxz'], 'g-', label='τxz')
        self.ax.plot(stress_data['angles'],
                     stress_data['τyz'], 'purple', label='τyz')

        # Mark current angle
        self.ax.axvline(x=current_angle, color='k', linestyle='--')

        self.ax.grid(True)
        self.ax.set_xlabel('θy (degrees)')
        self.ax.set_ylabel('Stress')
        self.ax.set_title('Stress vs. Rotation Angle θy')
        self.ax.legend()

    def add_stress_text(self, stress_text):
        """
        Add text box with current stress values

        Args:
            stress_text (str): Formatted text string with stress values
        """
        self.ax.text(0.05, 0.05, stress_text, transform=self.ax.transAxes,
                     bbox=dict(facecolor='white', alpha=0.8), fontsize=9)

class StressPlotter:
    """Class to handle the 2D stress plotting and analysis."""

    def __init__(self, ax):
        """Initialize with a matploblib axis"""
        self.ax = ax

    def plot_stress_variation(self, stress_data, current_angle):
        """Plot stress variation with rotation angle"""
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
        """Add text box with current stress values"""
        self.ax.text(0.05, 0.05, stress_text, transform=self.ax.transAxes,
                     bbox=dict(facecolor='white', alpha=0.8), fontsize=9)

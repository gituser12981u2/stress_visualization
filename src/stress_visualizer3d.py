import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider


class StressVisualizer:
    def __init__(self):
        # Initial stress tensor
        self.stress_tensor = np.array([
            [1.0, 0.5, 0.0],  # [σx, τxy, τxz]
            [0.5, 0.0, 0.0],  # [τxy, σy, τyz]
            [0.0, 0.0, 0.0]   # [τxz, τyz, σz]
        ])

        self.fig = plt.figure(figsize=(12, 8))
        self.ax = self.fig.add_subplot(121, projection='3d')
        self.ax2 = self.fig.add_subplot(122)

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

    def get_rotation_matrix(self, theta_x, theta_y, theta_z):
        """Create 3D rotation matrix from Euler angles"""
        theta_x, theta_y, theta_z = np.radians([theta_x, theta_y, theta_z])

        # Rotation matrices for each axis
        Rx = np.array([
            [1, 0, 0],
            [0, np.cos(theta_x), -np.sin(theta_x)],
            [0, np.sin(theta_x), np.cos(theta_x)]
        ])

        Ry = np.array([
            [np.cos(theta_y), 0, np.sin(theta_y)],
            [0, 1, 0],
            [-np.sin(theta_y), 0, np.cos(theta_y)]
        ])

        Rz = np.array([
            [np.cos(theta_z), -np.sin(theta_z), 0],
            [np.sin(theta_z), np.cos(theta_z), 0],
            [0, 0, 1]
        ])

        # Return combined rotation matrix
        return Rz @ Ry @ Rx

    def transform_stress_tensor(self, R):
        """Transform stress tensor using rotation matrix R"""
        return R @ self.stress_tensor @ R.T

    def draw_cube(self):
        theta_x = self.theta_x_slider.val
        theta_y = self.theta_y_slider.val
        theta_z = self.theta_z_slider.val

        self.ax.clear()
        self.ax2.clear()

        R = self.get_rotation_matrix(theta_x, theta_y, theta_z)
        transformed_stress = self.transform_stress_tensor(R)

        # Draw cube and stress visualization
        vertices = np.array([
            [0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
            [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]
        ])

        # Define edges
        edges = [
            [0, 1], [1, 2], [2, 3], [3, 0],
            [4, 5], [5, 6], [6, 7], [7, 4],
            [0, 4], [1, 5], [2, 6], [3, 7]
        ]

        # Rotate vertices
        rotated_vertices = vertices @ R.T

        # Draw edges
        for edge in edges:
            xs = [rotated_vertices[edge[0]][0], rotated_vertices[edge[1]][0]]
            ys = [rotated_vertices[edge[0]][1], rotated_vertices[edge[1]][1]]
            zs = [rotated_vertices[edge[0]][2], rotated_vertices[edge[1]][2]]
            self.ax.plot(xs, ys, zs, 'b')

        # Draw stress arrows
        center = np.mean(rotated_vertices, axis=0)
        scale = 0.3

        # Normal stresses--diagonal elements
        colors = ['r', 'g', 'b']
        for i in range(3):
            if transformed_stress[i, i] != 0:
                direction = np.zeros(3)
                direction[i] = scale * transformed_stress[i, i]
                self.ax.quiver(center[0], center[1], center[2],
                               direction[0], direction[1], direction[2],
                               color=colors[i], label=f'σ{["x", "y", "z"][i]}')

        # Plot stress evolution vs theta-y
        angles = np.linspace(0, 360, 360)
        normal_stresses = []
        shear_stresses = []

        for angle in angles:
            R = self.get_rotation_matrix(theta_x, angle, theta_z)
            transformed = self.transform_stress_tensor(R)
            normal_stresses.append(transformed[0, 0])  # σx
            shear_stresses.append(transformed[0, 1])   # τxy

        self.ax2.plot(angles, normal_stresses, label='σx')
        self.ax2.plot(angles, shear_stresses, label='τxy')
        self.ax2.axvline(x=theta_y, color='k', linestyle='--')
        self.ax2.set_xlabel('θy (degrees)')
        self.ax2.set_ylabel('Stress')
        self.ax2.legend()

        self.ax.set_box_aspect([1, 1, 1])
        self.ax.set_xlabel('X')
        self.ax.set_ylabel('Y')
        self.ax.set_zlabel('Z')
        self.ax.legend()

    def update(self, val):
        self.draw_cube()
        plt.draw()

    def show(self):
        self.draw_cube()
        plt.show()


if __name__ == "__main__":
    visualizer = StressVisualizer()
    visualizer.show()

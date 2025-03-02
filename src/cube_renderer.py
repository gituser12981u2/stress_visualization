import numpy as np
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


class CubeRenderer:
    """Class to handle the 3D rendering of the cube with stress arrows."""

    def __init__(self, ax):
        self.ax = ax
        self.size = 0.5  # Cube half size
        self.arrow_scale = 0.5  # Increased arrow size
        self.threshold = 0.01  # Threshold for displaying stress arrows

        # Define cube vertices (centered at origin)
        self.vertices = np.array([
            [-self.size, -self.size, -self.size],  # 0
            [self.size, -self.size, -self.size],   # 1
            [self.size, self.size, -self.size],    # 2
            [-self.size, self.size, -self.size],   # 3
            [-self.size, -self.size, self.size],   # 4
            [self.size, -self.size, self.size],    # 5
            [self.size, self.size, self.size],     # 6
            [-self.size, self.size, self.size]     # 7
        ])

        # Define faces by vertex indices
        self.faces = [
            [0, 1, 2, 3],  # -z face (bottom)
            [4, 5, 6, 7],  # +z face (top)
            [0, 1, 5, 4],  # -y face (front)
            [3, 2, 6, 7],  # +y face (back)
            [0, 3, 7, 4],  # -x face (left)
            [1, 2, 6, 5]   # +x face (right)
        ]

        self.face_colors = ['lightblue', 'lightblue', 'lightgreen',
                            'lightgreen', 'lightpink', 'lightpink']

    def draw_cube(self, R, transformed_stress):
        """Render the cube with stress arrows given the rotation and stress."""
        self.ax.clear()

        # Rotate vertices
        rotated_vertices = self.vertices @ R.T

        # Create faces collection
        face_collection = Poly3DCollection(
            [rotated_vertices[face] for face in self.faces], alpha=0.3)
        face_collection.set_facecolors(self.face_colors)
        face_collection.set_edgecolor('k')
        self.ax.add_collection3d(face_collection)

        # Define face centers and normals
        face_centers, face_normals = self._calculate_face_properties(
            rotated_vertices, R)

        # Draw normal stress arrows
        self._draw_normal_stress_arrows(
            face_centers, face_normals, transformed_stress)

        # Draw normal stress arrows
        self._draw_shear_stress_arrows(
            face_centers, R, transformed_stress)

        # Add legend
        self.ax.plot([], [], 'r-', label='Normal stress (σ)')
        self.ax.plot([], [], 'b-', label='Shear stress τxy')
        self.ax.plot([], [], 'g-', label='Shear stress τxz')
        self.ax.plot([], [], 'purple', label='Shear stress τyz')

        # Set plot properties
        self._set_plot_properties(rotated_vertices)

    def _calculate_face_properties(self, rotated_vertices, R):
        """Calculate face centers and normals."""
        face_centers = []
        face_normals = []

        # Calculate face centers and normals
        for i, face in enumerate(self.faces):
            face_center = np.mean(rotated_vertices[face], axis=0)
            face_centers.append(face_center)

            # Calculate face normal (for +/- x, y, z faces)
            if i == 0:
                normal = np.array([0, 0, -1])  # -z
            elif i == 1:
                normal = np.array([0, 0, 1])   # +z
            elif i == 2:
                normal = np.array([0, -1, 0])  # -y
            elif i == 3:
                normal = np.array([0, 1, 0])   # +y
            elif i == 4:
                normal = np.array([-1, 0, 0])  # -x
            elif i == 5:
                normal = np.array([1, 0, 0])   # +x

            # Rotate normal
            rotated_normal = normal @ R.T
            face_normals.append(rotated_normal)

        return face_centers, face_normals

    def _draw_normal_stress_arrows(self, face_centers, face_normals,
                                   transformed_stress):
        """Draw normal stress arrows on cube face."""

        # -x and +x faces (σx)
        if abs(transformed_stress[0, 0]) > self.threshold:
            sigma_x = transformed_stress[0, 0]
            # -x face
            self.ax.quiver(face_centers[4][0], face_centers[4][1],
                           face_centers[4][2],
                           -face_normals[4][0] * sigma_x * self.arrow_scale,
                           -face_normals[4][1] * sigma_x * self.arrow_scale,
                           -face_normals[4][2] * sigma_x * self.arrow_scale,
                           color='red', arrow_length_ratio=0.2, linewidth=2)
            # +x face
            self.ax.quiver(face_centers[5][0], face_centers[5][1],
                           face_centers[5][2],
                           face_normals[5][0] * sigma_x * self.arrow_scale,
                           face_normals[5][1] * sigma_x * self.arrow_scale,
                           face_normals[5][2] * sigma_x * self.arrow_scale,
                           color='red', arrow_length_ratio=0.2, linewidth=2)

        # -y and +y faces (σy)
        if abs(transformed_stress[1, 1]) > self.threshold:
            sigma_y = transformed_stress[1, 1]
            # -y face
            self.ax.quiver(face_centers[2][0], face_centers[2][1],
                           face_centers[2][2],
                           -face_normals[2][0] * sigma_y * self.arrow_scale,
                           -face_normals[2][1] * sigma_y * self.arrow_scale,
                           -face_normals[2][2] * sigma_y * self.arrow_scale,
                           color='red', arrow_length_ratio=0.2, linewidth=2)
            # +y face
            self.ax.quiver(face_centers[3][0], face_centers[3][1],
                           face_centers[3][2],
                           face_normals[3][0] * sigma_y * self.arrow_scale,
                           face_normals[3][1] * sigma_y * self.arrow_scale,
                           face_normals[3][2] * sigma_y * self.arrow_scale,
                           color='red', arrow_length_ratio=0.2, linewidth=2)

        # -z and +z faces (σz)
        if abs(transformed_stress[2, 2]) > self.threshold:
            sigma_z = transformed_stress[2, 2]
            # -z face
            self.ax.quiver(face_centers[0][0], face_centers[0][1],
                           face_centers[0][2],
                           -face_normals[0][0] * sigma_z * self.arrow_scale,
                           -face_normals[0][1] * sigma_z * self.arrow_scale,
                           -face_normals[0][2] * sigma_z * self.arrow_scale,
                           color='red', arrow_length_ratio=0.2, linewidth=2)
            # +z face
            self.ax.quiver(face_centers[1][0], face_centers[1][1],
                           face_centers[1][2],
                           face_normals[1][0] * sigma_z * self.arrow_scale,
                           face_normals[1][1] * sigma_z * self.arrow_scale,
                           face_normals[1][2] * sigma_z * self.arrow_scale,
                           color='red', arrow_length_ratio=0.2, linewidth=2)

    def _draw_shear_stress_arrows(self, face_centers, R, transformed_stress):
        """Draw shear stress arrows on cube faces."""

        # Define direction vectors for each face
        x_dir_rotated = np.array([1, 0, 0]) @ R.T
        y_dir_rotated = np.array([0, 1, 0]) @ R.T
        z_dir_rotated = np.array([0, 0, 1]) @ R.T

        # τxy component
        if abs(transformed_stress[0, 1]) > self.threshold:
            tau_xy = transformed_stress[0, 1]

            # On -x face (perpendicular to x, pointing in y direction)
            self.ax.quiver(face_centers[4][0], face_centers[4][1],
                           face_centers[4][2],
                           y_dir_rotated[0] * tau_xy * self.arrow_scale,
                           y_dir_rotated[1] * tau_xy * self.arrow_scale,
                           y_dir_rotated[2] * tau_xy * self.arrow_scale,
                           color='blue', arrow_length_ratio=0.2, linewidth=2)

            # On +x face (perpendicular to x, pointing in -y direction)
            self.ax.quiver(face_centers[5][0], face_centers[5][1],
                           face_centers[5][2],
                           -y_dir_rotated[0] * tau_xy * self.arrow_scale,
                           -y_dir_rotated[1] * tau_xy * self.arrow_scale,
                           -y_dir_rotated[2] * tau_xy * self.arrow_scale,
                           color='blue', arrow_length_ratio=0.2, linewidth=2)

            # On -y face (perpendicular to y, pointing in x direction)
            self.ax.quiver(face_centers[2][0], face_centers[2][1],
                           face_centers[2][2],
                           x_dir_rotated[0] * tau_xy * self.arrow_scale,
                           x_dir_rotated[1] * tau_xy * self.arrow_scale,
                           x_dir_rotated[2] * tau_xy * self.arrow_scale,
                           color='blue', arrow_length_ratio=0.2, linewidth=2)

            # On +y face (perpendicular to y, pointing in -x direction)
            self.ax.quiver(face_centers[3][0], face_centers[3][1],
                           face_centers[3][2],
                           -x_dir_rotated[0] * tau_xy * self.arrow_scale,
                           -x_dir_rotated[1] * tau_xy * self.arrow_scale,
                           -x_dir_rotated[2] * tau_xy * self.arrow_scale,
                           color='blue', arrow_length_ratio=0.2, linewidth=2)

        # τxz component
        if abs(transformed_stress[0, 2]) > self.threshold:
            tau_xz = transformed_stress[0, 2]

            # On -x face (perpendicular to x, pointing in z direction)
            self.ax.quiver(face_centers[4][0], face_centers[4][1],
                           face_centers[4][2],
                           z_dir_rotated[0] * tau_xz * self.arrow_scale,
                           z_dir_rotated[1] * tau_xz * self.arrow_scale,
                           z_dir_rotated[2] * tau_xz * self.arrow_scale,
                           color='green', arrow_length_ratio=0.2, linewidth=2)

            # On +x face (perpendicular to x, pointing in -z direction)
            self.ax.quiver(face_centers[5][0], face_centers[5][1],
                           face_centers[5][2],
                           -z_dir_rotated[0] * tau_xz * self.arrow_scale,
                           -z_dir_rotated[1] * tau_xz * self.arrow_scale,
                           -z_dir_rotated[2] * tau_xz * self.arrow_scale,
                           color='green', arrow_length_ratio=0.2, linewidth=2)

            # On -z face (perpendicular to z, pointing in x direction)
            self.ax.quiver(face_centers[0][0], face_centers[0][1],
                           face_centers[0][2],
                           x_dir_rotated[0] * tau_xz * self.arrow_scale,
                           x_dir_rotated[1] * tau_xz * self.arrow_scale,
                           x_dir_rotated[2] * tau_xz * self.arrow_scale,
                           color='green', arrow_length_ratio=0.2, linewidth=2)

            # On +z face (perpendicular to z, pointing in -x direction)
            self.ax.quiver(face_centers[1][0], face_centers[1][1],
                           face_centers[1][2],
                           -x_dir_rotated[0] * tau_xz * self.arrow_scale,
                           -x_dir_rotated[1] * tau_xz * self.arrow_scale,
                           -x_dir_rotated[2] * tau_xz * self.arrow_scale,
                           color='green', arrow_length_ratio=0.2, linewidth=2)

        # τyz component
        if abs(transformed_stress[1, 2]) > self.threshold:
            tau_yz = transformed_stress[1, 2]

            # On -y face (perpendicular to y, pointing in z direction)
            self.ax.quiver(face_centers[2][0], face_centers[2][1],
                           face_centers[2][2],
                           z_dir_rotated[0] * tau_yz * self.arrow_scale,
                           z_dir_rotated[1] * tau_yz * self.arrow_scale,
                           z_dir_rotated[2] * tau_yz * self.arrow_scale,
                           color='purple', arrow_length_ratio=0.2, linewidth=2)

            # On +y face (perpendicular to y, pointing in -z direction)
            self.ax.quiver(face_centers[3][0], face_centers[3][1],
                           face_centers[3][2],
                           -z_dir_rotated[0] * tau_yz * self.arrow_scale,
                           -z_dir_rotated[1] * tau_yz * self.arrow_scale,
                           -z_dir_rotated[2] * tau_yz * self.arrow_scale,
                           color='purple', arrow_length_ratio=0.2, linewidth=2)

            # On -z face (perpendicular to z, pointing in y direction)
            self.ax.quiver(face_centers[0][0], face_centers[0][1],
                           face_centers[0][2],
                           y_dir_rotated[0] * tau_yz * self.arrow_scale,
                           y_dir_rotated[1] * tau_yz * self.arrow_scale,
                           y_dir_rotated[2] * tau_yz * self.arrow_scale,
                           color='purple', arrow_length_ratio=0.2, linewidth=2)

            # On +z face (perpendicular to z, pointing in -y direction)
            self.ax.quiver(face_centers[1][0], face_centers[1][1],
                           face_centers[1][2],
                           -y_dir_rotated[0] * tau_yz * self.arrow_scale,
                           -y_dir_rotated[1] * tau_yz * self.arrow_scale,
                           -y_dir_rotated[2] * tau_yz * self.arrow_scale,
                           color='purple', arrow_length_ratio=0.2, linewidth=2)

    def _set_plot_properties(self, rotated_vertices):
        """Set properties for the 3D plot."""

        self.ax.set_box_aspect([1, 1, 1])
        self.ax.set_xlabel('X')
        self.ax.set_ylabel('Y')
        self.ax.set_zlabel('Z')
        self.ax.set_title('Stress Tensor Visualization')
        self.ax.legend(loc='upper right')

        # Set limits based on rotated vertices with extra margin
        margin = 0.7
        self.ax.set_xlim([np.min(rotated_vertices[:, 0])*margin,
                         np.max(rotated_vertices[:, 0])*margin])
        self.ax.set_ylim([np.min(rotated_vertices[:, 1])*margin,
                         np.max(rotated_vertices[:, 1])*margin])
        self.ax.set_zlim([np.min(rotated_vertices[:, 2])*margin,
                         np.max(rotated_vertices[:, 2])*margin])

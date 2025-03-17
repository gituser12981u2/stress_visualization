"""
Stress Tensor Module

This module provides the core functionality for stress tensor calculations
and transformations. It handles the mathematical operations needed to transform
a stress tensor under different coordinates rotations.

Classes:
    StressTensor: Manages stress tensor transformations and calculations
"""


import numpy as np


class StressTensor:
    """
    Class to handle stress tensor operations and transformations.

    This class creates and manipulates stress tensors, applying rotational
    transformations and calculating stress variations with angle. It provides
    methods to compute rotation matrices, transform stress tensors, and
    analyze stress components.

    Attributes:
        stress_tensor (numpy.ndarray): 3x3 symmetric stress tensor
    """

    def __init__(self, initial_tensor=None):
        """
        Initialize a stress tensor object.

        Args:
            initial_tensor (numpy.ndarray, optional): Initial 3x3 stress
            tensor.
                If None, a default tensor is used.
        """
        # Initial stress tensor
        if initial_tensor is None:
            self.stress_tensor = np.array([
                [1.0, 0.5, 0.0],  # [σx, τxy, τxz]
                [0.5, 1.0, 0.0],  # [τxy, σy, τyz]
                [0.0, 0.5, 1.0]   # [τxz, τyz, σz]
            ])
        else:
            self.stress_tensor = initial_tensor

    def get_rotation_matrix(self, theta_x, theta_y, theta_z):
        """
        Create 3D rotation matrix from Euler angles

        Args:
            theta_x (float): Rotation angle around x-axis in degrees
            theta_y (float): Rotation angle around y-axis in degrees
            theta_z (float): Rotation angle around z-axis in degrees

        Returns:
            numpy.ndarray: 3x3 rotation matrix (R = Rz*Ry*Rx)
        """
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

    def transform_stress_tensor(self, theta_x, theta_y, theta_z):
        """
        Transform stress tensor using rotation matrix R

        Args:
            theta_x (float): Rotation angle around x-axis in degrees
            theta_y (float): Rotation angle around y-axis  in degrees
            theta_z (float): Rotation angle around z-axis in degrees

        Returns:
            numpy.ndarray: Transformed 3x3 stress tensor
        """
        R = self.get_rotation_matrix(theta_x, theta_y, theta_z)
        return R @ self.stress_tensor @ R.T

    def calculate_stress_variation(self, fixed_theta_x, fixed_theta_z,
                                   angle_samples=360):
        """
        Calculate stress variation for rotation around Y axis

        Args:
            fixed_theta_x (float): Fixed rotation angle around x-axis in
            degrees
            fixed_theta_z (float): Fixed rotation angle around z-axis in
            degrees
            angle_samples (int, optional): Number of angle samples.
            Defaults to 60.

        Returns:
            dict: Dictionary containing stress component arrays for all angles
        """
        angles = np.linspace(0, 360, angle_samples)

        # Initialize arrays to store stress components
        normal_stresses_x = np.zeros(angle_samples)
        normal_stresses_y = np.zeros(angle_samples)
        normal_stresses_z = np.zeros(angle_samples)
        shear_stresses_xy = np.zeros(angle_samples)
        shear_stresses_xz = np.zeros(angle_samples)
        shear_stresses_yz = np.zeros(angle_samples)

        # Calculate stresses at each angle
        for i, angle in enumerate(angles):
            transformed = self.transform_stress_tensor(
                fixed_theta_x, angle, fixed_theta_z)
            normal_stresses_x[i] = transformed[0, 0]  # σx
            normal_stresses_y[i] = transformed[1, 1]  # σy
            normal_stresses_z[i] = transformed[2, 2]  # σz
            shear_stresses_xy[i] = transformed[0, 1]  # τxy
            shear_stresses_xz[i] = transformed[0, 2]  # τxz
            shear_stresses_yz[i] = transformed[1, 2]  # τyz

        return {
            'angles': angles,
            'σx': normal_stresses_x,
            'σy': normal_stresses_y,
            'σz': normal_stresses_z,
            'τxy': shear_stresses_xy,
            'τxz': shear_stresses_xz,
            'τyz': shear_stresses_yz
        }

    def get_stress_text(self, transformed_tensor, theta_x, theta_y, theta_z):
        """
        Generate formatted text for current stress values.

        Args:
            transformed_tensor (numpy.ndarray): Transformed 3x3 stress tensor
            theta_x (float): Rotation angle around x-axis in degrees
            theta_y (float): Rotation angle around y-axis in degrees
            theta_z (float): Rotation angle around z-axis in degrees

        Returns:
            str: Formatted string with stress component values
        """
        return (
            f"Current stresses at θx={theta_x:.1f}°, θy={theta_y:.1f}°, θz={theta_z:.1f}°:\n"
            f"σx = {transformed_tensor[0, 0]:.2f}\n"
            f"σy = {transformed_tensor[1, 1]:.2f}\n"
            f"σz = {transformed_tensor[2, 2]:.2f}\n"
            f"τxy = {transformed_tensor[0, 1]:.2f}\n"
            f"τxz = {transformed_tensor[0, 2]:.2f}\n"
            f"τyz = {transformed_tensor[1, 2]:.2f}"
        )

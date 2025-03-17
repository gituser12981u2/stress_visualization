# Stress Tensor Visualization Tool

This interactive tool provides comprehensive 2D and 3D visualizations for stress tensors, allowing users to understand the complex behavior of stress components under various coordinate transformations.

## Overview

This visualization tool helps in understanding how stress transforms on an element when the coordinate systems are rotated--which is often utilized in material failure analysis, structural design, and mechanics of materials education.

## Features

### 3D Stress Tensor Visualization (`stress_visualizer3d.py`)

- Interactive 3D visualization of stress tensors on a cubic element
- Real-time rotation around three axes ($\theta x$, $\theta y$, $\theta z$) using sliders
- Color-coded stress components:
  - Red arrows: Normal stresses ($\sigma x$, $\sigma y$, $\sigma z$)
  - Blue arrows: Shear stress $\tau xy$
  - Green arrows: Shear stress $\tau xz$
  - Purple arrows: Shear stress $\tau yz$
- Stress variation plots showing how stress components change with rotation angle

### 2D Stress Transformation Visualization (`stress_visualizer2d.py`)

- Interactive 2D visualization of plane stress
- Rotation of the element using a slider
- Color-coded stress components:
  - Red arrows: Normal stress in x-direction ($\sigma x$)
  - Green arrows: Normal stress in y-direction ($\sigma y$)
  - Blue arrows: Shear stress ($\tau xy$)
- Real-time plots of stress variation with rotation angle
- Principal stress and principal angle calculation and visualization

## Project Structure

- `stress_visualizer2d.py`: Main script for 2D stress visualization
- `stress_plotter.py`: Manages 3D plotting of stress variation
- `stress_visualizer3d.py`: Main script for 3D stress visualization
- `cube_renderer.py`: Handles 3D rendering of the cube with stress arrows
- `stress_tensor.py`: Contains the core stress tensor mathematics

## Getting Started

## Prerequisites

- Python 3.x

## Installation

1. Clone this repository
2. Install required packages:

```
pip install -r requirements.txt
```

### Running the Visualizations

To run the 3D visualization:

```
python stress_visualizer3d.py
```

```
python stress_visualizer2d.py
```

## Mathematical background

### Stress tensor

The stress tensor in 3D is represented as a 3x3 symmetric matrix:

$$\sigma = \begin{bmatrix}
\sigma_{xx} & \tau_{xy} & \tau_{xz} \\
\tau_{xy} & \sigma_{yy} & \tau_{yz} \\
\tau_{xz} & \tau_{yz} & \sigma_{zz}
\end{bmatrix}$$

### Stress Transformation

When a coordinate system is rotated, the stress tensor transforms according to:

$$\sigma' = \mathbf{R} \cdot \sigma \cdot \mathbf{R}^T$$

Where $\sigma'$ is the transformed stress tensor, $\sigma$ is the original stress tensor, $\mathbf{R}$ is the rotation matrix, and $\mathbf{R}^T$ is the transpose of the rotation matrix.

### Rotation Matrix

The 3D rotation matrix is created by combining rotations around the x, y, and z axes:

$$\mathbf{R}_x = \begin{bmatrix}
1 & 0 & 0 \\
0 & \cos(\theta_x) & -\sin(\theta_x) \\ 0 & \sin(\theta_x) & \cos(\theta_x) \end{bmatrix}$$

$$\mathbf{R}_y = \begin{bmatrix} \cos(\theta_y) & 0 & \sin(\theta_y) \\ 0 & 1 & 0 \\ -\sin(\theta_y) & 0 & \cos(\theta_y) \end{bmatrix}$$

$$\mathbf{R}_z = \begin{bmatrix} \cos(\theta_z) & -\sin(\theta_z) & 0 \\ \sin(\theta_z) & \cos(\theta_z) & 0 \\ 0 & 0 & 1 \end{bmatrix}$$

$$\mathbf{R} = \mathbf{R}_z \cdot \mathbf{R}_y \cdot \mathbf{R}_x$$

### 2D Stress Transformation

For 2d plane stress, the transformation equations are:

$$\sigma_x' = \frac{\sigma_x + \sigma_y}{2} + \frac{\sigma_x - \sigma_y}{2}\cos(2\theta) + \tau_{xy}\sin(2\theta)$$

$$\sigma_y' = \frac{\sigma_x + \sigma_y}{2} - \frac{\sigma_x - \sigma_y}{2}\cos(2\theta) - \tau_{xy}\sin(2\theta)$$

$$\tau_{xy}' = -\frac{\sigma_x - \sigma_y}{2}\sin(2\theta) + \tau_{xy}\cos(2\theta)$$

### Principal Stresses

Prinicpal stresses are the eigenvalues of the stress tensor, representing the maximum and minimum normal stresses at a point:

$$\sigma_{1, 2} = \frac{\sigma_x + \sigma_y}{2} \pm \sqrt{\left(\frac{\sigma_x - \sigma_y}{2}\right)^2 + \tau_{xy}^2}$$

The principal angle is given by:

$$2\theta_p = \arctan\left(\frac{2\tau_{xy}}{\sigma_x - \sigma_y}\right)$$

## Stress Invariants

The stress tensor has three invariants that remain unchanged under coordinate transformation:

$$I_1 = \sigma_{xx} + \sigma_{yy} + \sigma_{zz} = \text{tr}(\sigma)$$

where $\text{tr}$ is the trace of the matrix.

$$I_2 = \sigma_{xx}\sigma_{yy} + \sigma_{yy}\sigma_{zz} + \sigma_{zz}\sigma_{xx} - \tau_{xy}^2 - \tau_{yz}^2 - \tau_{xz}^2$$

$$I_3 = \det(\sigma)$$

## Pedagogical Applications

This visualization tool has several potential educational benefits:

1. **Spatial Reasoning**: Develops intuition about hoe stress components transform in 3D and 2D elements
2. **Problem Solving**: Allows verification of hand calculations and exploration of transformations on the element.
3. **Interactive Learning**: Encourages active exploration rather than passive learning through immediate visual feedback
4. **Connection to Physical Behavior**: Bridges the gap between mathematical equations and real-world material behavior, such as the principal stresses being $\pm$ 45 degrees from one another.
5. **Principal Stress Visualization**: Demonstrates the significance of principal stresses since it shows how shear stress gradually disappears and the element is moved towards the principal angles.
6. **Stress Trajectory Exploration**: Shows how stress components vary with rotation angle, which helps in understanding optimal material orientation.

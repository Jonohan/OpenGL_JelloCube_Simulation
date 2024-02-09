# Project Overview

This is a computer simulation project aimed at modeling the conversion of external forces and elastic potential energy of a jello cube, achieving physical effects similar to those in the real world.

![Alt text for my GIF](https://github.com/Jonohan/OpenGL_JelloCube_Simulation/blob/main/image/CubeAnime.gif)


**Development Environment:**
- Programmed on Visual Studio 2022 (compatible with other versions).

**Test Files:**
The following is the world file, which needs to be imported to initialize parameters in the simulation environment, such as the initial velocity of the object, the direction of the field, etc.

- Primary world file: `./world/rotate.w` (suitable for testing most functions).
- Alternate test file: `./test.w` (useful for testing mouse dragging without external forces interference).

There is a linked about how to add above world files as command line arguments in Visual Studio:
[condor.depaul.edu/glancast/393class/docs/helpdocs/cmdArgsInVS.html](https://condor.depaul.edu/glancast/393class/docs/helpdocs/cmdArgsInVS.html)

## Physical Simuluaiton Implementations

1. **Structural, Shear, and Bend Springs**
   - These springs are essential components in simulating the physical behavior of materials in a virtual environment. 
	   - Structural springs connect particles in a straight line, typically resembling the material's primary structure. 
	   - Shear springs connect particles diagonally, providing resistance against shearing forces. 
	   - Bend springs connect non-adjacent particles to simulate bending resistance. 
	
	Together, with calculating the elastic and damping forces:
	Hooke's Law (elastic force): $$\mathbf{F}_{\text{Hook}} = -k_{Hook} (\|\mathbf{x}_{A} - \mathbf{x}_{B}\| - RestLen) \frac{\mathbf{x}_{A} - \mathbf{x}_{B}}{\|\mathbf{x}_{A} - \mathbf{x}_{B}\|}$$
	Damping:
	 $$\mathbf{F}_{\text{d}} = -k_{d} \frac{(\mathbf{v}_{A} - \mathbf{v}_{B}) \cdot (\mathbf{x}_{A} - \mathbf{x}_{B})}{\|\mathbf{x}_{A} - \mathbf{x}_{B}\|} \frac{\mathbf{x}_{A} - \mathbf{x}_{B}}{\|\mathbf{x}_{A} - \mathbf{x}_{B}\|}$$
- These springs help in realistically modeling how jello cube deforms and reacts under various forces, such as stretching, compressing, and twisting.

2. **External Forces (Force Field)**
   - This refers to the simulation of forces that originate outside the object's immediate structure. To visually illustrate, below is the formula for linear interpolation in one dimension. 
$$
\mathbf{P} = (1 - a) \mathbf{P}_0 + a \mathbf{P}_1
$$
   - In order to transform equidistant field points in space to each approximate force application point on the jello cube, three-dimensional linear interpolation calculations are used.

   - By applying these forces, the simulation can mimic real-world interactions, such as an object falling due to gravity, swaying due to wind, or being attracted/repelled by magnetic fields. This adds a layer of realism to the simulation, allowing objects to interact with their surroundings in a believable manner.

3. **Bouncing off the Walls**
   - This functionality simulates the interaction of an object with its boundaries or walls, primarily using plane equations to detect collisions. The process involves determining whether a part of the object has exceeded the boundaries defined by these plane equations, indicating a collision with a wall. When such a collision is detected, the response is calculated based on the points of force application on the object and their projection on the plane of the wall. Upon the object entering the wall, a temporary spring force is applied to simulate the effect of bouncing back.
   
   - To accurately simulate the reaction of the object against the inclined plane, different method of collision detection is implemented. 
   - To determine the position of a point relative to this plane, substitute the point's coordinates into the plane equation. Assume a point $P$ with coordinates $(x_0, y_0, z_0)$. Substitute these values into the plane equation as follows: 
   $$F(x_0,y_0,z_0) = ax_0 + by_0 + cz_0 + d $$
	- If $F(x_0,y_0,z_0) > 0$, then point $P$ lies on the side of the plane to which the normal vector points
	- If $F(x_0,y_0,z_0) < 0$, then point $P$ lies on the side of the plane to which the normal vector points
   
4. **Mouse Dragging**: 
   - This feature enables intuitive manipulation of a 3D cube in a simulated environment through mouse interaction. It begins by determining the camera's direction in the world coordinates, which establishes the 'up' and 'right' directions relative to the camera's orientation. This information is crucial for accurately translating the 2D mouse drag movement into a 3D force vector. 
- To calculate the 'up' and 'right' vectors relative to the camera's orientation in a 3D environment, the cross product formula is used:

	1. **Calculate the Right Vector**:
	   The 'right' vector $(\mathbf{R})$ can be obtained by the cross product of the camera's direction vector $(\mathbf{C})$ and the world's 'up' vector $(\mathbf{U}_w = (0, 0, 1))$:
		$$ \mathbf{R} = \mathbf{C} \times \mathbf{U}_w $$
	
	2. **Calculate the Up Vector**:
	   The 'up' vector $(\mathbf{U})$ is then obtained by the cross product of the 'right' vector $(\mathbf{R})$ and the camera's direction vector $(\mathbf{C})$:
		$$ \mathbf{U} = \mathbf{R} \times \mathbf{C} $$
	
	These calculations ensure that the 'right' and 'up' vectors are always perpendicular to the camera's direction and to each other, maintaining orthogonality and consistency with the camera's orientation as it moves and rotates.

   
   - This force vector is then uniformly applied to all nodes of the cube, allowing the user to control and move the entire cube in a manner that mimics physical dragging or pushing. This not only enhances the realism of the interaction but also provides an immersive and intuitive user experience in manipulating the cube within the 3D space.

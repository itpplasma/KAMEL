import numpy as np

class WKB_Grid:
    
    def __init__(self, existing_r, r_res):
        self.existing_r = existing_r
        self.r_min = np.min(existing_r)
        self.r_max = np.max(existing_r)
        self.r_res = r_res

    def generate_non_equidistant_grid(self, start, end, num_points, center, spread):
        """
        Generates a non-uniform grid with higher point density around a given value.

        Parameters:
        - start: The starting value of the grid.
        - end: The ending value of the grid.
        - num_points: The number of points in the grid.
        - center: The value around which points will be more densely clustered.
        - spread: Controls the spread of the points around the center. Smaller values mean more clustering.

        Returns:
        - grid: A numpy array containing the non-uniform grid points.
        """
    
        # Generate a uniform grid in the range [0, 1]
        uniform_grid = np.linspace(0, 1, num_points)
    
        # Apply a transformation to create non-uniform spacing
        # Here we use a sigmoid function to cluster points around the center
        transformed_grid = 1 / (1 + np.exp(-spread * (uniform_grid - 0.5)))
    
        # Normalize the transformed grid to the range [0, 1]
        transformed_grid = (transformed_grid - transformed_grid.min()) / (transformed_grid.max() - transformed_grid.min())
    
        # Scale and translate to the desired range [start, end] and center it around the 'center' value
        grid = start + (end - start) * transformed_grid
    
        # Adjust grid to be centered around the 'center' value
        grid = center - (grid.mean() - grid)
    
        return grid
    
    def gen_grid(self, grid, center, spread, num_points):
            # Convert grid to a numpy array if it's not already
        grid_temp = np.linspace(np.min(grid), np.max(grid), num_points)
        grid_max = np.max(grid)
        
        hrmax = (np.max(grid) - np.min(grid)) / num_points
        new_grid = []
        new_grid.append(grid[0])
        i = 0
        while new_grid[i] < grid_max:
            new_grid.append(new_grid[i] + hrmax / ( 1.0 + 5 * np.exp(-((new_grid[i] - center) / spread)**2)))
            i += 1
    
        return new_grid
    
    def interpolate_on_grid(self):
        pass
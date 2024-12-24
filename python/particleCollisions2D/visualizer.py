import matplotlib.pyplot as plt
import numpy as np

class Visualizer:
    def __init__(self, world, figsize=(8,8)):
        self.world = world
        self.fig, self.ax = plt.subplots(figsize=figsize)
        self.ax.set_aspect('equal')
        self.lines = []          # for bonds
        self.points = None       # for particles
        self.trails = []         # for particle trails
        self._setup_plot()
        
    def _setup_plot(self):
        """Setup plot with initial data"""
        # Plot particles
        positions = []
        for mol in self.world.molecules:
            for p in mol.particles:
                positions.append(p.pos)
        positions = np.array(positions)
        
        self.points = self.ax.scatter(positions[:,0], positions[:,1], c='b', s=100)
        
        # Plot bonds
        for mol in self.world.molecules:
            for bond in mol.bonds:
                p1, p2 = mol.particles[bond.i], mol.particles[bond.j]
                line, = self.ax.plot([p1.pos[0], p2.pos[0]], [p1.pos[1], p2.pos[1]], 'k-')
                self.lines.append(line)
                
        # Set plot limits with some margin
        if len(positions) > 0:
            margin = 2.0
            xmin, xmax = positions[:,0].min(), positions[:,0].max()
            ymin, ymax = positions[:,1].min(), positions[:,1].max()
            dx = max(xmax - xmin, 1.0)
            dy = max(ymax - ymin, 1.0)
            self.ax.set_xlim(xmin - margin*dx, xmax + margin*dx)
            self.ax.set_ylim(ymin - margin*dy, ymax + margin*dy)
        else:
            self.ax.set_xlim(-5, 5)
            self.ax.set_ylim(-5, 5)
            
        plt.grid(True)
        
    def update(self):
        """Update visualization with current particle positions"""
        # Update particles
        positions = []
        for mol in self.world.molecules:
            for p in mol.particles:
                positions.append(p.pos)
        positions = np.array(positions)
        self.points.set_offsets(positions)
        
        # Update bonds
        i = 0
        for mol in self.world.molecules:
            for bond in mol.bonds:
                p1, p2 = mol.particles[bond.i], mol.particles[bond.j]
                self.lines[i].set_data([p1.pos[0], p2.pos[0]], [p1.pos[1], p2.pos[1]])
                i += 1
        
        self.fig.canvas.draw()
        self.fig.canvas.flush_events()
        plt.pause(0.01)  # Small pause to allow GUI to update

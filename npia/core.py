import mahotas as mh
import numpy as np

class AFMImage:
    def __init__(self, image_file, scale_bar=None, ub=None, lb=None):
        self.image = cv2.imread(image_file)
        self.original_image = None
        self.scale_bar = scale_bar
        self.upper_bound = ub
        self.lower_bound = lb
        
        self.n_particles = None
        self.particle_distances = []  # (x, y, radius, nearest_neighbor_dist)
        return

    def prepare(self):
        return

    def apply_threshold(self):
        return

    def watershed(self):
        return

    def get_distances(self):
        return

    def show_image(self):
        return

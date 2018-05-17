import mahotas as mh
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


class Image:
    def __init__(self, image_file, scale_bar=None, ub=None, lb=None):
        self.image = mh.imread(image_file)
        self.original_image = None
        self.scale_bar = scale_bar
        self.upper_bound = ub
        self.lower_bound = lb
        
        self.n_particles = None
        self.particle_distances = []  # (x, y, radius, nearest_neighbor_dist)
        return

    # Dami
    def prepare(self):
        return

    # Luke
    def apply_threshold(self):
        return

    def watershed(self):
        return

    # Alex
    def get_distances(self):
        return

    def show_image(self):
        return

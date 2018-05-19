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

	self.labels = None
	self.pix_to_micron = None
        self.n_particles = None
        self.particle_distances = []  # (x, y, radius, nearest_neighbor_dist)
        return

    # Dami

    def rgb2gray(self):
        return np.dot(self[...,:3], [0.299, 0.587, 0.114])

    def prepare(self, ipx=2015, ipy=2015):		# ipx & ipy are exported pixel counts from AFM
        self.image = rgb2gray(self)			#
        self.image = self.image[2:ipy,2:ipx] 		# Crop image
	scal = self.image[2030,0:2000]			#
	bar = np.count_nonzero(scal-255)		#
	self.pix_to_micron = self.scale_bar/bar 	# pixel to um scale conversion
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

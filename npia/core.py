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
        self.threshold = None
        
        self.pix_to_micron = None
        self.n_particles = None
        self.particle_distances = []  # (x, y, radius, nearest_neighbor_dist)
        return

    # Dami
    def prepare(self):
        return

    # Luke
    def apply_threshold(self, threshold=None):
        """
        Sets threshold for removing background noise.
        
        Parameters
        ----------
        threshold : int
        """
        if threshold:
            self.threshold = threshold
        else:
            # otsu bi-modal histogram thresholding
            # https://en.wikipedia.org/wiki/Otsu%27s_method
            self.threshold = mh.thresholding.otsu(self.image.astype(uint8))

        self.image = self.image[self.image > self.threshold]
        return

    def intensity_hist(self, x_lim=None, y_lim=None):
        plt.hist(pargc.ravel(),bins=int(255/2), fc='k', ec='k')
        plt.xlabel('Intensity')
        plt.ylabel('Pixel Count')

        if x_lim:
            pass
        else:
            plt.xlim(x_lim)

        if y_lim:
            pass
        else:
            plt.ylim(y_lim)

        plt.show()
        return
    
    def apply_gaussian_filter(self, sigma=5):
        """Applys a gaussian filter to image to smooth edges.

        Parameters
        ----------
        sigma : float or int
            Gaussian width"""

        self.image = mh.gaussian_filter(self.image, sigma).astype('uint8')
        return

    def watershed(self):
        Ta = mh.thresholding.otsu(self.image, 0)
        labels, self.n_particles = mh.label(self.image > 0.7 * Ta)
        dist = mh.distance(self.image > 0.05 * self.threshold)
        dist = dist.max() - dist
        dist -= dist.min() # inverting color
        dist = dist / float(dist.ptp()) * 255
        dist = dist.astype(np.uint8)
        dist = mh.stretch(dist, 0, 255)

        thresh = np.median(dist)

        #not accurate to particles detected(?) but matches dist graph well
        thresh = (dist < thresh)  
        areas = 0.9 * mh.cwatershed(dist, labeled)
        self.labels = areas * thresh
        return

    # Alex
    def get_distances(self):
        return

    def show_image(self):
        return

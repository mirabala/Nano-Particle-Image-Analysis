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
        self.threshold = None
        
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
        plt.hist(self.image.ravel(),bins=int(255/2), fc='k', ec='k')
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
    def get_distances(self,n_pix=9):
        Bc1 = np.ones((n_pix,n_pix))
	seeds,n_seeds = mh.label(areas, Bc=Bc1)
	sp = np.where(seeds)
	locg = mh.center_of_mass(self.image, seeds)
	locg = locg.astype(int)
	sg = mh.labeled.labeled_size(seeds)
	particle_radius_mean = np.sqrt(mean(sg[1:])/np.pi)*self.pix_to_micron
	sg = np.sqrt(sg[1:]/np.pi)*self.pix_to_micron
	pr_med = np.median(sg)
	pr_std = np.sqrt(std(sg[1:])/np.pi)*self.pix_to_micron
	
	return pr_med,pr_mean,pr_std,sg
    def Particle_Separation_Analysis
	xsg,ysg = rmaxg.shape
	x = xsg*shape.pix_to_micron # x image length | um
	y = ysg*shape.pix_to_micron # y image length | um
	rho = n_seeds/(y*x) # Particle Density
	xa = np.linspace(0,x,ipx)
	ya = np.linspace(0,y,ipy)
	rmaxp = np.where(rmaxg)
	dmi=np.empty(len(locg[:,0]-1)
	return
    def show_image(self):
        return

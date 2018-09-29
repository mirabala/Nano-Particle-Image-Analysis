import mahotas as mh
import matplotlib.pyplot as plt
import numpy as np
# from numpy.linalg import norm
import pandas as pd  # noqa


def get_distance(x1, y1, x2, y2):
    return np.sqrt((x2 - x1)**2 + (y2 - y1)**2)


class Image:
    def __init__(self, image_file, scale_bar, ub=None, lb=None):
        self.image = mh.imread(image_file)
        self.original_image = None
        self.scale_bar = scale_bar
        self.upper_bound = ub
        self.lower_bound = lb

        self.labels = None
        self.pix_to_micron = None
        self.threshold = None

        self.com = None
        self.sg = None
        self.n_pix = None
        self.ipx = None
        self.ipy = None
        self.locg = None
        self.pix_to_micron = None
        self.n_particles = None
        self.particle_distances = []  # (x, y, radius, nearest_neighbor_dist)
        self.seeds = None
        self.n_seeds = None
        self.spots = None
        
        self.Ta = None
        self.areas = None
        self.dist = None
        self.xloc = []
        self.p_dist = None
        self.da = None
        self.maxes = None
        self.T = None

        return

    # Dami

    def rgb2gray(self):
        """
        Converts the image from rgb to gray scale using the format provided at
        bit.ly/python-rgb-to-gray

        """
        return np.dot(self.image[..., :3], [0.299, 0.587, 0.114])

    def prepare(self, ipx=2015, ipy=2015):
        """
        Crops the image and calculates pixel to micron parameter for later use

        Parameters
        ----------
        ipx : int
            describes the AFM exported image pixel count in the x dimension
        ipy : int
            describes the AFM exported image pixel count in the y dimension
        """
        self.ipx = ipx
        self.ipy = ipy
        self.image = self.rgb2gray()
        scal = self.image[2030, 0:2000]
        self.image = self.image[2:self.ipy, 2:self.ipx]  # Crop image
        self.original_image = self.image
        bar = np.count_nonzero(scal - 255)
        # pixel to um scale conversion
        self.pix_to_micron = self.scale_bar / bar
        return

    # Luke
    def apply_threshold(self, upper_threshold=255, lower_threshold=20):
        """
        Sets threshold for removing background noise.

        Parameters
        ----------
        threshold : int
        """


        self.image[self.image > upper_threshold] = upper_threshold
        self.image[self.image < lower_threshold] = 0

        return
    
    def apply_otsu_threshold(self, s = 0):
        """
        Sets threshold for removing background noise.

        Parameters
        ----------
        threshold : int
        """
#        if threshold:
#            self.threshold = threshold
#        else:
#            # otsu bi-modal histogram thresholding
#            # https://en.wikipedia.org/wiki/Otsu%27s_method
        self.threshold = mh.thresholding.otsu(self.image.astype(np.uint8),s)

        return 
    
    
    def intensity_hist(self, bins=(255/10), x_lim=None, y_lim=None):
        """
        histogram for pixel intensity to help manual thresholding in the removal of background data
        
        Parameters
        ----------
        bins : int
            number of bins to show in histogram
        x_lim : int
            x limit of histogram for purpose of zooming in on a region
        y_lim : int
            y limit of histogram if desired different than largest bin count
        """
        plt.hist(self.image.ravel(), bins=int(bins), fc='k', ec='k')
        plt.xlabel('Intensity')
        plt.ylabel('Pixel Count')

        if x_lim:
            plt.xlim(x_lim)
        else:
            pass

        if y_lim:
            plt.ylim(y_lim)
        else:
            pass

        plt.show()
        return

    def apply_gaussian_filter(self, sigma=5):
        """
        Applys a gaussian filter to image to smooth edges.

        Parameters
        ----------
        sigma : float or int
            Gaussian width
        """

        self.image = mh.gaussian_filter(self.image, sigma).astype('uint8')
        return
    
    def label_particles(self, T = 0, f = 0.7):
        """
        Labels particles by connected pixels 
        
        Parameters
        -----------
        T : int
            threshold to identify particles
        f : int
            modification of threshold value to better capture entire particle region
        """
        #all pixels in "one" particle same value
        labeled,nr_objects = mh.label(self.image > f*T)  
        plt.imshow(labeled,extent=[0,self.pix_to_micron*self.ipx,0,self.pix_to_micron*self.ipy])
        plt.xlabel('Distance / $\mu$m')
        plt.ylabel('Distance / $\mu$m')
        plt.show()

    def get_com(self, scan_filter=(10, 10)):
        """
        Calculates center of mass of particle using regional maxima calculated over the entire matrix
        
        Parameters
        -----------
        scan_filter : int
            size of a weighted square region for regional maxima identification
        """
        self.maxes = mh.regmax(self.image, Bc=np.ones(scan_filter)).astype(int)
        self.spots, n_spots = mh.label(self.maxes, Bc=np.ones(scan_filter))
        com = mh.center_of_mass(self.image, self.spots)
        plt.imshow(self.spots)
        self.com = com
        return

    def watershed(self, Ta = 0):
        """
        Identification of particles through inverted slope comparisons
        
        Parameters
        -----------
        Ta : int
            Threshold value in which the particles will be identified by
        """
        self.Ta=Ta
        dist = mh.distance(self.image > 0.05 * self.Ta)
        dist1 = dist
        dist = dist.max() - dist
        dist -= dist.min()  # inverting color
        dist = dist / float(dist.ptp()) * 255
        dist = dist.astype(np.uint8)
        self.dist = mh.stretch(dist, 0, 255)
        self.labels, self.n_particles = mh.label(self.image > 0.7 * self.Ta)

        thr = np.median(dist)

        # not accurate to particles detected(?)
        # but matches dist graph well
        thresh = (dist < thr)
        areas = 0.9 * mh.cwatershed(dist, self.labels)
        self.areas = areas * thresh
        return

    # Alex
    def get_position(self, n_pix=9):
        """
        Gets particle position and size from watershed analysis

        Parameters
        ----------
        n_pix : float or int
            number of pixels in square array for peak labeling
        """
        
        self.n_pix = n_pix
        Bc1 = np.ones((self.n_pix, self.n_pix))
        self.seeds, self.n_seeds = mh.label(self.areas, Bc=Bc1)
        self.locg = mh.center_of_mass(self.original_image, self.seeds)
        self.locg = self.locg.astype(int)
        sg = mh.labeled.labeled_size(self.seeds)
        particle_radius_mean = np.sqrt(np.mean(sg[1:]) / np.pi) *\
            self.pix_to_micron
        sg = np.sqrt(sg[1:] / np.pi) * self.pix_to_micron
        pr_med = np.median(sg)
        pr_std = np.sqrt(np.std(sg[1:]) / np.pi) * self.pix_to_micron
        self.sg = sg
        return pr_med, particle_radius_mean, pr_std, sg
    
    def get_distances(self,cutoff = 3):
        """
        calculates minimum distance between one particle to all other particles for all particles
        
        Parameters
        -----------
        cutoff : int
            value of normalized particle separation which the desired particles must be greater than
        """
#       Initialize empty array of minimum
#       particle distance for each particle
        self.xloc=np.empty([1,2])
        d=np.empty([len(self.locg[:,0])-1,len(self.locg[:,0])-1]) 
        self.da=np.empty(len(self.locg[:,0])-1)
        a=0
        
        for s in range(0,len(self.locg[:,0])-1):
            
            for ss in range(0,len(self.locg[:,0])-1):
                
#               distance between particle S and every other particle
                d[s,ss]=np.linalg.norm(self.locg[s,:]-self.locg[ss,:])
            
#           distance to closest particle
            self.da[s]=np.amin(d[s,np.nonzero(d[s,:])])
            
            if self.da[s]*self.pix_to_micron/self.sg[s]>cutoff: #thresholded value
                self.xloc[a,:]=self.locg[s,:]
                self.xloc=np.pad(self.xloc,[(0,1),(0,0)],'constant')
                a=a+1

    def Particle_Separation_Analysis(self, cutoff=3):
        """Applys a gaussian filter to image to smooth edges."""
        Bc_ = np.ones((self.n_pix, self.n_pix))
        rmaxg = mh.regmax(self.image, Bc_)
        xsg, ysg = rmaxg.shape
        x = xsg * self.pix_to_micron  # x image length | um
        y = ysg * self.pix_to_micron  # y image length | um
        rho = self.n_seeds / (y * x)  # Particle Density
        
        # l, w = np.shape(xloc)
        # mask = np.ones((l, w), dtype=bool)
        # mask[0:int((l + 1) / 2), 2:w] = False
        # xloc = (np.reshape(xloc[mask, ...], (-1, 2)))
        # xloc = xloc[int(l / 2), :]
        # da_mean = np.mean(da) * self.pix_to_micron
        # da_std = np.std(da) * self.pix_to_micron
        # plt.matshow(d)
        # plt.jet()
        return rho  # , da_mean, da_std

    def new_particle_distances(self, cutoff=3):
        

        for i in range(len(self.com) - 1):
            self.p_dist = get_distance(self.com[i, 0], self.com[i, 1],
                                  self.com[np.arange(len(self.com)) != i, 0],
                                  self.com[np.arange(len(self.com)) != i, 1])
            if np.all(self.p_dist * self.pix_to_micron > cutoff):
                self.xloc.append((self.com[i, 0], self.com[i, 1]))
        return 

    def show_image(self, color_map, x_label='Distance / $\mu$m',
                   y_label='Distance / $\mu$m', overlay_com=None,
                   fig_size=(5, 5),filename= None):
        """Display image input as has been updated

        Parameters
        ----------
        color_map : matplotlib colormap
            Matplotlib colormap designation
        x_label : string
            x label of histogram
        y_label : string
            y label of histogram"""
        plt.figure(figsize=(fig_size))
        plt.imshow(self.image, cmap=color_map,
                   extent=[0, self.pix_to_micron * self.ipx,
                           0, self.pix_to_micron * self.ipy])
        if overlay_com is not None:
            for x, y in overlay_com:
                plt.plot(y*self.pix_to_micron, (self.ipy-x)*self.pix_to_micron, 'o', markerfacecolor='none',
                         markeredgecolor='r', markeredgewidth=2, ms=10)
        else:
            pass
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        if filename is not None:
            plt.savefig(filename)
        plt.show()
        return

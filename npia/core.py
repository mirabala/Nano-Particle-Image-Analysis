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
        self.seeds = None
        return

    # Dami

    def rgb2gray(self):
        return np.dot(self[...,:3], [0.299, 0.587, 0.114])

    def prepare(self, ipx=2015, ipy=2015):		# ipx & ipy are exported pixel counts from AFM
        self.image = rgb2gray(self)
        self.image = self.image[2:ipy,2:ipx] 		# Crop image
	scal = self.image[2030,0:2000]
	bar = np.count_nonzero(scal-255)
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
        """Gets particle position and size from watershed analysis
        
        Parameters
        ----------
        n_pix : float or int
            number of pixels in square array for peak labeling"""
        Bc1 = np.ones((n_pix,n_pix))
        self.seeds,n_seeds = mh.label(areas, Bc=Bc1)
        sp = np.where(seeds)
        locg = mh.center_of_mass(self.image, seeds)
        locg = locg.astype(int)
        sg = mh.labeled.labeled_size(seeds)
        particle_radius_mean = np.sqrt(mean(sg[1:])/np.pi)*self.pix_to_micron
        sg = np.sqrt(sg[1:]/np.pi)*self.pix_to_micron
        pr_med = np.median(sg)
        pr_std = np.sqrt(std(sg[1:])/np.pi)*self.pix_to_micron
        return pr_med,pr_mean,pr_std,sg

    def Particle_Separation_Analysis(self):
        """Calculates closest neighbor to each particle

        Parameters
        ----------
        """
        Bc_=np.ones((n_pix,n_pix))
        rmaxg = mh.regmax(pargcfs,Bc_)
        xsg,ysg = rmaxg.shape
        x = xsg*shape.pix_to_micron # x image length | um
        y = ysg*shape.pix_to_micron # y image length | um
        rho = n_seeds/(y*x) # Particle Density
        xa = np.linspace(0,x,ipx)
        ya = np.linspace(0,y,ipy)
        rmaxp = np.where(rmaxg)
        dmi=np.empty(len(locg[:,0])-1) #initialize empty array of minimum particle distance for each particle
        d=np.empty([len(locg[:,0])-1,len(locg[:,0])-1])
        da=np.empty(len(locg[:,0])-1)
        a=0
        xloc=np.empty([1,2])
        for s in range(0,len(locg[:,0])-1):
            dm=np.empty(len(locg[:,0])-1)   #distance between particle S and every other particle
            for ss in range(0,len(locg[:,0])-1):
                d[s,ss]=norm(locg[s,:]-locg[ss,:])
                dm[ss]=np.sqrt(np.square(xa[locg[s,0]]-xa[locg[ss,0]])+np.square(ya[locg[s,1]]-ya[locg[ss,1]]))
            dmi[s]=np.amin(dm[np.nonzero(dm)])
            da[s]=np.amin(d[s,np.nonzero(d[s,:])])
            if da[s]*self.pix_to_micron/sg[s]>3: #thresholded value
                xloc[a,[0,1]]=locg[s,:]
                xloc=np.pad(xloc,[0,2],'constant',constant_values=[0])
                a=a+1
        print(xloc)
        l,w=np.shape(xloc)
        mask=np.ones((l,w),dtype=bool)
        mask[0:int((l+1)/2),2:w]=False
        xloc=(np.reshape(xloc[mask,...],(-1,2)))
        xloc=xloc[int(l/2),:]
        da_mean = np.mean(da)*self.pix_to_micron
        da_std = np.std(da)*self.pix_to_micron
        matshow(d);
        jet()
        print(xloc)
        return da_mean, da_std

    def show_image(self,self.seeds,color_map,x_label='Distance / $\mu$m',y_label='Distance / $\mu$m'):
        """Display image input as has been updated

        Parameters
        ----------
        color_map : matplotlib colormap
            Matplotlib colormap designation
        x_label : string
            x label of histogram
        y_label : string
            y label of histogram"""
        imshow(self.seeds, cmap=color_map, extent=[0,self.pix_to_micron*ipx,0,self.pix_to_micron*ipy])
        xlabel('x_label')
        ylabel('y_label)
        show()
        return

    def show_hist(self, data, n_bins=255/10,x_label='Intensity',y_label='Pixel Count'):
        """Shows histogram of image pixel intensity

        Parameters
        ----------
        n_bins : float or int
            number of histogram bins
        x_label : string
            x label of histogram
        y_label : string
            y label of histogram"""
        h=hist(data.ravel(),bins=n_bins,fc='k', ec='k')
        xlabel(x_label)
        ylabel(y_label)
	show()
        return

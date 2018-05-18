import mahotas as mh
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd  # noqa


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

        self.sg
        self.n_pix = None
        self.ipx = None
        self.ipy = None
        self.locg = None
        self.pix_to_micron = None
        self.n_particles = None
        self.particle_distances = []  # (x, y, radius, nearest_neighbor_dist)
        self.seeds = None
        self.n_seeds = None
        return

    # Dami

    def rgb2gray(self):
        return np.dot(self.image[..., :3], [0.299, 0.587, 0.114])

    def prepare(self, ipx=2015, ipy=2015):
        self.ipx = ipx
        self.ipy = ipy
        self.image = self.rgb2gray()
        self.image = self.image[2:self.ipy, 2:self.ipx]  # Crop image
        scal = self.image[2030, 0:2000]
        bar = np.count_nonzero(scal - 255)
        # pixel to um scale conversion
        self.pix_to_micron = self.scale_bar / bar
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
            self.threshold = mh.thresholding.otsu(self.image.astype(np.uint8))

        self.image = self.image[self.image > self.threshold]
        return

    def intensity_hist(self, bins=(255/10), x_lim=None, y_lim=None):
        plt.hist(self.image.ravel(), bins=int(bins), fc='k', ec='k')
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
        self.labels, self.n_particles = mh.label(self.image > 0.7 * Ta)
        dist = mh.distance(self.image > 0.05 * self.threshold)
        dist = dist.max() - dist
        dist -= dist.min()  # inverting color
        dist = dist / float(dist.ptp()) * 255
        dist = dist.astype(np.uint8)
        dist = mh.stretch(dist, 0, 255)

        thresh = np.median(dist)

        # not accurate to particles detected(?)
        # but matches dist graph well
        thresh = (dist < thresh)
        areas = 0.9 * mh.cwatershed(dist, self.labels)
        self.areas = areas * thresh
        return

    # Alex
    def get_distances(self, n_pix=9):
        """Gets particle position and size from watershed analysis

        Parameters
        ----------
        n_pix : float or int
            number of pixels in square array for peak labeling"""
        self.n_pix = n_pix
        Bc1 = np.ones((self.n_pix, self.n_pix))
        self.seeds, self.n_seeds = mh.label(self.areas, Bc=Bc1)
        self.locg = mh.center_of_mass(self.image, self.seeds)
        self.locg = self.locg.astype(int)
        sg = mh.labeled.labeled_size(self.seeds)
        particle_radius_mean = np.sqrt(np.mean(sg[1:]) / np.pi) *\
            self.pix_to_micron
        sg = np.sqrt(sg[1:] / np.pi) * self.pix_to_micron
        pr_med = np.median(sg)
        pr_std = np.sqrt(np.std(sg[1:]) / np.pi) * self.pix_to_micron
        self.sg = sg
        return pr_med, particle_radius_mean, pr_std, sg

    def Particle_Separation_Analysis(self):
        """Applys a gaussian filter to image to smooth edges."""
        Bc_ = np.ones((self.n_pix, self.n_pix))
        rmaxg = mh.regmax(self.image, Bc_)
        xsg, ysg = rmaxg.shape
        x = xsg * self.pix_to_micron  # x image length | um
        y = ysg * self.pix_to_micron  # y image length | um
        rho = self.n_seeds / (y * x)  # Particle Density
        xa = np.linspace(0, x, self.ipx)
        ya = np.linspace(0, y, self.ipy)
        # Initialize empty array of minimum
        # particle distance for each particle
        dmi = np.empty(len(self.locg[:, 0]) - 1)
        d = np.empty([len(self.locg[:, 0]) - 1, len(self.locg[:, 0]) - 1])
        da = np.empty(len(self.locg[:, 0]) - 1)
        a = 0
        xloc = np.empty([1, 2])
        for s in range(0, len(self.locg[:, 0]) - 1):
            # distance between particle S and every other particle
            dm = np.empty(len(self.locg[:, 0]) - 1)
            for ss in range(0, len(self.locg[:, 0]) - 1):
                d[s, ss] = np.norm(self.locg[s, :] - self.locg[ss, :])
                dm[ss] = np.sqrt(np.square(xa[self.locg[s, 0]] -
                                           xa[self.locg[ss, 0]])
                                 + np.square(ya[self.locg[s, 1]] -
                                             ya[self.locg[ss, 1]]))
            dmi[s] = np.amin(dm[np.nonzero(dm)])
            da[s] = np.amin(d[s, np.nonzero(d[s, :])])
            # thresholded value
            if da[s] * self.pix_to_micron / self.sg[s] > 3:
                xloc[a, [0, 1]] = self.locg[s, :]
                xloc = np.pad(xloc, [0, 2], 'constant', constant_values=[0])
                a = a + 1
        print(xloc)
        l, w = np.shape(xloc)
        mask = np.ones((l, w), dtype=bool)
        mask[0:int((l + 1) / 2), 2:w] = False
        xloc = (np.reshape(xloc[mask, ...], (-1, 2)))
        xloc = xloc[int(l / 2), :]
        da_mean = np.mean(da) * self.pix_to_micron
        da_std = np.std(da) * self.pix_to_micron
        plt.matshow(d)
        plt.jet()
        print(xloc)
        return da_mean, da_std, rho

    def show_image(self, color_map, x_label='Distance / $\mu$m',
                   y_label='Distance / $\mu$m'):
        """Display image input as has been updated

        Parameters
        ----------
        color_map : matplotlib colormap
            Matplotlib colormap designation
        x_label : string
            x label of histogram
        y_label : string
            y label of histogram"""
        plt.imshow(self.seeds, cmap=color_map,
                   extent=[0, self.pix_to_micron * self.ipx,
                           0, self.pix_to_micron * self.ipy])
        plt.xlabel('x_label')
        plt.ylabel('y_label')
        plt.show()
        return

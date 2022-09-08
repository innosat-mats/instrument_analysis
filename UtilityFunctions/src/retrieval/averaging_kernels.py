#%%
import numpy as np
import plotly.graph_objects as go
from matplotlib import pyplot as plt
from joblib import Parallel, delayed

#%%
class Grid:
    """
    A class to represent an n-dimensional (n<4) grid.

    ...

    Attributes
    ----------
    dimension : int
        dimensionality of grid
    x : ndarray
        x-grid
    dx : float
        grid distance
    y  : ndarray
        y-grid
    dy : float
        grid distance
    z  : ndarray 
        z-grid
    dz : float
        grid distance

    Methods
    -------
    meshgrid():
        Returns the meshgrid of suitable dimension for the Grid
    """

    def __init__(self,x,y=np.array([]),z=np.array([])):
        """
        Constructs a Grid object

        Parameters
        ----------
        x : ndarray
            x-grid 
        y : ndarray, optional
            y-grid 
        z : ndarray, optional
            z-grid 
        """
        
        if (y.size == 0) & (z.size == 0):
            self.dimension = 1 
            self.x = x
            self.dx = x[1]-x[0]
        elif z.size == 0:
            self.dimension = 2
            self.x = x
            self.y = y
            self.dx = x[1]-x[0]
            self.dy = y[1]-y[0]
        else:
            self.dimension = 3
            self.x = x
            self.y = y
            self.z = z
            self.dx = x[1]-x[0]
            self.dy = y[1]-y[0]
            self.dz = z[1]-z[0]


    def meshgrid(self):
        """
        Function for returning the meshgrid of the Grid object

        Parameters
        ----------
        None

        Returns
        ----------
        XX, YY, ZZ: ndarray
            Returns a meshgrid of same dimensionality of the Grid in order x,y,z (i.e. indexing='ij')
        """
        
        if self.dimension == 1:
            return self.x
        elif self.dimension == 2:
            return np.meshgrid(self.x,self.y,indexing='ij')
        elif self.dimension == 3:
            return np.meshgrid(self.x,self.y,self.z,indexing='ij')
        else:
            raise ValueError

class Kernel:
    """
    A class to represent an n-dimensional (n<4) averaging kernel.

    ...

    Attributes
    ----------
    dimensions : int
        dimensionality of AVK
    sigma : ndarray of len(dimension)
        sigma in each dimension the gaussian AVK case

    Methods
    -------
    plot_kernel():
        Plots the averaging kernel
    
    get_kernel(grid: Grid, x0: ndarray):
        Plots the averaging kernel at position x0
    """

    def __init__(self,widths,type='gaussian'):
        """
        Constructs a Kernel object

        Parameters
        ----------
        widths : ndarray
            FWHM of averaging kernels defined in an nd array equal to the Kernel dimensions
        type : str
            functional form of the AVK (currently only 'gaussan' implemented)
        """

        if type != 'gaussian':
            raise NotImplementedError('Only gaussian kernels are implemented')
        
        if type == 'gaussian':
            self.dimensions = len(widths)
            self.sigma = np.array(widths)/np.sqrt(2*np.log2(2))
            if np.any(widths <= 0):
                raise NotImplementedError('AVK must have finite width in all dimensions')

    def plot_kernel(self,plot=True):
        """
        Function for plotting (or returning) the AVK on a suitable grid

        Parameters
        ----------
        plot: bool, optional
            Whether to plot the kernel or not (default is True)

        Returns
        ----------
        grid: Grid
            The grid which the kernel is plotted on
        
        AVK: nparray
            The AVK on the plotted grid
        """
        
        gridwidth=10
        if self.dimensions==1:
            grid = Grid(np.arange(-gridwidth*self.sigma[0],gridwidth*self.sigma[0],0.1*self.sigma[0]))
            AVK = self.get_kernel(grid,np.array([0]))
            if plot:
                plt.plot(grid.x,AVK)

        elif self.dimensions==2:
            grid = Grid(np.arange(-gridwidth*self.sigma[0],gridwidth*self.sigma[0],0.1*self.sigma[0]),np.arange(-gridwidth*self.sigma[1],gridwidth*self.sigma[1],0.1*self.sigma[1]))
            AVK = self.get_kernel(grid,np.array([0,0]))
            if plot:
                plt.pcolor(grid.x,grid.y,AVK.T)
                plt.gca().set_aspect('equal', adjustable='box')

        elif self.dimensions==3:
            gridwidth=5
            grid = Grid(np.arange(-gridwidth*self.sigma[0],gridwidth*self.sigma[0],0.2*self.sigma[0]),np.arange(-gridwidth*self.sigma[1],gridwidth*self.sigma[1],0.2*self.sigma[1]),np.arange(-gridwidth*self.sigma[2],gridwidth*self.sigma[2],0.2*self.sigma[2]))
            AVK = self.get_kernel(grid,np.array([0,0,0]))
            AVK_peak = AVK.max()

            fig = go.Figure(data=go.Isosurface(
            x=grid.meshgrid()[0].flatten(),
            y=grid.meshgrid()[1].flatten(),
            z=grid.meshgrid()[2].flatten(),
            value=AVK.flatten()/AVK_peak,
            opacity=0.6,
            isomin=1/4,
            isomax=1,
            surface_count=3))
            fig.show()


        else:
            raise ValueError('Dimensions cannot be other than 1,2,3')

        return grid,AVK

    def get_kernel(self,grid,x0):
        """
        Returns the AVK for given grid and position

        Parameters
        ----------
        grid: Grid
            Grid to return the AVK on

        x0: nparray of length dimension
            Position of kernel in the grid (value of row index of the AVK matrix)

        Returns
        ----------
        
        AVK: nparray
            The AVK on the given grid
        """
        
        if self.dimensions != grid.dimension:
            raise NotImplementedError('Grid needs to be same dimension as AVK')

        if grid.dimension == 1:
            AVK =  np.exp(-(((grid.meshgrid()-x0[0])**2/(2*self.sigma[0]**2)))) 
            C = 1/(grid.dx*np.sum(AVK))
            return C*AVK
        if grid.dimension == 2:
            AVK =  np.exp(-(((grid.meshgrid()[0]-x0[0])**2/(2*self.sigma[0]**2))+((grid.meshgrid()[1]-x0[1])**2/(2*self.sigma[1]**2))))  
            C = 1/(grid.dx*grid.dy*np.sum(np.sum(AVK)))
            return C*AVK
        if grid.dimension == 3:
            AVK =  np.exp(-(((grid.meshgrid()[0]-x0[0])**2/(2*self.sigma[0]**2))+((grid.meshgrid()[1]-x0[1])**2/(2*self.sigma[1]**2))+((grid.meshgrid()[2]-x0[2])**2/(2*self.sigma[2]**2))))  
            C = 1/(grid.dx*grid.dy*grid.dz*np.sum(np.sum(np.sum(AVK))))
            return C*AVK
        else:
            raise NotImplementedError   

    def apply_kernel(self,field,grid,index):
        """
        Smoothes a field on a given grid with the AVK applied at an index

        Parameters
        ----------
        field: nparray
            the field to apply the AVK to, ndarray of correct dimensionality

        grid: Grid
            the grid of the field

        index: ndarray
            grid index of the AVK to apply

        Returns
        ----------
        
        field_smooth: nparray
            The smoothed field
        """
        
        if grid.dimension != self.dimensions:
            raise NotImplementedError('field and AVK needs to be same dimension')

        if grid.dimension == 1:
            x0 = np.array([grid.x[index]])
            AVK = self.get_kernel(grid,x0)
            field_smooth = np.sum(field*AVK*grid.dx)

            return field_smooth

        if grid.dimension == 2:
            x0 = np.array([grid.x[index[0]],grid.y[index[1]]])
            AVK_analytical = self.get_kernel(grid,x0)
            field_smooth = np.sum(np.sum(field*AVK_analytical*grid.dx*grid.dy))

            return field_smooth

        if grid.dimension == 3:
            x0 = np.array([grid.x[index[0]],grid.y[index[1]],grid.z[index[2]]])
            AVK_analytical = self.get_kernel(grid,x0)
            field_smooth = np.sum(np.sum(np.sum(field*AVK_analytical*grid.dx*grid.dy*grid.dz)))

            return field_smooth
        else:
            raise NotImplementedError


def apply_3d_kernel(field, x, y, z, fwhm, only_kernel=True,
                    parallel=True, n_jobs=16):
    """
    Applies averaging kernels to a 3D field

    Parameters
    ----------
    field : array-like
        three dimensional field to apply the AVK to
    x : array
        coordinates in first dim
    y : array
        coordinates in second dim
    z : array
        coordinates in second dim
    fwhm : array
        full width half mean of kernels [fwhm_x, fwhm_y, fwhm_z)
    only_kernel : bool
        return only kernel matrix (for debugging)
    parallel : bool
        if parallel processing should be used
    n_jobs : int
        number of parallel processes

    Returns
    ----------

    averaged_field: array
        3d field with applied kernels (or AVK on grid if only_kernel == True)

    """

    def parallel_loop(avg_field, field, grid_3, i, ny, nz):
        for j in np.arange(ny):
            for k in np.arange(nz):
                avg_field[j, k] = AVK.apply_kernel(field, grid_3,
                                                   np.array([i, j, k]))
        return avg_field

    grid_3 = Grid(x, y, z)
    avg_field = np.zeros((len(y), len(z)))
    AVK = Kernel(np.array([fwhm[0], fwhm[1], fwhm[2]]))

    if only_kernel:
        AVK.plot_kernel()

    else:
        if parallel:
            # parallel processing (watch max_nbytes / n_jobs)
            avg_field = (Parallel(n_jobs=n_jobs, max_nbytes=None)
                                 (delayed(parallel_loop)(avg_field,
                                  field, grid_3, i, len(y), len(z))
                                  for i in range(0, len(x))))

        else:
            # no parallel processing (lengthy!)
            for i in np.arange(len(x)):
                for j in np.arange(len(y)):
                    for k in np.arange(len(z)):
                        avg_field[i, j, k] = AVK.apply_kernel(field, grid_3, np.array([i, j, k]))

    return avg_field


# %%

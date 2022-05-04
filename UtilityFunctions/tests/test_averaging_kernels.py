import numpy as np
from matplotlib import pyplot as plt
from retrieval.averaging_kernels import Kernel, Grid

def test_1d_avk():

    #%%
    # #1D example
    x = np.arange(0,100,0.1)
    y = np.sin(x*1.5)
    grid_1 = Grid(x)

    AVK = Kernel(np.array([2]))
    AVK.plot_kernel()

    #%%
    y_smooth=np.zeros(len(x))
    for index in np.arange(len(x)):
        y_smooth[index] = AVK.apply_kernel(y,grid_1,index)

    plt.plot(x,y,x,y_smooth)
    plt.show()

    return True

def test_2d_avk():
        
    #%%
    #2D example with constant AVK. Signal in x direction
    x = np.arange(0,40,0.2)
    y = np.arange(0,80,0.2)
    [xx,yy] = np.meshgrid(x,y,indexing='ij') #use ij indexing to have x along first axis
    z = np.sin(xx/2)

    grid_2 = Grid(x,y)

    AVK = Kernel(np.array([2,8]))
    AVK.plot_kernel()

    #%%
    z_smooth=np.zeros((len(x),len(y)))
    for i in np.arange(len(x)):
        for j in np.arange(len(y)):
            z_smooth[i,j] = AVK.apply_kernel(z,grid_2,np.array([i,j]))

    plt.pcolor(x,y,z.T)
    plt.title('Orginal field')
    plt.clim([-1,1])
    plt.show()
    plt.pcolor(x,y,z_smooth.T)
    plt.title('Smoothed field')
    plt.clim([-1,1])


    #%%
    #2D example with constant AVK. Signal in y direction


    x = np.arange(0,40,0.2)
    y = np.arange(0,80,0.2)
    [xx,yy] = np.meshgrid(x,y,indexing='ij') #use ij indexing to have x along first axis
    z = np.sin(yy/2)

    grid_2 = Grid(x,y)

    AVK = Kernel(np.array([2,8]))
    AVK.plot_kernel()

    #%%
    z_smooth=np.zeros((len(x),len(y)))
    for i in np.arange(len(x)):
        for j in np.arange(len(y)):
            z_smooth[i,j] = AVK.apply_kernel(z,grid_2,np.array([i,j]))

    plt.pcolor(x,y,z.T)
    plt.title('Orginal field')
    plt.clim([-1,1])
    plt.show()
    plt.pcolor(x,y,z_smooth.T)
    plt.title('Smoothed field')
    plt.clim([-1,1])


    #%%
    #2D example with varying AVK

    x = np.arange(0,40,0.2)
    y = np.arange(0,80,0.2)
    [xx,yy] = np.meshgrid(x,y,indexing='ij') #use ij indexing to have x along first axis
    z = np.sin(xx/2)

    grid_2 = Grid(x,y)

    #%%
    z_smooth=np.zeros((len(x),len(y)))
    for i in np.arange(len(x)):
        for j in np.arange(len(y)):
            AVK = Kernel(np.array([max(0.01,y[j]/10),2]))
            z_smooth[i,j] = AVK.apply_kernel(z,grid_2,np.array([i,j]))

    plt.pcolor(x,y,z.T)
    plt.title('Orginal field')
    plt.clim([-1,1])
    plt.show()
    plt.pcolor(x,y,z_smooth.T)
    plt.title('Smoothed field')
    plt.clim([-1,1])

    return True


def test_3d_AVK():

    #%%
    #3D example with constant AVK. Signal in x direction

    x = np.arange(0,20,1)
    y = np.arange(0,30,0.5)
    z = np.arange(0,40,0.2)

    [xx,yy,zz] = np.meshgrid(x,y,z,indexing='ij') #use ij indexing to have x along first axis
    grid_3 = Grid(x,y,z)

    AVK = Kernel(np.array([2,4,8]))
    AVK.plot_kernel()

    #%%
    x = np.arange(0,20,1)
    y = np.arange(0,30,0.5)
    z = np.arange(0,40,0.2)
    [xx,yy,zz] = np.meshgrid(x,y,z,indexing='ij') #use ij indexing to have x along first axis
    grid_3 = Grid(x,y,z)

    AVK = Kernel(np.array([2,4,8]))


    field = np.sin(xx/2)
    field_smooth=np.zeros((len(x),len(y),len(z)))
    for i in np.arange(len(x)):
        for j in np.arange(len(y)):
            for k in np.arange(len(z)):
                field_smooth[i,j,k] = AVK.apply_kernel(field,grid_3,np.array([i,j,k]))

    #%%
    z_sel = int(len(z)/2)
    plt.pcolor(x,y,field[:,:,z_sel].T)
    plt.title('Orginal field')
    plt.clim([-1,1])
    plt.show()
    plt.pcolor(x,y,field_smooth[:,:,z_sel].T)
    plt.title('Smoothed field')
    plt.clim([-1,1])

    #%%
    #3D example with constant AVK. Signal in y direction

    field = np.sin(yy/2)
    field_smooth=np.zeros((len(x),len(y),len(z)))
    for i in np.arange(len(x)):
        for j in np.arange(len(y)):
            for k in np.arange(len(z)):
                field_smooth[i,j,k] = AVK.apply_kernel(field,grid_3,np.array([i,j,k]))

    #%%
    z_sel = int(len(z)/2)
    plt.pcolor(x,y,field[:,:,z_sel].T)
    plt.title('Orginal field')
    plt.clim([-1,1])
    plt.show()
    plt.pcolor(x,y,field_smooth[:,:,z_sel].T)
    plt.title('Smoothed field')
    plt.clim([-1,1])


    #%%
    #3D example with constant AVK. Signal in z direction

    field = np.sin(zz/2)
    field_smooth=np.zeros((len(x),len(y),len(z)))
    for i in np.arange(len(x)):
        for j in np.arange(len(y)):
            for k in np.arange(len(z)):
                field_smooth[i,j,k] = AVK.apply_kernel(field,grid_3,np.array([i,j,k]))

    #%%
    x_sel = int(len(x)/2)
    plt.pcolor(y,z,field[x_sel,:,:].T)
    plt.title('Orginal field')
    plt.clim([-1,1])
    plt.show()
    plt.pcolor(y,z,field_smooth[x_sel,:,:].T)
    plt.title('Smoothed field')
    plt.clim([-1,1])

    return True

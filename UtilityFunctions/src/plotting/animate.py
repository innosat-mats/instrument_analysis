import imageio
import glob
import os


def generate_gif(directory, outfile):
    """Generate gif from CCD image figures.
    For now image names are assumed to end with indeces.
    See plotCCD scripts for naming convention.

    Parameters
    ----------
    directory : str
        input directory
    outfile : str
        full path to output file (.gif)
    """

    filenames = glob.glob(f"{directory}/*")

    filenames = sorted(filenames, key=lambda x:
                       int(os.path.split(x)[1].split('_')[-1].split('.')[0]))

    with imageio.get_writer(outfile, mode='I') as writer:
        for filename in filenames:
            print(filename)
            image = imageio.imread(filename)
            writer.append_data(image)

    return

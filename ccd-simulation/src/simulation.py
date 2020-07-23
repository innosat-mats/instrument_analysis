import numpy as np


class CCD:
    """CCD class"""

    def __init__(self, rows, columns):
        self.rows = rows
        self.columns = columns
        self.image = np.zeros([rows, columns])
        self.shift_register = np.zeros([1, columns])
        self.well = 0

    def size(self):
        return self.rows, self.columns

    def clock_row(self):
        self.shift_register = self.shift_register + self.image[0, :]
        self.image[0:-1, :] = self.image[1:, :]
        self.image[-1, :].fill(0)

    def clock_column(self):
        self.well = self.well + self.shift_register[0, 0]
        self.shift_register[0, 0:-1] = self.shift_register[0, 1:]
        self.shift_register[0, -1] = 0

    def reset_shift_register(self):
        self.shift_register.fill(0)

    def read_well(self):
        value = self.well
        self.well = 0
        return value

    def get_image(
        self, nrow=None, nrskip=0, nrbin=1, ncol=None, ncskip=0, ncbin=1,
    ):

        if nrow is None:
            nrow = self.rows
        if ncol is None:
            ncol = self.columns

        out_image = np.zeros([nrow, ncol])

        for i in range(nrow):
            self.reset_shift_register()

            for _ in range(nrbin):
                self.clock_row()

            for j in range(ncol):

                for _ in range(ncbin):
                    self.clock_column()

                value = self.read_well()
                out_image[i, j] = value

        return out_image

    def set_image(self, image):
        if image.shape != (self.rows, self.columns):
            raise ValueError("Image size must mach CCD size")
        else:
            self.image = image

    def reset_ccd(self):
        self.image = np.zeros([self.rows, self.columns])
        self.shift_register = np.zeros([1, self.columns])
        self.well = 0

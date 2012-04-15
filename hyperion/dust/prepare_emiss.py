from __future__ import print_function, division

import atpy
import numpy as np

from ..util.constants import c

variable = {}
variable['T'] = 'temperature'
variable['E'] = 'specific_energy'


def prepare_emiss(files, values, emissvar, filename_out):
    '''
    Prepare an emissivity FITS file to read in to prepare_dust

    Required arguments:

        *files*: [ list of strings ]
            List of two column files containing wavelength and emissivity. All
            files should have the same number of lines and be defined on the
            same wavelength grid. The wavelength should be given in microns,
            and the emissivity should be in cgs.

        *values*: [ list of floats ]
            List of the independent variable values corresponding to files

        *emissvar*: [ string ]
            Independent variable name. This can be:
                - temperature
                - specific_energy

        *filename_out*: [ string ]
            Filename for the output file
    '''

    # Create empty emissivity array
    n_j = len(files)
    n_lines = file(files[0]).read().count('\n')
    jnu = np.zeros((n_lines, n_j))

    for i, filename in enumerate(files):

        # Read in emissivity file
        data = np.loadtxt(filename, dtype=[('wav', float), ('jnu', float)])

        # Check the file length
        if len(data) != n_lines:
            raise Exception("Files with different lengths")

        # Save wavelength grid
        if i == 0:
            wav = data['wav']
        else:
            if np.any(data['wav'] != wav):
                raise Exception("Wavelength grids do not agree")

        # Save the emissivity
        jnu[:, i] = data['jnu']

    # Order the emissivities by the independent variable
    order = np.argsort(values)
    jnu = jnu[:, order]

    # Create a table set
    ts = atpy.TableSet()

    # Create a table with frequency and emissivity
    t = atpy.Table(name='Emissivities')
    t.add_column('nu', c * 1.e4 / wav)
    t.add_column('jnu', jnu)
    t.sort('nu')
    ts.append(t)

    # Create a table with the independent variable
    t = atpy.Table(name='Emissivity variable')
    t.add_column(variable[emissvar], values)
    ts.append(t)

    # Write emissivities to disk
    ts.write(filename_out, overwrite=True)

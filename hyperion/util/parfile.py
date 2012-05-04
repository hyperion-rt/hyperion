from __future__ import print_function, division

import os


def smart_type(value):
    '''
    Return value as an int or float if possible
    '''

    try:
        value = int(value)
        return value
    except:
        try:
            value = float(value)
            return value
        except:
            value = value.replace("'", "")
            if value.lower() in ['yes', 'no']:
                return value.lower() == 'yes'
            else:
                return value


def parse(filename):
    '''
    Parse a parfile with each line following the format:

    value = key = comment
    '''

    if not os.path.exists(filename):
        raise Exception('No such file or directory: ' + filename)

    f = open(filename, 'rb')

    parameters = {}

    for line in f.readlines():
        if '=' in line:
            cols = line.split('=')
            value, key = cols[0].strip(), cols[1].strip()
            parameters[key.lower()] = smart_type(value)

    f.close()

    return parameters

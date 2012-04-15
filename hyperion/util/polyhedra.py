from __future__ import print_function, division

from math import sqrt, atan2, degrees

polyhedra = {}

# Tetrahedron
polyhedra[4] = [(+1., +1., +1.),
               (-1., -1., +1.),
               (-1., +1., -1.),
               (+1., -1., -1.)]

# Octahedron
polyhedra[6] = [(+1., +0., +0.),
                (-1., +0., +0.),
                (+0., +1., +0.),
                (+0., -1., +0.),
                (+0., +0., +1.),
                (+0., +0., -1.)]

# Cube
polyhedra[8] = [(-1., -1., -1.),
                (-1., -1., +1.),
                (-1., +1., -1.),
                (-1., +1., +1.),
                (+1., -1., -1.),
                (+1., -1., +1.),
                (+1., +1., -1.),
                (+1., +1., +1.)]

phi = (+1. + sqrt(5.)) / 2.

# Icosahedron
polyhedra[12] = [(+0., +1., +phi),
                 (+0., -1., +phi),
                 (+0., +1., -phi),
                 (+0., -1., -phi),
                 (+1., +phi, +0.),
                 (-1., +phi, +0.),
                 (+1., -phi, +0.),
                 (-1., -phi, +0.),
                 (+1., +0., +phi),
                 (-1., +0., +phi),
                 (+1., +0., -phi),
                 (-1., +0., -phi)]

# Dodecahedron
polyhedra[20] = [(+1., +1., +1.),
                (-1., +1., +1.),
                (+1., -1., +1.),
                (+1., +1., -1.),
                (-1., -1., +1.),
                (-1., +1., -1.),
                (+1., -1., -1.),
                (-1., -1., -1.),
                (+0., +1. / phi, +phi),
                (+0., -1. / phi, +phi),
                (+0., +1. / phi, -phi),
                (+0., -1. / phi, -phi),
                (+1. / phi, +phi, +0.),
                (-1. / phi, +phi, +0.),
                (+1. / phi, -phi, +0.),
                (-1. / phi, -phi, +0.),
                (+phi, +0., +1. / phi),
                (-phi, +0., +1. / phi),
                (+phi, +0., -1. / phi),
                (-phi, +0., -1. / phi)]

# Generate pre-defined viewing angles
regular_viewing_angles = {}
for p in polyhedra:
    theta, phi = [], []
    for v in polyhedra[p]:
        r = v[0] * v[0] + v[1] * v[1]
        R = v[0] * v[0] + v[1] * v[1] + v[2] * v[2]
        theta.append(degrees(atan2(r / R, v[2] / R)))
        if v[0] == 0.:
            if v[1] == 0.:
                phi.append(0.)
            else:
                phi.append(degrees(atan2(v[1] / r, 0.)))
        elif v[1] == 0.:
            phi.append(degrees(atan2(0., v[0] / r)))
        else:
            phi.append(degrees(atan2(v[1] / r, v[0] / r)))

    regular_viewing_angles[p] = (theta, phi)

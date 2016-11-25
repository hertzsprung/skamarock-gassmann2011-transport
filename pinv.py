#!/usr/bin/env python3
import numpy as np
import numpy.linalg as la
import itertools
from fractions import Fraction
from flux_divergence import split_advective_coefficients

np.set_printoptions(precision=5, linewidth=120, suppress=True)

#central = 1024
central = 1

B_fluxdiv = [[1,-2, 1,4,-2,1,-8, 4,-2], \
             [1,-1, 1,1,-1,1,-1, 1,-1], \
             [1, 0, 1,0, 0,1,0,  0,0], \
             [1, 1, 1,1, 1,1,1,  1,1], \
             [1,-2, 0,4, 0,0,-8, 0,0], \
             [1,-1, 0,1, 0,0,-1, 0,0], \
             [1, 0, 0,0, 0,0,0,  0,0], \
             [1, 1, 0,1, 0,0,1,  0,0], \
             [1,-2,-1,4, 2,1,-8,-4,-2], \
             [1,-1,-1,1, 1,1,-1,-1,-1], \
             [1, 0,-1,0, 0,1,0,  0,0], \
             [1, 1,-1,1,-1,1,1, -1,1]]

B_fluxdiv_inv = la.pinv(B_fluxdiv)
#print("B_fluxdiv_inv", B_fluxdiv_inv)

B =    [[1,-5/2, 1,25/4,-5/2,1,-125/8, 25/4,-5/2], \
        [1,-3/2, 1, 9/4,-3/2,1, -27/8,  9/4,-3/2], \
        [1,-1/2, 1, 1/4,-1/2,1,  -1/8,  1/4,-1/2], \
        [1, 1/2, 1, 1/4, 1/2,1,   1/8,  1/4, 1/2], \
        [1,-5/2, 0,25/4, 0,  0,-125/8,  0,   0], \
        [1,-3/2, 0, 9/4, 0,  0, -27/8,  0,   0], \
        central*np.array([1,-1/2, 0, 1/4, 0,  0,  -1/8,  0,   0]), \
        central*np.array([1, 1/2, 0, 1/4, 0,  0,   1/8,  0,   0]), \
        [1,-5/2,-1,25/4, 5/2,1,-125/8,-25/4,-5/2], \
        [1,-3/2,-1, 9/4, 3/2,1, -27/8, -9/4,-3/2], \
        [1,-1/2,-1, 1/4, 1/2,1,  -1/8, -1/4,-1/2], \
        [1, 1/2,-1, 1/4,-1/2,1,   1/8, -1/4, 1/2]]

Binv = la.pinv(B)
for i in range(np.shape(Binv)[0]):
    Binv[i][6] *= central
    Binv[i][7] *= central

first_deriv_x = lambda Binv, x, y: (Binv[1] + 2*Binv[3]*x + Binv[4]*y + 3*Binv[6]*x*x + 2*Binv[7]*x*y + Binv[8]*y*y).reshape((3,4))
deriv_xy = lambda Binv, x, y: (Binv[4] + 2*Binv[7] + 2*Binv[8]*y).reshape((3,4))
second_deriv_x = lambda Binv, x, y: (2*Binv[3] + 6*Binv[6]*x + 2*Binv[7]*y).reshape((3,4))
third_deriv_x = lambda Binv, x, y: (6*Binv[6]).reshape((3,4))
first_deriv_y = lambda Binv, x, y: (Binv[2] + Binv[4]*x + 2*Binv[5]*y + Binv[7]*x*x + 2*Binv[8]*x*y).reshape((3, 4))
second_deriv_y = lambda Binv, x, y: (2*Binv[5] + 2*Binv[8]*x).reshape((3, 4))
frac = np.vectorize(lambda x: Fraction(x).limit_denominator(1000))

#print("Binv\n", Binv)

def r_minus_l(coeffs):
    l = np.insert(coeffs, 4, 0, axis=1)
    r = np.insert(coeffs, 0, 0, axis=1)
    return r - l 

cubicfit_coeffs = np.reshape(Binv[0], (3,4))
cubicfit_r_minus_l = r_minus_l(cubicfit_coeffs)

print("cubicfit\n", cubicfit_coeffs)
#print("cubicfit r-l\n", cubicfit_r_minus_l)

cubic_fluxdiv = np.zeros((3,5))
cubic_fluxdiv[1] = 1/6 * np.array([0, 1, -6, 3, 2])
#print("cubic_fluxdiv\n", cubic_fluxdiv)

#quartic_fluxdiv = np.repeat([1/36 * np.array([1, -8, 0, 8, -1])], 3, axis=0)
#print("quartic_fluxdiv\n", quartic_fluxdiv)

#print("1st deriv x upupupwind\n", first_deriv_x(Binv, -5/2, 0))
#print("1st deriv x upupwind\n", first_deriv_x(Binv, -3/2, 0))
first_deriv_x_up = first_deriv_x(Binv, -1/2, 0)
print("1st deriv x upwind\n", first_deriv_x_up)
first_deriv_split = np.array([ \
        split_advective_coefficients(first_deriv_x_up[0]), \
        split_advective_coefficients(first_deriv_x_up[1]), \
        split_advective_coefficients(first_deriv_x_up[2]) \
])
first_deriv_split = np.insert(first_deriv_split, 0, 0, axis=1)
print("1st deriv x upwind split\n", first_deriv_split)

corr = first_deriv_split - cubicfit_coeffs
print("corr\n", corr)
#print("1st deriv x downwind\n", first_deriv_x(Binv, 1/2, 0))
print()
#print("2nd deriv x upupupwind\n", second_deriv_x(Binv, -5/2, 0))
#print("2nd deriv x upupwind\n", second_deriv_x(Binv, -3/2, 0))
print("2nd deriv x upwind\n", second_deriv_x(Binv, -1/2, 0))
print("2nd deriv x downwind\n", second_deriv_x(Binv, 1/2, 0))
print()
print("2nd deriv y upwind\n", second_deriv_y(Binv, -1/2, 0))
print("2nd deriv y downwind\n", second_deriv_y(Binv, 1/2, 0))
print()
#print("3rd deriv x\n", third_deriv_x(Binv, 0, 0))
#print()
#print("1st deriv y upwind\n", first_deriv_y(Binv, -1/2, 0))
#print("1st deriv y downwind\n", first_deriv_y(Binv, 1/2, 0))
#print("deriv xy\n", first_deriv_y(Binv, 0, 0))

#correction_matrix = np.transpose([second_deriv_x(Binv, -1/2, 0)[1], \
#        second_deriv_x(Binv, 1/2, 0)[1]])
correction_matrix = np.transpose([[0,1,-2,1], [1,-2,1,0]])
print("correction_matrix\n", correction_matrix)



correction_coeffs = la.lstsq(correction_matrix, corr[0])[0]
print("correction_coeffs\n", correction_coeffs)
print(np.dot(correction_matrix, correction_coeffs))

#cubicfit_correction = np.zeros((3,4))
#for i in range(np.shape(cubicfit_correction)[0]):
#    for j in range(np.shape(cubicfit_correction)[1]):
#        cubicfit_correction[i][j] = correction_coeffs[0] * second_deriv_x(Binv, -1/2, 0)[i][j] + \
#                     correction_coeffs[1] * second_deriv_x(Binv, 1/2, 0)[i][j]

#print("cubicfit_correction\n", cubicfit_correction)
#print("corrected cubicfit\n", cubicfit_coeffs + cubicfit_correction)
#print("corrected cubicfit r-l\n", r_minus_l(cubicfit_coeffs + cubicfit_correction))

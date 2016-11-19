import numpy as np
import numpy.linalg as la

class Upwind:
    def __init__(self, mesh):
        self._mesh = mesh

    def __call__(self, u, T, i):
        T_left = T[(i-1) % T.size]
        T_right = T[i % T.size]

        return -u * (T_right - T_left) / self._mesh.dx[i % T.size]

class Centred:
    def __init__(self, mesh):
        self._mesh = mesh

    def __call__(self, u, T, i):
        T_right = self.approximate(T, i)
        T_left = self.approximate(T, i+1)
        
        return -u * (T_right - T_left) / self._mesh.dx[i]

    def approximate(self, T, i):
        Cf = self._mesh.Cf[i]

        C_left = self._mesh.C[(i-1) % T.size]
        if C_left > Cf:
            C_left -= self._mesh.width

        C_right = self._mesh.C[i % T.size]
        if C_right < Cf:
            C_right += self._mesh.width

        left_weight = (C_right - Cf)/(C_right - C_left)

        return left_weight * T[(i-1) % T.size] + (1 - left_weight) * T[i % T.size]

class SkamarockGassmann:
    def __init__(self, mesh):
        self._mesh = mesh

    def __call__(self, u, T, i):
        T_right = 0.5*(T[i] + T[(i+1) % T.size]) - 1/6 * self._second_derivative(T, i)
        T_left = 0.5*(T[(i-1) % T.size] + T[i]) - 1/6 * self._second_derivative(T, i-1)

        return -u * (T_right - T_left) / self._mesh.dx[i]

    def _second_derivative(self, T, i):
        return T[(i+1)%T.size] - 2*T[i] + T[(i-1)%T.size]

class SkamarockGassmannNonUniformCentring:
    def __init__(self, mesh):
        self._centring = Centred(mesh)
        self._mesh = mesh

    def __call__(self, u, T, i):
        T_right = self._centring.approximate(T, i+1) - 1/6 * self._second_derivative(T, i)
        T_left = self._centring.approximate(T, i) - 1/6 * self._second_derivative(T, i-1)

        return -u * (T_right - T_left) / self._mesh.dx[i]

    def _second_derivative(self, T, i):
        return T[(i+1)%T.size] - 2*T[i] + T[(i-1)%T.size]

class CubicFit:
    def __init__(self, mesh):
        self._mesh = mesh
        self._Cf_coefficients = []

        for i in range(mesh.C.size):
            origin = mesh.Cf[i]

            downwind_i = i
            downwind_C = mesh.C[downwind_i]
            if downwind_C < origin:
                downwind_C += mesh.width

            upwind_i = (i-1) % mesh.C.size
            upwind_C = mesh.C[upwind_i]
            if upwind_C > origin:
                upwind_C -= mesh.width

            upupwind_i = (i-2) % mesh.C.size
            upupwind_C = mesh.C[upupwind_i]
            if upupwind_C > origin:
                upupwind_C -= mesh.width

            upupupwind_i = (i-3) % mesh.C.size
            upupupwind_C = mesh.C[upupupwind_i]
            if upupupwind_C > origin:
                upupupwind_C -= mesh.width

            stencil_C = (np.array([upupupwind_C, upupwind_C, upwind_C, downwind_C]) - origin)/mesh.mean_dx()

            B = []
            m = [1, 1, 1, 1]
            for i, C in enumerate(stencil_C):
                B.append(np.multiply(m[i], [1, C, C**2, C**3]))

            Binv = la.pinv(B)
            self._Cf_coefficients.append(Binv[0]*m[0])

    def __call__(self, u, T, i):
        T_left = self._approximate(T, i)
        T_right = self._approximate(T, i+1)
        return -u * (T_right - T_left) / self._mesh.dx[i]

    def _approximate(self, T, i):
        downwind_i = i % T.size
        upwind_i = (i-1) % T.size
        upupwind_i = (i-2) % T.size
        upupupwind_i = (i-3) % T.size

        stencil_T = [T[upupupwind_i], T[upupwind_i], T[upwind_i], T[downwind_i]]

        return np.dot(self._Cf_coefficients[i % T.size], stencil_T)

# not flux-form, so possibly not conservative?
# computes dT/dx at cell centre using an upwind-biased four-point stencil
class LeastSquaresDerivative:
    def __init__(self, mesh, polynomial_degree=3, stencil_start=-2, stencil_end=1):
        self._mesh = mesh
        self._stencil_start = stencil_start
        self._stencil_end = stencil_end
        self._advective_coeffs = []
        self._flux_coeffs_left = []
        self._flux_coeffs_right = []

        max_flux_coeff_error = 0

        for i in range(mesh.C.size):
            origin = self._mesh.C[i]

            stencil_C = []
            for stencil_i in range(stencil_start, stencil_end+1):
                C = self._mesh.C[(i+stencil_i) % mesh.C.size]
                if stencil_i > 0 and C < origin:
                    C += self._mesh.width
                elif stencil_i < 0 and C > origin:
                    C -= self._mesh.width

                stencil_C += [C]

            stencil_C = np.array(stencil_C) - origin

            B = []
            for C in stencil_C:
                terms = []
                for exponent in range(polynomial_degree+1):
                    terms += [C**exponent]
                B.append(terms)

            Binv = la.pinv(B)
            w_advective = Binv[1]

            flux_matrix = np.zeros(shape=(len(stencil_C), len(stencil_C)-1))
            for x in range(len(stencil_C)-1):
                flux_matrix[x][x] = -1
                flux_matrix[x+1][x] = 1

            w_flux = np.multiply(self._mesh.dx[i], np.dot(la.pinv(flux_matrix), w_advective))
    
            self._advective_coeffs.append(w_advective)
            self._flux_coeffs_left.append(w_flux)
            self._flux_coeffs_right.append(w_flux)

    def __call__(self, u, T, i):
        T_left = self._approximate(T, i, self._flux_coeffs_left[i])
        T_right = self._approximate(T, i+1, self._flux_coeffs_right[i])

        return -u * (T_right - T_left) / self._mesh.dx[i]

    def _approximate(self, T, i, coeffs):
        downwind_i = i % T.size
        upwind_i = (i-1) % T.size
        upupwind_i = (i-2) % T.size

        indices = []
        for index in range(i+self._stencil_start, i+self._stencil_end):
            indices += [index % T.size]

        return np.dot(coeffs, T[indices])

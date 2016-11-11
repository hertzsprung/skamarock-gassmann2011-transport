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
    def __init__(self, mesh):
        self._mesh = mesh
        self._advective_coeffs = []
        self._flux_coeffs_left = []
        self._flux_coeffs_right = []
        self._flux_coeffs_corrected = []

        max_flux_coeff_error = 0

        for i in range(mesh.C.size):
            origin = self._mesh.C[i]

            downwind_i = (i+1) % mesh.C.size
            downwind_C = self._mesh.C[downwind_i]
            if downwind_C < origin:
                downwind_C += self._mesh.width

            upwind_i = (i-1) % mesh.C.size
            upwind_C = self._mesh.C[upwind_i]
            if upwind_C > origin:
                upwind_C -= self._mesh.width

            upupwind_i = (i-2) % mesh.C.size
            upupwind_C = self._mesh.C[upupwind_i]
            if upupwind_C > origin:
                upupwind_C -= self._mesh.width

            stencil_C = np.array([upupwind_C, upwind_C, origin, downwind_C]) - origin

            B = []
            for C in stencil_C:
                B.append([1, C, C**2, C**3])

            Binv = la.pinv(B)
            w_advective = Binv[1]

            a = -w_advective[0]
            c = w_advective[3]
            b = a - w_advective[1]
            w_flux = np.multiply(self._mesh.dx[i], [a, b, c])

            max_flux_coeff_error = max(abs(b - c - w_advective[2]), max_flux_coeff_error)
    
            self._advective_coeffs.append(w_advective)
            self._flux_coeffs_left.append(w_flux)
            self._flux_coeffs_right.append(w_flux)

        assert max_flux_coeff_error < 1e-12

        for l, r in zip(self._rotate(self._flux_coeffs_left, 1), self._flux_coeffs_right):
            # TODO: can we match adjacent r/l flux coeff pairs and adjust them so that they are equal?
            self._flux_coeffs_corrected.append(0.5 * (l + r))
#            print(r, l, 0.5 * (l + r), la.norm(r - l))

    def _rotate(self, l, n):
        return l[n:] + l[:n]

    def __call__(self, u, T, i):
        T_left = self._approximate(T, i, self._flux_coeffs_left[i])
        T_right = self._approximate(T, i+1, self._flux_coeffs_right[i])

        return -u * (T_right - T_left) / self._mesh.dx[i]

    def _approximate(self, T, i, coeffs):
        downwind_i = i % T.size
        upwind_i = (i-1) % T.size
        upupwind_i = (i-2) % T.size

        stencil_T = [T[upupwind_i], T[upwind_i], T[downwind_i]]

        return np.dot(coeffs, stencil_T)

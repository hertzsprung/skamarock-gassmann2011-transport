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
        # FIXME: this should account for nonuniform mesh geometries
        T_right = 0.5*(T[i] + T[(i+1) % T.size])
        T_left = 0.5*(T[(i-1) % T.size] + T[i])
        
        return -u * (T_right - T_left) / self._mesh.dx[i]

class SkamarockGassmann:
    def __init__(self, mesh):
        self._mesh = mesh

    def __call__(self, u, T, i):
        T_right = 0.5*(T[i] + T[(i+1) % T.size]) - 1/6 * self._second_derivative(T, i)
        T_left = 0.5*(T[(i-1) % T.size] + T[i]) - 1/6 * self._second_derivative(T, i-1)

        return -u * (T_right - T_left) / self._mesh.dx[i]

    def _second_derivative(self, T, i):
        return T[(i+1)%T.size] - 2*T[i] + T[(i-1)%T.size]

# not flux-form, so possibly not conservative?
# computes dT/dx at cell centre using an upwind-biased four-point stencil
class LeastSquaresDerivative:
    def __init__(self, mesh):
        self._mesh = mesh
        self._advective_coeffs = []
#        self._flux_coeffs = []

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
#            a = -w_advective[0]
#            c = w_advective[3]
#            b = a - w_advective[1]
#            w_flux = [a, b, c]
#            assert abs(b - c - w_advective[2]) < 1e-10

            self._advective_coeffs.append(w_advective)

    def __call__(self, u, T, i):
        origin = self._mesh.C[i]

        downwind_i = (i+1) % T.size
        upwind_i = (i-1) % T.size
        upupwind_i = (i-2) % T.size

        stencil_T = [T[upupwind_i], T[upwind_i], T[i], T[downwind_i]]

        first_derivative = np.dot(self._advective_coeffs[i], stencil_T)
        return -u * first_derivative

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
        T_right = self._approximate(T, i+1)
        T_left = self._approximate(T, i)
        return -u * (T_right - T_left) / self._mesh.dx[i]

    def _approximate(self, T, i):
        downwind_i = i % T.size
        upwind_i = (i-1) % T.size
        upupwind_i = (i-2) % T.size
        upupupwind_i = (i-3) % T.size

        stencil_T = [T[upupupwind_i], T[upupwind_i], T[upwind_i], T[downwind_i]]

        return np.dot(self._Cf_coefficients[i % T.size], stencil_T)

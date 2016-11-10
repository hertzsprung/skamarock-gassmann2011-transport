import numpy as np
import numpy.linalg as la

class Centred:
    def __call__(self, mesh, u, T, i):
        T_right = 0.5*(T[i] + T[(i+1) % T.size])
        T_left = 0.5*(T[(i-1) % T.size] + T[i])

        return -u * (T_right - T_left) / mesh.dx[i]

class SkamarockGassmann:
    def __call__(self, mesh, u, T, i):
        T_right = 0.5*(T[i] + T[(i+1) % T.size]) - 1/6 * self._second_derivative(T, i)
        T_left = 0.5*(T[(i-1) % T.size] + T[i]) - 1/6 * self._second_derivative(T, i-1)

        return -u * (T_right - T_left) / mesh.dx[i]

    def _second_derivative(self, T, i):
        return T[(i+1)%T.size] - 2*T[i] + T[(i-1)%T.size]

# not flux-form, so possibly not conservative?
# computes dT/dx at cell centre using an upwind-biased four-point stencil
class LeastSquaresDerivative:
    def __call__(self, mesh, u, T, i):
        origin = mesh.C[i]

        downwind_i = (i+1) % T.size
        downwind_C = mesh.C[downwind_i]
        if downwind_C < origin:
            downwind_C += mesh.width

        upwind_i = (i-1) % T.size
        upwind_C = mesh.C[upwind_i]
        if upwind_C > origin:
            upwind_C -= mesh.width

        upupwind_i = (i-2) % T.size
        upupwind_C = mesh.C[upupwind_i]
        if upupwind_C > origin:
            upupwind_C -= mesh.width

        stencil_C = np.array([upupwind_C, upwind_C, origin, downwind_C]) - origin
        stencil_T = [T[upupwind_i], T[upwind_i], T[i], T[downwind_i]]

        B = []
        for C in stencil_C:
            B.append([1, C, C**2, C**3])

        Binv = la.pinv(B)
        first_derivative = np.dot(Binv[1], stencil_T)
        return -u * first_derivative

class CubicFit:
    def __call__(self, mesh, u, T, i):
        T_right = self._approximate(mesh, T, i)
        T_left = self._approximate(mesh, T, i-1)
        return -u * (T_right - T_left) / mesh.dx[i]

    def _approximate(self, mesh, T, i):
        origin = mesh.Cf[i+1]

        downwind_i = (i+1) % T.size
        downwind_C = mesh.C[downwind_i]
        if downwind_C < origin:
            downwind_C += mesh.width

        upwind_i = i
        upwind_C = mesh.C[upwind_i]
        if upwind_C > origin:
            upwind_C -= mesh.width

        upupwind_i = (i-1) % T.size
        upupwind_C = mesh.C[upupwind_i]
        if upupwind_C > origin:
            upupwind_C -= mesh.width

        upupupwind_i = (i-2) % T.size
        upupupwind_C = mesh.C[upupupwind_i]
        if upupupwind_C > origin:
            upupupwind_C -= mesh.width

        stencil_C = (np.array([upupupwind_C, upupwind_C, upwind_C, downwind_C]) - origin)/mesh.mean_dx()
        stencil_T = [T[upupupwind_i], T[upupwind_i], T[upwind_i], T[downwind_i]]

        B = []
        m = [1, 1, 1, 1]
        for i, C in enumerate(stencil_C):
            B.append(np.multiply(m[i], [1, C, C**2, C**3]))

        Binv = la.pinv(B)
        return np.dot(Binv[0]*m[0], stencil_T)

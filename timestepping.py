import numpy as np

class ForwardEuler:
    def __init__(self, flux_divergence):
        self._flux_divergence = flux_divergence

    def __call__(self, T_old, dt, u):
        T = np.zeros_like(T_old)

        for i in range(T.size):
            T[i] = T_old[i] + dt*self._flux_divergence(u, T_old, i)

        return T

class RungeKutta3Stage2ndOrder:
    def __init__(self, flux_divergence):
        self._flux_divergence = flux_divergence

    def __call__(self, T_old, dt, u):
        T_star = self._stage(dt, T_old, T_old, u)
        T_star_star = self._stage(dt, T_old, T_star, u)
        T = self._stage(dt, T_old, T_star_star, u)
        return T

    def _stage(self, dt, T_old, T_fractional, u):
        T = np.zeros_like(T_old)

        for i in range(T.size):
            T[i] = T_old[i] + 0.5*dt*(self._flux_divergence(u, T_old, i) + self._flux_divergence(u, T_fractional, i))

        return T

class RungeKutta3:
    def __init__(self, flux_divergence):
        self._flux_divergence = flux_divergence

    def __call__(self, T_old, dt, u):
        T_star = self._stage(dt/3, T_old, T_old, u)
        T_star_star = self._stage(dt/2, T_old, T_star, u)
        T = self._stage(dt, T_old, T_star_star, u)
        return T

    def _stage(self, dt_fractional, T_old, T_fractional, u):
        T = np.zeros_like(T_old)

        for i in range(T.size):
            T[i] = T_old[i] + dt_fractional*self._flux_divergence(u, T_fractional, i)

        return T

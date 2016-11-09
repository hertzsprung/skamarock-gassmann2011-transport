import math
import numpy as np
import os

class SimulationSpec:
    def __init__(self, initial, dt = 0.05, nx = 10):
        self._initial = initial
        self._dt = dt
        self._nx = nx

        self._width = 1.0
        self._end_time = 1.0
        self._u = 1
        self.times = np.arange(0, self._end_time+np.finfo(float).eps, step=self._dt)
        self.C = np.linspace(0.5*self._width/self._nx, 1-0.5*self._width/self._nx, num=self._nx)
        self._Cf = np.linspace(0, 1, num=self._nx+1)
        self.dx = [r - l for l,r in zip(self._Cf[:-1], self._Cf[1:])]
        self.T = np.array(np.array([initial(C) for C in self.C]))

    def advect(self):
        simulation = Simulation(self)

        for t_old, t in zip(self.times[:-1], self.times[1:]):
            simulation.append(self._advect(simulation, t - t_old))

        return simulation

    def _advect(self, simulation, dt):
        T_old = simulation.T_latest()

        T_star = self._runge_kutta_stage(dt, T_old, T_old)
        T_star_star = self._runge_kutta_stage(dt, T_old, T_star)
        T = self._runge_kutta_stage(dt, T_old, T_star_star)

        return T

    def _runge_kutta_stage(self, dt, T_old, T_fractional):
        T = np.zeros_like(T_old)

        for i in range(T.size):
            T[i] = T_old[i] + 0.5*dt*(self._flux_divergence(T_old, i) + self._flux_divergence(T_fractional, i))

        return T

    def _flux_divergence(self, T, i):
        T_right = 0.5*(T[i] + T[(i+1) % T.size])
        T_left = 0.5*(T[(i-1) % T.size] + T[i])

        return -self._u * (T_right - T_left) / self.dx[i]

    def mean_dx(self):
        return np.mean(self.dx)

    def refine(self):
        return SimulationSpec(self._initial, self._dt/2, self._nx*2)

    def max_courant(self):
        return np.max(self._u * self._dt / self._dx)

class Simulation:
    def __init__(self, spec):
        self._spec = spec
        self.T = [spec.T]

    def append(self, T):
        self.T = np.append(self.T, [T], axis=0)

    def T_initial(self):
        return self.T[0]

    def T_latest(self):
        return self.T[-1]

    def error_l2(self):
        return math.sqrt(np.sum((self.T_latest() - self.T_initial())**2 * self._spec.dx)) / np.sum(self.T_initial()**2 * self._spec.dx)

    def dump(self, directory):
        for i, T in enumerate(self.T):
            with open(os.path.join(directory, str(self._spec.times[i])) + ".dat", 'w') as f:
                for C, T_value in zip(self._spec.C, T):
                    print(C, T_value, file=f)

class Convergence:
    def __init__(self, spec):
        self._spec = spec
        self._dxs = []
        self._errors_l2 = []

    def converge(self):
        simulation = self._spec.advect()
        self._dxs.append(self._spec.mean_dx())
        self._errors_l2.append(simulation.error_l2())
        self._spec = self._spec.refine()

    def order(self):
        log_dxs = np.log(self._dxs)
        log_errors_l2 = np.log(self._errors_l2)

        A = np.vstack([log_dxs, np.ones(len(log_dxs))]).T
        m, c = np.linalg.lstsq(A, log_errors_l2)[0]
        return m

class InitialConditions:
    @staticmethod
    def sine_wave():
        return lambda x: math.sin(2*math.pi*x)

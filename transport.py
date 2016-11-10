import math
import numpy as np
import os
import itertools
import random

class InitialConditions:
    @staticmethod
    def sine_wave():
        return lambda x: math.sin(2*math.pi*x)

class Mesh:
    @staticmethod
    def uniform(nx = 10):
        return lambda width: np.linspace(0.5*width/nx, 1-0.5*width/nx, num=nx)

    @staticmethod
    def nonuniform(nx = 10, nonuniformity=0.5):
        return lambda width: Mesh._nonuniform(width, nx, nonuniformity)

    @staticmethod
    def _nonuniform(width, nx, nonuniformity, seed=0):
        random.seed(seed)
        return [C + random.uniform(-1,1)*nonuniformity*width/nx for C in Mesh.uniform(nx)(width)]

    def __init__(self, C = uniform.__func__(), width = 1.0):
        self.width = width
        self.C = C(width)
        self.Cf = [0]
        self.Cf += [0.5*(l + r) for l, r in zip(self.C[:-1], self.C[1:])]
        self.Cf += [width]
        self.dx = [r - l for l,r in zip(self.Cf[:-1], self.Cf[1:])]

    def refine(self):
        refined_Cf = list(self.roundrobin(self.Cf, self.C))
        refined_C = [0.5*(f0+f1) for f0, f1 in zip(refined_Cf[:-1], refined_Cf[1:])]

        return Mesh(C = lambda width: refined_C, width = self.width)

    def roundrobin(self, *iterables):
        "roundrobin('ABC', 'D', 'EF') --> A D E B F C"
        # Recipe credited to George Sakkis
        pending = len(iterables)
        nexts = itertools.cycle(iter(it).__next__ for it in iterables)
        while pending:
            try:
                for next in nexts:
                    yield next()
            except StopIteration:
                pending -= 1
                nexts = itertools.cycle(itertools.islice(nexts, pending))

class SimulationSpec:
    def __init__(self, initial = InitialConditions.sine_wave(), dt = 0.05, mesh = Mesh()):
        self._initial = initial
        self.mesh = mesh
        self._dt = dt

        self._width = 1.0
        self._end_time = 1.0
        self._u = 1
        self.times = np.arange(0, self._end_time+np.finfo(float).eps, step=self._dt)
        self.T = np.array(np.array([initial(C) for C in mesh.C]))

    def advect(self):
        simulation = Simulation(self)

        for t_old, t in zip(self.times[:-1], self.times[1:]):
            simulation.append(self._advect(simulation, t - t_old))

        return simulation

    def _advect(self, simulation, dt):
        T_old = simulation.T_latest()

        T_star = self._runge_kutta_stage(dt/3, T_old, T_old)
        T_star_star = self._runge_kutta_stage(dt/2, T_old, T_star)
        T = self._runge_kutta_stage(dt, T_old, T_star_star)

        return T

    def _runge_kutta_stage(self, dt_fractional, T_old, T_fractional):
        T = np.zeros_like(T_old)

        for i in range(T.size):
            T[i] = T_old[i] + dt_fractional*self._flux_divergence(T_fractional, i)

        return T

    def _flux_divergence(self, T, i):
        T_right = 0.5*(T[i] + T[(i+1) % T.size]) - 1/6 * self._second_derivative(T, i)
        T_left = 0.5*(T[(i-1) % T.size] + T[i]) - 1/6 * self._second_derivative(T, i-1)

        return -self._u * (T_right - T_left) / self.mesh.dx[i]

    def _second_derivative(self, T, i):
        return T[(i+1)%T.size] - 2*T[i] + T[(i-1)%T.size]

    def mean_dx(self):
        return np.mean(self.mesh.dx)

    def refine(self):
        return SimulationSpec(self._initial, self._dt/2, self.mesh.refine())

    def max_courant(self):
        return np.max(self._u * self._dt / self.mesh.dx)

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
        return math.sqrt(np.sum((self.T_latest() - self.T_initial())**2 * self._spec.mesh.dx)) / np.sum(self.T_initial()**2 * self._spec.mesh.dx)

    def error_linf(self):
        return np.max(np.absolute(self.T_latest() - self.T_initial())) / np.max(np.absolute(self.T_initial()))

    def dump(self, directory):
        for i, T in enumerate(self.T):
            with open(os.path.join(directory, str(self._spec.times[i])) + ".dat", 'w') as f:
                for C, T_value in zip(self._spec.mesh.C, T):
                    print(C, T_value, file=f)

class Convergence:
    def __init__(self, spec):
        self._spec = spec
        self._dxs = []
        self.errors_l2 = []
        self.errors_linf = []

    def converge(self):
        simulation = self._spec.advect()
        self._dxs.append(self._spec.mean_dx())
        self.errors_l2.append(simulation.error_l2())
        self.errors_linf.append(simulation.error_linf())
        self._spec = self._spec.refine()
        return simulation

    def order(self, errors = None):
        if errors is None:
            errors = self.errors_l2

        log_dxs = np.log(self._dxs)
        log_errors = np.log(errors)

        A = np.vstack([log_dxs, np.ones(len(log_dxs))]).T
        m, c = np.linalg.lstsq(A, log_errors)[0]
        return m

import math
import numpy as np
import os
import itertools
import random
import flux_divergence
import timestepping

class InitialConditions:
    @staticmethod
    def sine_wave(k=1):
        return lambda x: math.sin(2*math.pi*k*x)

    def quadratic():
        return lambda x: (x-0.5)**2

class Mesh:
    @staticmethod
    def uniform(nx = 10):
        return lambda width: np.linspace(0.0, 1.0, num=nx+1)

    @staticmethod
    def nonuniform(nx = 10, nonuniformity=0.5):
        return lambda width: Mesh._nonuniform(width, nx, nonuniformity)

    @staticmethod
    def _nonuniform(width, nx, nonuniformity, seed=0):
        random.seed(seed)
        internal_faces = [Cf + random.uniform(-1,1)*nonuniformity*width/nx for Cf in Mesh.uniform(nx)(width)[1:-1]]
        return list(itertools.chain([0.0], internal_faces, [width]))

    def __init__(self, Cf = uniform.__func__(), width = 1.0):
        self.width = width
        self.Cf = np.array(Cf(width))
        self.C = np.array([0.5*(l + r) for l,r in zip(self.Cf[:-1], self.Cf[1:])])
        self.dx = [r - l for l,r in zip(self.Cf[:-1], self.Cf[1:])]

    def mean_dx(self):
        return np.mean(self.dx)

    def refine(self):
        refined_Cf = list(self.roundrobin(self.Cf, self.C))

        return Mesh(Cf = lambda width: refined_Cf, width = self.width)

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
    def __init__(self, initial = InitialConditions.sine_wave(), dt = 0.025, end_time = 1.0, mesh = Mesh(), flux_divergence=flux_divergence.SkamarockGassmann, timestepping = timestepping.RungeKutta3):
        self._initial = initial
        self.mesh = mesh
        self._dt = dt
        self._flux_divergence = flux_divergence(mesh)
        self._timestepping = timestepping(self._flux_divergence)

        self._width = 1.0
        self._end_time = end_time
        self._u = 1.0
        self.times = np.arange(0.0, self._end_time+np.finfo(float).eps, step=self._dt)
        self.T = np.array(np.array([initial(C) for C in mesh.C]))

    def advect(self):
        simulation = Simulation(self)

        for t_old, t in zip(self.times[:-1], self.times[1:]):
            simulation.append(self._timestepping(simulation.T_latest(), t - t_old, self._u))

        return simulation

    def refine(self):
        return SimulationSpec(self._initial, self._dt/2, self._end_time, self.mesh.refine(), self._flux_divergence.__class__, self._timestepping.__class__)

    def max_courant(self):
        return self._u * self._dt / np.min(self.mesh.dx)

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

    def mass_change(self):
        return np.sum(self._spec.mesh.dx * self.T_latest()) - np.sum(self._spec.mesh.dx*self.T_initial())

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
        self._dxs.append(self._spec.mesh.mean_dx())
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

    def dump(self, filename):
        with open(filename, 'w') as f: 
            for mean_dx, l2, linf in zip(self._dxs, self.errors_l2, self.errors_linf):
                print(mean_dx, l2, linf, file=f)

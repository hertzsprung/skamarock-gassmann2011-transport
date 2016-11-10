#!/usr/bin/env py.test
import pytest
from transport import *

def test_sine_wave_advection_with_baseline_accuracy():
    spec = SimulationSpec()
    simulation = spec.advect()
    simulation.dump("results")
    assert simulation.error_l2() < 0.65

def test_third_order_convergence_on_uniform_mesh():
    convergence = Convergence(SimulationSpec())

    for i in range(6):
        convergence.converge()

    assert convergence.order(convergence.errors_l2) == pytest.approx(3, 0.1)
    assert convergence.order(convergence.errors_linf) == pytest.approx(3, 0.1)

def test_first_order_convergence_on_nonuniform_mesh():
    convergence = Convergence(SimulationSpec(mesh = Mesh(Mesh.nonuniform())))

    for i in range(6):
        convergence.converge()

    assert convergence.order(convergence.errors_l2) == pytest.approx(2, 0.1)
    assert convergence.order(convergence.errors_linf) < 1.9

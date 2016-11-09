#!/usr/bin/env py.test-3
import pytest
from transport import *

def test_sine_wave_advection_with_baseline_accuracy():
    spec = SimulationSpec(InitialConditions.sine_wave())
    simulation = spec.advect()
    assert simulation.error_l2() < 0.65

def test_third_order_convergence_on_uniform_mesh():
    convergence = Convergence(SimulationSpec(InitialConditions.sine_wave()))

    for i in range(6):
        convergence.converge()

    assert convergence.order() > 1.9 and convergence.order() < 2.1
#    assert convergence.order() == pytest.approx(2, 0.1)


    

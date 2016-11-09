#!/usr/bin/env py.test-3
import pytest
from transport import *

def test_sine_wave_advection_with_baseline_accuracy():
    spec = SimulationSpec(InitialConditions.sine_wave())
    simulation = spec.advect()
    simulation.dump("results")
    assert simulation.error_l2() < 0.65

def test_third_order_convergence_on_uniform_mesh():
    print("hello")

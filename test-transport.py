import pytest
from transport import *
from flux_divergence import *

def test_sine_wave_advection_with_skamarock_gassmann():
    spec = SimulationSpec(mesh = Mesh(Mesh.nonuniform()))
    simulation = spec.advect()
    simulation.dump("results")
    assert simulation.error_l2() < 0.25

def test_sine_wave_advection_with_least_squares():
    spec = SimulationSpec(mesh = Mesh(Mesh.nonuniform()), flux_divergence=LeastSquaresDerivative())
    simulation = spec.advect()
    simulation.dump("results")
    assert simulation.error_l2() < 0.2

def test_sine_wave_advection_with_cubic_fit():
    spec = SimulationSpec(mesh = Mesh(Mesh.nonuniform()), flux_divergence=CubicFit())
    simulation = spec.advect()
    simulation.dump("results")
    assert simulation.error_l2() < 0.2

def test_third_order_convergence_on_uniform_mesh_with_skamarock_gassmann():
    convergence = Convergence(SimulationSpec())

    for i in range(6):
        convergence.converge()

    convergence.dump("results/sk-uniform.dat")
    assert convergence.order(convergence.errors_l2) == pytest.approx(3, 0.1)
    assert convergence.order(convergence.errors_linf) == pytest.approx(3, 0.1)

def test_first_order_convergence_on_nonuniform_mesh_with_skamarock_gassmann():
    convergence = Convergence(SimulationSpec(mesh = Mesh(Mesh.nonuniform())))

    for i in range(6):
        convergence.converge()

    convergence.dump("results/sk-nonuniform.dat")
    assert convergence.order(convergence.errors_l2) == pytest.approx(2, 0.1)
    assert convergence.order(convergence.errors_linf) < 1.9

def test_third_order_convergence_on_uniform_mesh_with_least_squares():
    convergence = Convergence(SimulationSpec(flux_divergence=LeastSquaresDerivative()))

    for i in range(6):
        convergence.converge()

    convergence.dump("results/ls-uniform.dat")
    assert convergence.order(convergence.errors_l2) == pytest.approx(3, 0.1)
    assert convergence.order(convergence.errors_linf) == pytest.approx(3, 0.1)

def test_first_order_convergence_on_nonuniform_mesh_with_least_squares():
    convergence = Convergence(SimulationSpec(mesh = Mesh(Mesh.nonuniform()), flux_divergence=LeastSquaresDerivative()))

    for i in range(6):
        convergence.converge()

    convergence.dump("results/ls-nonuniform.dat")
    assert convergence.order(convergence.errors_l2) == pytest.approx(2, 0.2)
    assert convergence.order(convergence.errors_linf) < 1.9

def test_third_order_convergence_on_uniform_mesh_with_cubic_fit():
    convergence = Convergence(SimulationSpec(flux_divergence=CubicFit()))

    for i in range(6):
        convergence.converge()

    convergence.dump("results/cf-uniform.dat")
    assert convergence.order(convergence.errors_l2) == pytest.approx(3, 0.1)
    assert convergence.order(convergence.errors_linf) == pytest.approx(3, 0.1)

def test_first_order_convergence_on_nonuniform_mesh_with_cubic_fit():
    convergence = Convergence(SimulationSpec(mesh = Mesh(Mesh.nonuniform()), flux_divergence=CubicFit()))

    for i in range(6):
        convergence.converge()

    convergence.dump("results/cf-nonuniform.dat")
    assert convergence.order(convergence.errors_l2) == pytest.approx(2, 0.2)
    assert convergence.order(convergence.errors_linf) < 1.9

import pytest
import itertools
from transport import *
from flux_divergence import *
from timestepping import *

def check_l2_error_below(threshold, mesh_type, flux_divergence):
    spec = SimulationSpec(mesh = Mesh(mesh_type), flux_divergence = flux_divergence)
    simulation = spec.advect()
    assert simulation.error_l2() < threshold

def check_convergence(mesh_type, flux_divergence, l2_error, linf_error, dump_file = None):
    convergence = Convergence(SimulationSpec(mesh = Mesh(mesh_type), flux_divergence = flux_divergence))

    for i in range(6):
        convergence.converge()

    if dump_file:
        convergence.dump(dump_file)

    assert convergence.order(convergence.errors_l2) == pytest.approx(l2_error, 0.2)
    assert convergence.order(convergence.errors_linf) == pytest.approx(linf_error, 0.2)

def check_conservation(flux_divergence):
    spec = SimulationSpec(mesh = Mesh(Mesh.nonuniform()), flux_divergence = flux_divergence)
    simulation = spec.advect()
    assert simulation.mass_change() == pytest.approx(0)

def test_sine_wave_advection_with_skamarock_gassmann():
    check_l2_error_below(0.35, Mesh.nonuniform(), SkamarockGassmann)

def test_sine_wave_advection_with_least_squares():
    check_l2_error_below(0.2, Mesh.nonuniform(), LeastSquaresDerivative)

def test_sine_wave_advection_with_skamarock_gassmann_nonuniform_centring():
    check_l2_error_below(0.25, Mesh.nonuniform(), SkamarockGassmannNonUniformCentring)

def test_sine_wave_advection_with_cubic_fit():
    check_l2_error_below(0.2, Mesh.nonuniform(), CubicFit)

def test_third_order_convergence_on_uniform_mesh_with_skamarock_gassmann():
    check_convergence(Mesh.uniform(), SkamarockGassmann, l2_error=3, linf_error=3, dump_file="results/sk-uniform.dat")

def test_less_than_second_order_convergence_on_nonuniform_mesh_with_skamarock_gassmann():
    check_convergence(Mesh.nonuniform(), SkamarockGassmann, l2_error=2, linf_error=1.5, dump_file="results/sk-nonuniform.dat")

def test_third_order_convergence_on_uniform_mesh_with_skamarock_gassmann_nonuniform_centring():
    check_convergence(Mesh.uniform(), SkamarockGassmannNonUniformCentring, l2_error=3, linf_error=3, dump_file="results/sknuc-uniform.dat")

def test_less_than_second_order_convergence_on_nonuniform_mesh_with_skamarock_gassmann_nonuniform_centring():
    check_convergence(Mesh.nonuniform(), SkamarockGassmannNonUniformCentring, l2_error=2, linf_error=1.5, dump_file="results/sknuc-nonuniform.dat")

def test_third_order_convergence_on_uniform_mesh_with_least_squares():
    check_convergence(Mesh.uniform(), LeastSquaresDerivative, l2_error=3, linf_error=3, dump_file="results/ls-uniform.dat")

def test_third_order_convergence_on_nonuniform_mesh_with_least_squares():
    check_convergence(Mesh.nonuniform(), LeastSquaresDerivative, l2_error=3, linf_error=3, dump_file="results/ls-nonuniform.dat")

def test_second_order_convergence_on_uniform_mesh_with_cubic_fit():
    check_convergence(Mesh.uniform(), CubicFit, l2_error=2, linf_error=2, dump_file="results/cf-uniform.dat")

def test_second_order_convergence_on_nonuniform_mesh_with_cubic_fit():
    check_convergence(Mesh.nonuniform(), CubicFit, l2_error=2, linf_error=2, dump_file="results/cf-nonuniform.dat")

def test_second_order_convergence_on_uniform_mesh_with_centred():
    check_convergence(Mesh.uniform(), Centred, l2_error=2, linf_error=2, dump_file="results/centred-uniform.dat")

def test_second_order_convergence_on_nonuniform_mesh_with_centred():
    check_convergence(Mesh.nonuniform(), Centred, l2_error=2, linf_error=2, dump_file="results/centred-nonuniform.dat")

def test_skamarock_gassmann_is_conservative():
    check_conservation(SkamarockGassmann)

def test_cubic_fit_is_conservative():
    check_conservation(CubicFit)

def test_least_squares_is_nonconservative():
    spec = SimulationSpec(mesh = Mesh(Mesh.nonuniform()), flux_divergence = LeastSquaresDerivative)
    simulation = spec.advect()
    assert simulation.mass_change() > 1e-12

def test_centred_is_conservative():
    check_conservation(Centred)

def test_nonuniform_mesh_geometry():
    mesh = Mesh(Cf = lambda width: [0, 0.8, 1])

    assert mesh.C[0] == pytest.approx(0.4)
    assert mesh.C[1] == pytest.approx(0.9)

    assert mesh.dx[0] == pytest.approx(0.8)
    assert mesh.dx[1] == pytest.approx(0.2)

def test_cubic_fit_coefficients_mid_way_between_third_and_fourth_order():
    mesh = Mesh(Mesh.uniform(nx=5))
    cubic_fit = CubicFit(mesh)
    cubic_fit_flux_l = np.array(list(itertools.chain(cubic_fit._Cf_coefficients[3], [0])))
    cubic_fit_flux_r = np.array(list(itertools.chain([0], cubic_fit._Cf_coefficients[4])))
    cubic_fit_flux_div = cubic_fit_flux_r - cubic_fit_flux_l
    print(cubic_fit_flux_div)

    third_order = LeastSquaresDerivative(mesh)
    third_order_flux_l = np.array(list(itertools.chain(third_order._flux_coeffs_left[2], [0])))
    third_order_flux_r = np.array(list(itertools.chain([0], third_order._flux_coeffs_right[2])))
    third_order_flux_div = np.array(list(itertools.chain([0], third_order_flux_r - third_order_flux_l)))
    print(third_order_flux_div)

    fourth_order = LeastSquaresDerivative(mesh, polynomial_degree=4, stencil_start=-2, stencil_end=2)
    fourth_order_flux_l = np.array(list(itertools.chain(fourth_order._flux_coeffs_left[2], [0])))
    fourth_order_flux_r = np.array(list(itertools.chain([0], fourth_order._flux_coeffs_right[2])))
    fourth_order_flux_div = fourth_order_flux_r - fourth_order_flux_l
    print(fourth_order_flux_div)

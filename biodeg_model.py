"""
A model for a a typical biodegradation experiment where droplets are being
mixed and kept in suspension.
"""
import numpy as np
from IPython.terminal.debugger import set_trace as debug
from pint import UnitRegistry
unit = UnitRegistry()


def nr_droplets(diameter, a=2.5, k=13):
    """
    Diameter in micrometer
    """

    return k * diameter ** (-a)


def droplet_area(D):
    """
    A = pi r ** 2 = 0.25 * pi * D ** 2
    """
    return 0.25 * np.pi * D ** 2


def droplet_volume(D):
    """
    V = 4/3 pi ** r ** 3 = 4/3 pi * (D/2) ** 3 = 4/3 * 1/8 D ** 3
    V = 1/6 * pi * D ** 3
    """

    return np.pi * D ** 3 / 6


def diam_from_vol(V):
    """
    V = np.pi * D ** 3 / 6

    D = (V * 6 / np.pi) ** (1/3)
    """

    return (V * 6 / np.pi) ** (1 / 3)


def main():
    """
    Question: do I need to use a solver? I think that the system has an analytic solution.

    m_droplet(t) = np.exp((-k1 - k2) * t )
    m_water(t) = np.exp(k2 * t)
    m_biotransf_and_biological = np.exp(k1 * t)
    D_droplet = ... something.
    """

    nr_pseudocomponents = 2
    nr_sizebins = 3

    # Droplet sizes (um) and number
    droplet_sizes = np.logspace(1, 7, num=nr_sizebins, base=3.3) * unit('micrometer')
    droplet_number = nr_droplets(droplet_sizes.magnitude, k=2e9)
    droplet_volumes = droplet_volume(droplet_sizes)
    droplet_area = droplet_volume(droplet_sizes)

    # Rates
    droplet_rates = np.random.random(nr_pseudocomponents) / unit('second')
    waf_rates = droplet_rates * 3 / unit('second')
    dissolution_propensities = np.random.random(nr_pseudocomponents) / unit('second')

    ## Initial states ##

    # Mass fraction of each component in each droplet
    composition = np.random.random(nr_pseudocomponents)
    composition = composition / sum(composition)

    # Component density
    rho_min = 700
    rho_max = 1010
    component_densities = rho_max - np.random.random(nr_pseudocomponents) * (rho_max - rho_min)
    component_densities = component_densities * unit('kg') / unit('m') ** 3
    oil_density = np.sum(component_densities * composition)

    # Components in water
    water_components = np.zeros(nr_pseudocomponents)

    # secondary metabolites in water
    water_seconday_metabolites = np.zeros(nr_pseudocomponents)

    # Oil concentration
    oil_conc_ppm = 5 * 1e-6

    # Water volume
    tank_volume = 0.2 * unit('litre')
    water_rho = 1027 * unit('kg') / unit('m') ** 3
    tank_mass = tank_volume * water_rho
    tank_mass = tank_mass.to_base_units()

    # Oil mass (from ppm)
    oil_mass = oil_conc_ppm * tank_mass

    # Initial distribution of component mass per droplet
    oil_component_mass = oil_mass * composition

    # Need to pre-build array for mass of pseudocomponent per droplet size bin.
    # That is equal to mass per size bin * pseudocompoent fraction
    mass_per_size_class = droplet_volumes * droplet_number * oil_density
    debug()

    setup_yvector(droplet_sizes, mass_per_size_class, component_fractions)


def dissolution_rate(droplet_masses, droplet_densities):
    """
    Input: droplet masses per size bin and densities.

    Calculate the area of each droplet, find rate by Area * dissolution_propensity.
    """


def setup_yvector(bin_diameters, mass_per_size_class, component_fractions):
    """
    y = [bin1_comp1, bin1_comp2, ... bin1_compN, bin2_comp1, ... binM_compN, wat_comp1, ...
         wat_compN, wat_second_metabolit1, ...,  wat_second_metabolitN, bin1_diam, ..., binM_diam]
    """

    y_vector = []

    # First the mass of each component in each droplet size bin
    for bin_mass in mass_per_size_class:
        for comp_fraction in nr_components:
            y_vector.append(bin_mass * comp_fraction)

    # Then the mass of each component in water (initially 0)
    for i_comp in range(len(component_fractions)):
        y_vector.append(0)

    # Then the diameter of droplets in each bin
    for bin_diameter in bin_diameters:
        y_vector.append(bin_diameter)


def generate_rhs(nr_sizebins, nr_components):
    """
    Generate code for scipy.ode differential equation

    y = [bin1_comp1, bin1_comp2, ... bin1_compN, bin2_comp1, ... binM_compN, wat_comp1, ...
         wat_compN, wat_second_metabolit1, ...,  wat_second_metabolitN, bin1_diam, ..., binM_diam]

    The rate for the components is, where di,j means size bin i, component j.
    dmdrop_i,j/dt = (-kbio_i,j - kdiss_i,j) * y[i + j]

    The rate for transfer of components from dissolution to water components is
    dmwat_j/dt = sum(kdiss_:,j * y[: + j])

    The rate for the bin diameters is
    dD_i/dt = 1/3 * sum(y[i: i + nr_components]) ** (-2/3) * K,
    where K = ( (rho_i) * 6 / np.pi) ** 1/3

           0 ,   1 ,  2  ,  3  ,  4  ,  5  ,  6  ,  7  , 8
    y = [m1,1, m1,2, m1,3, m2,1, m2,2, m2,3, m3,1, m3,2, ...]
    y_idx [i, j] = nr_components * (i-1) + j

    y_idx [3,2 ] = 3 * (2) + 2 = 7, correct.
    """
    code = ''

    # Start
    code += '['

    # Mass transfer from droplets to water and microorganisms
    for bin_nr in range(nr_sizebins):
        for comp_nr in range(nr_components):
            y_idx = nr_components * bin_nr + nr_components 
            code += '(-kbio_drop[{0}, {1}] - kdiss[{0}, {1}) * y[{2}]'.format(bin_nr, comp_nr, y_idx)

    # Mass transfer to the water in the experiment
    for comp_nr in range(nr_components):
        ode_string = 'kdiss[{0}, {1}] * y[{0} + {1}]'
        code += ' '.join([ode_string.format(i, comp_nr) for i in range(nr_sizebins)])

    # Mass transfer from water components to biotransformed
    for comp_nr in range(nr_components):
        y_idx = 'y[{0} + {1} + {2}]'.format(nr_components)
        ode_string = 'kbio_waf[{0}] * y[{0}, {1}]'
        code += ' '.join([ode_string.format(i, comp_nr) for i in range(nr_sizebins)])


    # Finish
    code += ']'


def rhs(y, t, kbio_drop, kbio_waf, kdiss):
    """
    y is oil droplet mass for the different droplet sizes

    The size of y is nr_bins * nr_components + nr_components_water + nr_components_water

    y = [bin1_comp1, bin1_comp2, ... bin1_compN, bin2_comp1, ... binM_compN, wat_comp1, ...
         wat_compN, wat_second_metabolit1, ...,  wat_second_metabolitN, bin1_diam, ..., binM_diam]

    k_degr are degradation rates per compoent

    # calculate anew dissoluation rate each timestep.
    # dissolution rate = dissolution_propensity * droplet_bin_area (= nr_droplets * droplet_area)

    k_diss are dissolution rates per component per droplet size

    Can we find a rate for the change in diameter? Yes:
    D = (m * rho * 6 / np.pi) ** (1/3)
    D = m ** 1/3 * K -> dD/dt = d(m**1/3) / dt, here m = sum(y[i:i+nr_components])

    dx*a/dt = a * x ** a - 1

    dD/dt = 1/3 * m ** (-2/3) * K

    The rate for the components is, where di,j means size bin i, component j.
    dmdrop_i,j/dt = -kbio_i,j - kdiss_i,j

    The rate for the water components is
    dmwat_j/dt = sum(kdiss_:,j)

    The rate for the bin diameters is
    dD_i/dt = 1/3 * sum(y[i: i + nr_components]) ** (-2/3) * K,
    where K = ( (rho_i) * 6 / np.pi) ** 1/3
    """


if __name__ == '__main__':
    main()

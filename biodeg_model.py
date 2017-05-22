"""
A model for a a typical biodegradation experiment where droplets are being
mixed and kept in suspension.

Concepts:
* Oil represented through a droplet size distribution
* The size distribution changes as mass is dissolved and biodegraded

The model should capture the following processes for v.0.1:

* Oil dissolution from droplets to water column
* Biodegradation in droplet biofilm
* Biodegradation of dissolved compounds
* Dynamic change of droplet size distribution

The following should be achieved for v.0.2:
* Diffusion of compounds inside the droplet (slower dissolution of certain compounds in larger droplets)

The following should be achieved for v.0.3:
* Growth of a biofilm around the droplets
* Modelling the succession of bacterial species/families
** Scavenger species; species with preferences for different components; etc

Stretch goal for v.0.4:
* Biodegradation rates calculated as a function of microbiological ecosystem

In order to vary certain parameters such as number of components and number of
droplet size bins, it is necessary to generate the ode rhs dynamically. Solve the generation
and execution in the same script as such:

import rhs_code
<generate a new rhs_code.py based on settings>
reload(rhs_code)

"""
import numpy as np
from IPython.terminal.debugger import set_trace as debug  # NOQA
from pint import UnitRegistry
import scipy.integrate as scint
import importlib
import rhs_code
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
    m_droplet(t) = np.exp((-k1 - k2) * t )
    m_water(t) = np.exp(k2 * t)
    m_biotransf_and_biological = np.exp(k1 * t)
    D_droplet = ... something.

    You need to use a solver since the reaction rates keep changing throughout the simulation.
    You also need to use a solver since
    """

    nr_pseudocomponents = 2
    nr_sizebins = 3

    # Droplet sizes (um) and number
    droplet_sizes = np.logspace(1, 7, num=nr_sizebins, base=3.3) * unit('micrometer')
    droplet_number = nr_droplets(droplet_sizes.magnitude, k=2e9)
    droplet_volumes = droplet_volume(droplet_sizes)
    droplet_areas = droplet_volume(droplet_sizes)

    # Rates
    droplet_rates = np.random.random(nr_pseudocomponents) / unit('second')
    waf_rates = droplet_rates * 3 / unit('second')
    dissolution_propensities = np.random.random(nr_pseudocomponents) / unit('second')

    # Mass fraction of each component in each droplet
    composition = np.random.random(nr_pseudocomponents)
    composition = composition / sum(composition)

    # Component density
    rho_min = 700
    rho_max = 1010
    component_densities = rho_max - np.random.random(nr_pseudocomponents) * (rho_max - rho_min)
    component_densities = component_densities * unit('kg') / unit('m') ** 3
    oil_density = np.sum(component_densities * composition)

    # Oil concentration
    oil_conc_ppm = 5 * 1e-6

    # Water volume
    tank_volume = 0.2 * unit('litre')
    water_rho = 1027 * unit('kg') / unit('m') ** 3
    tank_mass = tank_volume * water_rho
    tank_mass = tank_mass.to_base_units()

    # Oil mass (from ppm)
    oil_mass = oil_conc_ppm * tank_mass

    # Need to pre-build array for mass of pseudocomponent per droplet size bin.
    # That is equal to mass per size bin * pseudocompoent fraction
    mass_per_size_class = droplet_volumes * droplet_number * oil_density
    mass_per_size_class = mass_per_size_class.to_base_units()
    sizebin_densities = np.ones(3) * oil_density.magnitude  #

    y0 = setup_yvector(droplet_sizes, mass_per_size_class, composition)

    # Regenerate model
    generate_rhs(nr_sizebins, nr_pseudocomponents, 'rhs_code.py')
    importlib.reload(rhs_code)

    # fix these rates to be of the correct dimensions
    kbio_drop, kbio_waf, kdiss = calc_parameters(droplet_rates, waf_rates,
                                                 dissolution_propensities,
                                                 droplet_areas)
    # Start integrating
    # Load generated module
    r = scint.ode(rhs_code.rhs).set_integrator('vode', method='adams')
    r.set_initial_value(y0, t=0)
    r.set_f_params(kbio_drop, kbio_waf, kdiss, sizebin_densities)

    t1 = 1000
    dt = 1
    while r.successful() and r.t < t1:
        # In between, update the reaction rate parameters and rho
        y_next = r.integrate(r.t + dt)

    debug()


def calc_parameters_kbio_drop(droplet_rates, droplet_areas):
    """
    Do something smarter here in the future
    # kbio_droplets is x rates per n size class [droplet_index, component_index]
    """
    kdrop_bio_pre = []
    for area in droplet_areas:
        kdrop_bio_pre.append(area.magnitude * droplet_rates.magnitude)

    kbio = np.asarray(kdrop_bio_pre)

    return kbio


def calc_parameters_kbio_waf(waf_rates):
    """
    Do something smarter here in the future
    """

    return waf_rates.magnitude


def calc_parameters_dissolution(dissolution_propensities, droplet_areas):
    """
    Do something smarter here in the future
    """
    kdiss_pre = []
    for area in droplet_areas:
        kdiss_pre.append(area.magnitude * dissolution_propensities.magnitude)

    kdiss = np.asarray(kdiss_pre)
    
    return kdiss


def calc_parameters(droplet_rates, waf_rates, dissolution_propensities, droplet_areas):
    """
    Calculate reaction parameters
    """

    # kbio_droplets is x rates per n size class [droplet_index, component_index]
    kbio_droplets = calc_parameters_kbio_drop(droplet_rates, droplet_areas)

    # kbio_waf is one per component
    kbio_waf = calc_parameters_kbio_waf(waf_rates)

    # kdiss is [droplet_index, component_index]
    kdiss = calc_parameters_dissolution(dissolution_propensities, droplet_areas)

    return kbio_droplets, kbio_waf, kdiss


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

    y_initial = []

    # First the mass of each component in each droplet size bin
    for bin_mass in mass_per_size_class:
        for comp_fraction in component_fractions:
            y_initial.append(bin_mass.magnitude * comp_fraction)

    # Then the mass of each component in water (initially 0)
    for i_comp in range(len(component_fractions)):
        y_initial.append(0)

    # Then the diameter of droplets in each bin
    for bin_diameter in bin_diameters:
        y_initial.append(bin_diameter.magnitude)

    return y_initial


def generate_rhs(nr_sizebins, nr_components, filename):
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
    code = 'import numpy as np\n\ndef rhs(t, y, kbio_drop, kbio_waf, kdiss, rho):\n'

    # Start
    code += '\t return ['

    # Droplets: -biodegradation and -dissolution from droplets
    for bin_nr in range(nr_sizebins):
        for comp_nr in range(nr_components):
            y_idx = nr_components * bin_nr + comp_nr
            if comp_nr > 0:
                code += '\t'
            code += '(-kbio_drop[{0}, {1}] - kdiss[{0}, {1}]) * y[{2}]'.format(bin_nr, comp_nr, y_idx)
            code += ',\n'

    # Dissolved mass: +dissolution -biodegradation
    for comp_nr in range(nr_components):
        # Dissolution
        growth = ''
        for bin_nr in range(nr_sizebins):
            y_idx_droplet = nr_components * bin_nr + comp_nr
            growth += ' + kdiss[{0}, {1}] * y[{2}]'.format(bin_nr, comp_nr, y_idx_droplet)

        # Biodegradation
        y_idx_dissolved = nr_components * nr_sizebins + comp_nr
        ode_string_biodeg = ' - kbio_waf[{0}] * y[{1}]'
        loss = ode_string_biodeg.format(comp_nr, y_idx_dissolved)
        code += '\t'
        code += growth + loss
        code += ',\n'

    # Droplet diameter: change with mass
    for bin_nr in range(nr_sizebins):

        start = bin_nr * nr_components
        end = start + nr_components
        ode_string = '1 / 3 * (sum(y[{0}:{1}]) ** (-2 / 3)) * (rho[{2}] * 6 / np.pi) ** 1 / 3'
        code += '\t'
        code += ode_string.format(start, end, bin_nr)
        code += ',\n'

    # Finish
    code += '\t]'

    with open(filename, 'wt') as f:
        f.write(code)

    return code


if __name__ == '__main__':
    main()

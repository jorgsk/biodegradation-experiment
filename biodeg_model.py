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

    return (V * 6 / np.pi) ** (1/3)


def main():

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


def dissolution_rate(droplet_masses, droplet_densities):
    """
    Input: droplet masses per size bin and densities.

    Calculate the area of each droplet, find rate by Area * dissolution_propensity.
    """


def rhs(y, t, k_degr, k_diss):
    """
    y is oil droplet mass for the different droplet sizes

    The size of y is nr_bins * nr_components + nr_components_water + nr_components_water

    y = [bin1_comp1, bin1_comp2, ... bin1_compN, bin2_comp1, ... binM_compN, wat_comp1, ...
         wat_compN, wat_second_metabolit1, ...,  wat_second_metabolitN, bin1_diam, ..., binM_diam]

    k_degr are degradation rates per compoent

    # calculate anew dissoluation rate each timestep.
    # dissolution rate = dissolution_propensity * droplet_bin_area (= nr_droplets * droplet_area)

    k_diss are dissolution rates per component per droplet size

    diameters are shrunk each timestep

    you lose some mass -> some volume

    Problem: Sphere with (D, V) now as V2. What is D2? Just go in reverse, no? f(V) = D

    for drop_bin in nr_binsizes():
    drop_mass = y[i: i + nr_components]
    drop_vol = drop_mass / drop_rho
    d_new = d_old - 

    Can we find a rate for the change in diameter?
    """


if __name__ == '__main__':
    main()

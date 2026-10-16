"""Absolute particle and energy bookkeeping for a directional emission model.

This does not validate a field or fit. It tests whether an already fitted
directional photon yield would require more power than the pulsar supplies.
The geometry must already include both poles with globally normalized seed
weights. No isotropic 4*pi factor belongs in the directional normalization.
"""
import numpy as np

PC_CM = 3.0856775814913673e18
C_CM_S = 29979245800.


def field_scale_from_source_power(period_s, spin_down_erg_s, dimensionless_source_power):
    """Common Gaussian E/B scale using a measured numerical source power.

    The input power must be the time-averaged outward integral of (E x B)
    dot dA/(4*pi) in c=RLC=1 units, measured outside the source and before
    artificial absorption. It must not be an unverified distant flux after
    numerical dissipation. L_phys=c*RLC^2*B_scale^2*L_code.
    Source stability and grid/domain/source-size sensitivity remain required.
    This replaces the point-dipole torque approximation only when those
    measurements exist; it does not remove the assumed moment of inertia
    used to obtain the observed spin-down luminosity.
    """
    values=np.asarray([period_s,spin_down_erg_s,dimensionless_source_power],float)
    if np.any(~np.isfinite(values)) or np.any(values<=0):
        raise ValueError('Positive finite period, spin-down and numerical source power required')
    rlc=C_CM_S*period_s/(2*np.pi)
    return float(np.sqrt(spin_down_erg_s/(C_CM_S*rlc**2*dimensionless_source_power)))


def population_power_budget(normalization, distance_pc, spin_down_erg_s,
                            radiation_erg_per_primary, *,
                            work_erg_per_primary=None,
                            injected_kinetic_erg_per_primary=None,
                            escaping_kinetic_erg_per_primary=None):
    """Convert F(E)=normalization*Y(E)/sr to a steady primary supply.

    Y is the phase-averaged photon yield per injected primary per steradian;
    F is photons/(cm2 s). Thus normalization=Ndot/d^2, NOT Ndot/(4*pi*d^2).
    The normalization is a single spectrum-wide fit, never one per band.
    Work may be signed. Injection/escape energies exclude rest mass. All
    supplied energies must be integrated over every modeled direction.

    A sub-spin-down budget is necessary, not sufficient. The calculation
    does not infer pair multiplicity or include unseen photon components,
    pair creation/rest energy, heating or the wind beyond the trace boundary.
    """
    required = np.asarray([normalization, distance_pc, spin_down_erg_s,
                           radiation_erg_per_primary], float)
    if (np.any(~np.isfinite(required)) or normalization < 0 or distance_pc <= 0
            or spin_down_erg_s <= 0 or radiation_erg_per_primary < 0):
        raise ValueError('Require finite normalization/radiation >=0 and distance/power >0')
    optional = [work_erg_per_primary, injected_kinetic_erg_per_primary,
                escaping_kinetic_erg_per_primary]
    for i, value in enumerate(optional):
        if value is not None and (not np.isfinite(value) or (i > 0 and value < 0)):
            raise ValueError('Work must be finite; injection/escape kinetic energy nonnegative')
    rate = float(normalization * (distance_pc * PC_CM)**2)
    output = dict(distance_pc=float(distance_pc), normalization=float(normalization),
                  primaries_per_second=rate,
                  spin_down_erg_s=float(spin_down_erg_s),
                  radiation_erg_per_primary=float(radiation_erg_per_primary),
                  radiation_power_erg_s=rate*float(radiation_erg_per_primary),
                  radiation_fraction_of_spin_down=rate*float(radiation_erg_per_primary)/spin_down_erg_s,
                  normalization_convention='F=Ndot*Y_per_primary_per_sr/distance^2; both poles already included; no 4*pi',
                  budget_complete=False, vela_validation=False,
                  caveat='Necessary energy check only; no pair supply solution, rest-mass creation cost, unmodeled radiation, heating or complete wind budget')
    for key, value in zip(('work', 'injected_kinetic', 'escaping_kinetic'), optional):
        if value is not None:
            output[key+'_erg_per_primary'] = float(value)
            output[key+'_power_erg_s'] = rate*float(value)
    if work_erg_per_primary is not None and injected_kinetic_erg_per_primary is not None:
        supplied = rate*(work_erg_per_primary+injected_kinetic_erg_per_primary)
        output['modeled_supply_power_erg_s'] = supplied
        output['modeled_supply_fraction_of_spin_down'] = supplied/spin_down_erg_s
    if all(value is not None for value in optional):
        residual = (injected_kinetic_erg_per_primary + work_erg_per_primary
                    - escaping_kinetic_erg_per_primary - radiation_erg_per_primary)
        output['energy_ledger_residual_erg_per_primary'] = float(residual)
    output['radiation_exceeds_spin_down'] = output['radiation_fraction_of_spin_down'] > 1
    if 'modeled_supply_fraction_of_spin_down' in output:
        output['modeled_supply_exceeds_spin_down'] = output['modeled_supply_fraction_of_spin_down'] > 1
    return output

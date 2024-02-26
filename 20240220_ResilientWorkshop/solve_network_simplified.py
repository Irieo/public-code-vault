# -*- coding: utf-8 -*-
# IR 2024-02-20 simplified version of solve_network.py

import logging
import numpy as np
import pandas as pd
import pypsa
import xarray as xr
from _benchmark import memory_logger

logger = logging.getLogger(__name__)
pypsa.pf.logger.setLevel(logging.WARNING)


def add_battery_constraints(n):
    """
    Add constraint ensuring that charger = discharger, i.e.
    1 * charger_size - efficiency * discharger_size = 0
    """
    if not n.links.p_nom_extendable.any():
        return

    discharger_bool = n.links.index.str.contains("battery discharger")
    charger_bool = n.links.index.str.contains("battery charger")

    dischargers_ext = n.links[discharger_bool].query("p_nom_extendable").index
    chargers_ext = n.links[charger_bool].query("p_nom_extendable").index

    eff = n.links.efficiency[dischargers_ext].values
    lhs = (
        n.model["Link-p_nom"].loc[chargers_ext]
        - n.model["Link-p_nom"].loc[dischargers_ext] * eff
    )

    n.model.add_constraints(lhs == 0, name="Link-charger_ratio")


def NECP_constraint(n, config) -> None:
    """
    Set a minimum share of RES based on NECP.
    for simplicity onwind and solar only
    """

    res_techs = ["onwind", "solar"]
    weights = n.snapshot_weightings["generators"]

    country_buses = n.buses.index[(n.buses.index.str[:2] == "DE")]

    country_loads = n.loads.index[n.loads.bus.isin(country_buses)]
    country_res_gens = n.generators.index[
        n.generators.bus.isin(country_buses) & n.generators.carrier.isin(res_techs)
    ]

    gens = n.model["Generator-p"].loc[:, country_res_gens] * weights
    lhs = gens.sum()

    target = snakemake.config["DE_NECP"]["RES_share"] / 100
    total_load = (n.loads_t.p_set[country_loads].sum(axis=1) * weights).sum()

    print(
        f"RES constraint for DE is {target} and total load {round(total_load/1e6, 2)} TWh"
    )

    n.model.add_constraints(lhs >= target * total_load, name=f"DE_res_constraint")


def no_coal_policy(n: pypsa.Network) -> None:
    """
    remove "coal" and "lignite" power plant fleet

    Args:
    - n: The network object to be modified.

    Returns:
    - None
    """
    n.generators.loc[n.generators.index.str.contains("coal|lignite"), "p_nom"] = 0


def extra_functionality(n, config):
    """
    Add extra constraints to the mathematical model.
    """
    # Add battery constraints
    add_battery_constraints(n)

    # conditionally add RES constraint
    if config["DE_NECP"]["switcher"]:
        NECP_constraint(n, config)


def solve_network(n, config, log_fn=None):
    n.optimize.create_model()

    extra_functionality(n, config)

    n.optimize.solve_model(
        solver_name="gurobi",
        log_fn=snakemake.log.solver,
    )

    return n


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "solve_sector_network",
            configfiles="config/config-workshop.yaml",
            simpl="",
            opts="3H",
            clusters="1",
            ll="vopt",
            sector_opts="",
            planning_horizons="2030",
        )

    n = pypsa.Network("../resources/example-workshop/networks/elec_s_1_ec_lvopt_3H.nc")

    # small_horison(n, nhours=240)

    if snakemake.config["no_coal_policy"]:
        no_coal_policy(n)

    with memory_logger(
        filename=getattr(snakemake.log, "memory", None), interval=30.0
    ) as mem:
        n = solve_network(
            n,
            config=snakemake.config,
            log_fn=snakemake.log.solver,
        )

    # n.export_to_netcdf(snakemake.output[0])

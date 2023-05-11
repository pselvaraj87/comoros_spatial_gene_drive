#!/usr/bin/env python3

import pathlib  # for a join

# idmtools ...
from idmtools.builders import SimulationBuilder
from idmtools.core.platform_factory import Platform
from idmtools.entities.experiment import Experiment
from idmtools.utils.filter_simulations import FilterItem

# emodpy
from emodpy.emod_task import EMODTask
from emodpy_malaria.reporters.builtin import ReportVectorGenetics, ReportVectorMigration

from comoros_exploration.helpers import *
import comoros_exploration.params as params
import comoros_exploration.manifest as manifest

from functools import partial
import pandas as pd


def get_serialization_paths(platform, serialization_exp_id):
    exp = Experiment.from_id(serialization_exp_id, children=False)
    exp.simulations = platform.get_children(exp.id, exp.item_type,
                                            children=["tags", "configuration", "files", "hpc_jobs"])

    sim_dict = {'Larval_Capacity': [], 'Outpath': []}
    for simulation in exp.simulations:
        # if simulation.tags['Run_Number'] == 0:
        string = simulation.get_platform_object().hpc_jobs[0].working_directory.replace('internal.idm.ctr', 'mnt')
        string = string.replace('\\', '/')
        string = string.replace('IDM2', 'idm2')

        sim_dict['Larval_Capacity'] += [float(simulation.tags['Larval_Capacity'])]
        sim_dict['Outpath'] += [string]

    df = pd.DataFrame(sim_dict)
    return df


def general_sim(serialization=0, serialized_exp_id=None):
    """
    This function is designed to be a parameterized version of the sequence of things we do
    every time we run an emod experiment.
    """

    # Create a platform
    # Show how to dynamically set priority and node_group
    platform = Platform("SLURM")

    # create EMODTask
    print("Creating EMODTask (from files)...")

    task = EMODTask.from_default2(
        config_path="my_config.json",
        eradication_path=manifest.eradication_path,
        ep4_custom_cb=None,
        campaign_builder=None,
        schema_path=manifest.schema_file,
        param_custom_cb=set_param_fn,
        demog_builder=build_demographics
    )
    task.set_sif(manifest.sif_path)

    # Create simulation sweep with builder
    builder = SimulationBuilder()

    # Add asset
    task.common_assets.add_asset("/Users/pselvaraj/Github/testing/emodpy-vector_genetics/download/schema.json")
    task.common_assets.add_asset("/Users/pselvaraj/Github/testing/emodpy-vector_genetics/comoros_exploration/vector_migration.bin")
    task.common_assets.add_asset(
        "/Users/pselvaraj/Github/testing/emodpy-vector_genetics/comoros_exploration/vector_migration.bin.json")

    if serialized_exp_id:
        serialized_population_path_df = get_serialization_paths(platform=platform,
                                                                serialization_exp_id=serialized_exp_id)

        # Sweep larval capacity
        func = partial(update_serialize, serialization=serialization, sim_duration=6 * 365,
                       serialized_population_path_df=serialized_population_path_df)
        # builder.add_sweep_definition(func, [8.15, 8.45, 8.65])
        builder.add_sweep_definition(func, [8.25])

        # Sweep run number
        builder.add_sweep_definition(update_sim_random_seed, range(params.nSims))

        # Sweep campaign type
        func = partial(update_camp_type, serialize=serialization, sim_duration=6 * 365)
        builder.add_sweep_definition(func, [True, False])

        exp_name = params.exp_name

        # Sweep transmission to human
        builder.add_sweep_definition(update_transmission_to_human, [0.2])
        #
        # Sweep infected progress
        builder.add_sweep_definition(update_infected_progress, [0.3])

        # Sweep mortality
        # builder.add_sweep_definition(update_mortality, [1.1, 1.2, 1.3, 1.4])
        builder.add_sweep_definition(update_mortality, [1.05])


        # # Single roundtrips
        #     config.parameters.Enable_Regional_Migration = 1
        #     DT._set_migration_pattern_srt(config)
        #     DT._set_regional_migration_filenames(config, "regional_migration.bin")
        #     DT._set_regional_migration_roundtrip_probability(config, 1)Add reporter - vector genetics allele frequency
        reporter = ReportVectorGenetics()  # Create the reporter
        reporter.config(rvg_config_builder, manifest)  # Config the reporter
        task.reporters.add_reporter(reporter)  # Add thre reporter

        # # Add vector genome reporter
        # reporter = ReportVectorGenetics()  # Create the reporter
        # reporter.config(rvg_config_builder_genome, manifest)  # Config the reporter
        # task.reporters.add_reporter(reporter)  # Add thre reporter
        # reporter = ReportVectorMigration()
        # reporter.config(rvmg_config_builder, manifest)
        # task.reporters.add_reporter(reporter)

    else:
        func = partial(update_serialize, serialization=serialization, sim_duration=40 * 365,
                       serialized_population_path_df=None)
        builder.add_sweep_definition(func, [8, 8.25])
        builder.add_sweep_definition(update_sim_random_seed, [0])
        func = partial(update_camp_type, serialize=serialization, sim_duration=40 * 365)
        builder.add_sweep_definition(func, [1])
        exp_name = params.exp_name + '_serialization'

        # reporter = ReportVectorMigration()
        # reporter.config(rvmg_config_builder, manifest)
        # task.reporters.add_reporter(reporter)

    # create experiment from builder
    print(f"Prompting for COMPS creds if necessary...")
    experiment = Experiment.from_builder(builder, task, name=exp_name)

    # The last step is to call run() on the ExperimentManager to run the simulations.
    experiment.run(wait_until_done=True, platform=platform)

    # Check result
    if not experiment.succeeded:
        print(f"Experiment {experiment.uid} failed.\n")
        exit()

    print(f"Experiment {experiment.uid} succeeded.")

    # Save experiment id to file
    with open("COMPS_ID", "w") as fd:
        fd.write(experiment.uid.hex)
    print()
    print(experiment.uid.hex)


if __name__ == "__main__":
    # TBD: user should be allowed to specify (override default) erad_path and input_path from command line
    # plan = EradicationBambooBuilds.MALARIA_LINUX
    # print("Retrieving Eradication and schema.json from Bamboo...")
    # get_model_files( plan, manifest )
    # print("...done.")

    serialization = 0
    serialization_experiment_id = 'c9f344d2-7953-ed11-a9ff-b88303911bc1'
    # serialization = 1
    # serialization_experiment_id = None
    general_sim(serialization=serialization, serialized_exp_id=serialization_experiment_id)

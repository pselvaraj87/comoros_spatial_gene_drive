import os
import json
import pandas as pd
from functools import \
    partial  # for setting Run_Number. In Jonathan Future World, Run_Number is set by dtk_pre_proc based on generic param_sweep_value...
import sporozoite_delay.manifest as manifest

import emodpy_malaria.demographics.MalariaDemographics as Demographics  # OK to call into emod-api
from emod_api.demographics import DemographicsTemplates as DT
from helpers_demographics import *

from emodpy_malaria import malaria_config as malconf
from emodpy_malaria import vector_config as vecconf
from emodpy_malaria.interventions.treatment_seeking import add_treatment_seeking
from emodpy_malaria.interventions.mosquitorelease import add_scheduled_mosquito_release
import emod_api.migration as migration
from emod_api.migration.vector_migration import Migration as vector_Migration
from emod_api.migration.vector_migration import from_demog_and_param_gravity as vector_from_demog_and_param_gravity

import emod_api.campaign as camp


def build_demographics():

    # YOU NEED TO BE ON VPN (IDM internal network) to be able to access this server
    # input_file = '/Users/pselvaraj/Github/testing/emodpy-vector_genetics/comoros_exploration/data/comoros_population_com_2019-07-01.csv'
    input_file = '/Users/pselvaraj/Github/testing/emodpy-vector_genetics/comoros_exploration/data/comoros_population_com_temp.csv'
    demographics = Demographics.from_pop_csv(input_file, site='Comoros')
    birthrate_file = '/Users/pselvaraj/Github/testing/emodpy-vector_genetics/comoros_exploration/data/WB_birth_rates.csv'
    birthrate_df = pd.read_csv(birthrate_file)
    birthrate = float(birthrate_df[birthrate_df['Country Name'] == 'Comoros']['2020 [YR2020]'].values[0])

    DT.MortalityRateByAge(demographics, age_bins=[0], mort_rates=[birthrate])
    DT.InitPrevUniform(demographics, low_prev=0.13, high_prev=0.15)
    DT.FullMigrationHeterogeneity(demographics)

    migration_partial = partial(migration.from_demog_and_param_gravity,
                                gravity_params=[7.50395776e-06, 9.65648371e-01, 9.65648371e-01, -1.10305489e+00],
                                id_ref='comoros', migration_type=migration.Migration.REGIONAL)

    return demographics, migration_partial



def find_genome_index_in_trait_modifiers(vsp, genome, trait_of_interest):
    flat_gen = [i for j in genome for i in j]
    flat_gen.sort()

    trait_modifiers = vsp['Gene_To_Trait_Modifiers']

    for i, t in enumerate(trait_modifiers):

        config_gen = [m for j in t['Allele_Combinations'] for m in j]
        config_gen.sort()
        traits = [j['Trait'] for j in t['Trait_Modifiers']]
        t_index = None
        if trait_of_interest in traits:
            t_index = traits.index(trait_of_interest)

        if (t_index is not None) and (config_gen == flat_gen):
            return i, t_index

    return None


def update_sim_bic(simulation, value):
    simulation.task.config.parameters.Base_Infectivity_Constant = value * 0.1
    return {"Base_Infectivity": value}


def update_sim_random_seed(simulation, value):
    simulation.task.config.parameters.Run_Number = value
    return {"Run_Number": value}


def update_sim_larval_capacity(simulation, value):
    simulation.task.config.parameters.Vector_Species_Params[0].Habitats[0]['Max_Larval_Capacity'] = value
    return {"Larval_Capacity": value}


def update_mortality(simulation, value):
    vsp = simulation.task.config.parameters.Vector_Species_Params[0]

    genome = [["b1", "b1"]]
    trait = "MORTALITY"
    idx1, idx2 = find_genome_index_in_trait_modifiers(vsp=vsp, genome=genome, trait_of_interest=trait)
    vsp.Gene_To_Trait_Modifiers[idx1].Trait_Modifiers[idx2].Modifier = value

    return {"Mortality": value}


def update_infected_progress(simulation, value):
    vsp = simulation.task.config.parameters.Vector_Species_Params[0]

    genome = [["X", "X"], ["b1", "b1"]]
    trait = "INFECTED_PROGRESS"
    idx1, idx2 = find_genome_index_in_trait_modifiers(vsp=vsp, genome=genome, trait_of_interest=trait)
    vsp.Gene_To_Trait_Modifiers[idx1].Trait_Modifiers[idx2].Modifier = value

    # genome = [["X", "X"], ["b1", "b0"]]
    # trait = "INFECTED_PROGRESS"
    # idx1, idx2 = find_genome_index_in_trait_modifiers(vsp=vsp, genome=genome, trait_of_interest=trait)
    # vsp.Gene_To_Trait_Modifiers[idx1].Trait_Modifiers[idx2].Modifier = 1-(1-value)/2

    # genome = [["X", "X"], ["b1", "b2"]]
    # trait = "INFECTED_PROGRESS"
    # idx1, idx2 = find_genome_index_in_trait_modifiers(vsp=vsp, genome=genome, trait_of_interest=trait)
    # vsp.Gene_To_Trait_Modifiers[idx1].Trait_Modifiers[idx2].Modifier = 1-(1-value)/2
    #
    # genome = [["X", "X"], ["b1", "b3"]]
    # trait = "INFECTED_PROGRESS"
    # idx1, idx2 = find_genome_index_in_trait_modifiers(vsp=vsp, genome=genome, trait_of_interest=trait)
    # vsp.Gene_To_Trait_Modifiers[idx1].Trait_Modifiers[idx2].Modifier = 1-(1-value)/2

    return {"Infected_Progress": value}


def update_transmission_to_human(simulation, value):
    vsp = simulation.task.config.parameters.Vector_Species_Params[0]

    genome = [["X", "X"], ["b1", "b1"]]
    trait = "TRANSMISSION_TO_HUMAN"
    idx1, idx2 = find_genome_index_in_trait_modifiers(vsp=vsp, genome=genome, trait_of_interest=trait)
    vsp.Gene_To_Trait_Modifiers[idx1].Trait_Modifiers[idx2].Modifier = value

    # genome = [["X", "X"], ["b1", "b0"]]
    # trait = "TRANSMISSION_TO_HUMAN"
    # idx1, idx2 = find_genome_index_in_trait_modifiers(vsp=vsp, genome=genome, trait_of_interest=trait)
    # vsp.Gene_To_Trait_Modifiers[idx1].Trait_Modifiers[idx2].Modifier = 1 - (1 - value) / 2
    #
    # genome = [["X", "X"], ["b1", "b2"]]
    # trait = "TRANSMISSION_TO_HUMAN"
    # idx1, idx2 = find_genome_index_in_trait_modifiers(vsp=vsp, genome=genome, trait_of_interest=trait)
    # vsp.Gene_To_Trait_Modifiers[idx1].Trait_Modifiers[idx2].Modifier = 1 - (1 - value) / 2
    #
    # genome = [["X", "X"], ["b1", "b3"]]
    # trait = "TRANSMISSION_TO_HUMAN"
    # idx1, idx2 = find_genome_index_in_trait_modifiers(vsp=vsp, genome=genome, trait_of_interest=trait)
    # vsp.Gene_To_Trait_Modifiers[idx1].Trait_Modifiers[idx2].Modifier = 1 - (1 - value) / 2

    return {"Transmission_To_Human": value}


def update_serialize(simulation, larval_multiplier, serialization=0, sim_duration=40 * 365,
                     serialized_population_path_df=None):
    if serialization:
        simulation.task.config.parameters.Simulation_Duration = sim_duration
        simulation.task.config.parameters.Serialization_Time_Steps = [sim_duration]
        simulation.task.config.parameters.Serialized_Population_Reading_Type = 'NONE'
        simulation.task.config.parameters.Serialized_Population_Writing_Type = 'TIMESTEP'
        update_sim_larval_capacity(simulation, value=pow(10, larval_multiplier))

    else:
        serialized_population_path = serialized_population_path_df[serialized_population_path_df['Larval_Capacity']
                                                                   == larval_multiplier]['Outpath'].values[0]
        simulation.task.config.parameters.Simulation_Duration = sim_duration
        simulation.task.config.parameters.Serialization_Mask_Node_Read = 0
        simulation.task.config.parameters.Serialization_Mask_Node_Write = 0
        # simulation.task.config.parameters.Serialization_Type = 'NONE'
        simulation.task.config.parameters.Serialized_Population_Path = os.path.join(serialized_population_path,
                                                                                    'output')
        # simulation.task.config.parameters.Serialized_Population_Path = '//mnt/idm2/home/pselvaraj/output/sporozoite_delay_gene_drive_seriali_20210705_023006/c8d/f21/f23/c8df21f2-38dd-eb11-a9ec-b88303911bc1/output'
        simulation.task.config.parameters.Serialized_Population_Reading_Type = 'READ'
        simulation.task.config.parameters.Serialized_Population_Writing_Type = 'NONE'
        simulation.task.config.parameters.Serialized_Population_Filenames = ['state-14600-%03d.dtk' % i for i in range(48)]
        update_sim_larval_capacity(simulation, value=pow(10, larval_multiplier))

    return {"Serialization": serialization, "Larval_Capacity": larval_multiplier}


def set_param_fn(config):
    """
    This function is a callback that is passed to emod-api.config to set parameters The Right Way.
    """
    config = malconf.set_team_defaults(config, manifest)
    # config = set_config.set_config(config)

    config.parameters.Age_Initialization_Distribution_Type = "DISTRIBUTION_OFF"
    config.parameters.Base_Rainfall = 150
    config.parameters.Climate_Model = "CLIMATE_CONSTANT"
    config.parameters.Enable_Disease_Mortality = 0
    config.parameters.Enable_Demographics_Risk = 0
    config.parameters.Enable_Vector_Species_Report = 0
    config.parameters.Enable_Initial_Prevalence = 1
    config.parameters.Enable_Vector_Aging = 1
    config.parameters.Enable_Natural_Mortality = 1
    config.parameters.Demographics_Filenames = ['demographics.json']
    config.parameters.Enable_Demographics_Birth = 1
    # config.parameters.pop("Serialized_Population_Filenames")

    config.parameters.Enable_Spatial_Output = 1
    config.parameters.Spatial_Output_Channels = ['Prevalence', 'Population']

    # Single roundtrips
    config.parameters.Enable_Regional_Migration = 1
    DT._set_migration_pattern_srt(config)
    DT._set_regional_migration_filenames(config, "regional_migration.bin")
    DT._set_regional_migration_roundtrip_probability(config, 1)

    # Vector migration
    config.parameters.Enable_Vector_Migration = 1
    config.parameters.Enable_Vector_Migration_Regional = 1
    config.parameters.Vector_Migration_Filename_Regional = 'vector_migration.bin'
    config.parameters.Vector_Sampling_Type = "TRACK_ALL_VECTORS"
    config.parameters.x_Vector_Migration_Regional = 10

    # --------- Define integral gene drive functions
    # - Initial resistance frequencies (defined in set_integral_genes)
    # rr10 = 0  # initial frequency of resistance (locus 1)
    # rr20 = 0  # initial frequency of resistance (locus 2)
    # - Mosquito release params (defined in add_integral_release)
    # NOTE: Compared to Table 1 in Nash et al. (2019), cc0, dd0, & ee0 are #s, not freqs
    # cc0 = 100  # initial release NUMBER of homozygous drive and effector construct
    # dd0 = 0  # initial release NUMBER of homozygous drive construct
    # ee0 = 0  # initial release NUMBER of homozygous effector construct
    # - Gene drive params (defined in add_integral_gene_drives)
    # d1 = 0.99  # transmission rate of drive
    # p_nhej = 0.5  # prob of NHEJ at each locus, given no drive
    # p_r_nhej = 1 / 3  # prob resistance arising from NHEJ at each locus
    # p_ihdr = 1e-4  # prob of incomplete HDR at each locus, given drive
    # p_r_ihdr = 1 / 3  # prob of resistance arising from incomplete HDR at each locus
    # - Fitness params (defined in add_integral_fitness_costs)
    # NOTE: Can check Fig. 2/IGD - Drive Plot (Locus 1) 2 Locus.nb from the Nash paper github
    # for accuracy of the below fitnesses; see section titled Solve Differential Equations Model parameters
    # sd1 = 0  # cost of hijacking target locus 1 (drive)
    # hd1 = 0.5  # dominance coefficient for hijacking at locus 1
    # sd2 = 0  # cost of hijacking target locus 2 (effector)
    # hd2 = 0.5  # dominance coefficient for hijacking at locus 2
    # sn = 0.05  # cost of expressing nuclease at locus 1
    # hn = 0.5  # dominance coefficient for expressing nuclease
    # se2 = 0.1  # cost of expressing effector at locus 2
    # he2 = 0.5  # dominance coefficient for expressing effector (locus 2)
    # sm = 1  # cost of loss of gene function
    # hm = 0.2  # dominance coefficient for loss of gene function
    # - Refractoriness (defined in add_integral_fitness_costs)
    # rc = 1  # homozygous degree of refractoriness
    # hrc1 = 1  # dominance coefficient for refractoriness (one effector allele)

    # Vector Genetics
    malconf.add_species(config, manifest, ['gambiae'])
    habitats = vecconf.configure_linear_spline(manifest, max_larval_capacity=pow(10, 8),
                                               capacity_distribution_number_of_years=1,
                                               capacity_distribution_over_time={
                                                                   "Times": [0.0, 30.417, 60.833, 91.25,
                                                                             121.667, 152.083,
                                                                             182.5, 212.917, 243.333,
                                                                             273.75, 304.167, 334.583],
                                                                   "Values": [i/35 for i in [3, 0.8, 1.25, 0.1, 2.7, 10, 6,
                                                                              35, 2.8, 1.5, 1.6, 2.1]]
                                                               })
    vecconf.set_species_param(config, 'gambiae', 'Habitats', [habitats], overwrite=True)

    vecconf.add_genes_and_alleles(config, manifest, 'gambiae', [('a0', 1.0), ('a1', 0), ('a2', 0), ('a3', 0)])
    vecconf.add_genes_and_alleles(config, manifest, 'gambiae', [('b0', 1.0), ('b1', 0), ('b2', 0), ('b3', 0)])

    # vecconf.add_trait(config, manifest, [["X", "X"], ["b1", "b1"]], "INFECTED_PROGRESS", 0.5)
    infected_progress = ("INFECTED_PROGRESS", 0.3)
    transmission_to_human = ("TRANSMISSION_TO_HUMAN", 0.2)
    vecconf.add_trait(config, manifest, 'gambiae', [["X", "X"], ["b1", "b1"]],
                      [transmission_to_human, infected_progress])

    # Fitness cost  parameters
    sd1 = 0
    hd1 = 0.5
    sd2 = 0
    hd2 = 0.5
    sn = 0.05
    hn = 0.5
    se2 = 0.1
    he2 = 0.5
    sm = 1
    hm = 0.2
    max_mf = 100
    hfrn = 0.5
    frn = 0
    hfre2 = 0.5
    fre2 = 0
    fre3 = 0.15  # Including fitness costs for homozygotes with effector only

    # Driver fitness
    vecconf.add_trait(config, manifest, 'gambiae', [["a0", "a1"]], [("MORTALITY",
                                                                     1 / ((1 - hd1 * sd1) * (1 - hn * sn))),
                                                                    ("FECUNDITY",
                                                                     1 - hfrn * frn)])
    vecconf.add_trait(config, manifest, 'gambiae', [["a0", "a3"]], [("MORTALITY",
                                                                     1 / (1 - hm * sm))])
    vecconf.add_trait(config, manifest, 'gambiae', [["a1", "a1"]], [("MORTALITY",
                                                                     1 / ((1 - sd1) * (1 - sn))),
                                                                    ("FECUNDITY",
                                                                     1 - frn)])
    vecconf.add_trait(config, manifest, 'gambiae', [["a1", "a2"]], [("MORTALITY",
                                                                     1 / ((1 - hd1 * sd1) * (1 - hn * sn))),
                                                                    ("FECUNDITY",
                                                                     1 - hfrn * frn)])
    vecconf.add_trait(config, manifest, 'gambiae', [["a1", "a3"]], [("MORTALITY",
                                                                     1 / ((1 - hd1 * sd1) * (1 - hn * sn) * (
                                                                             1 - hm * sm))),
                                                                    ("FECUNDITY",
                                                                     1 - hfrn * frn)])
    vecconf.add_trait(config, manifest, 'gambiae', [["a2", "a3"]],
                      [("MORTALITY", 1 / (1 - hm * sm))])
    vecconf.add_trait(config, manifest, 'gambiae', [["a3", "a3"]], [("MORTALITY", max_mf)])

    # Effector fitness
    vecconf.add_trait(config, manifest, 'gambiae', [["b0", "b3"]], [("MORTALITY", 1 / (1 - hm * sm))])
    vecconf.add_trait(config, manifest, 'gambiae', [["b1", "b1"]], [("MORTALITY", 1 / ((1 - sd2) * (1 - se2))),
                                                                    ("FECUNDITY", 1 - fre3)])
    vecconf.add_trait(config, manifest, 'gambiae', [["b1", "b3"]], [("MORTALITY", 1 / (
            (1 - hd2 * sd2) * (1 - he2 * se2) * (
            1 - hm * sm)))])
    vecconf.add_trait(config, manifest, 'gambiae', [["b2", "b3"]], [("MORTALITY", 1 / (1 - hm * sm))])
    vecconf.add_trait(config, manifest, 'gambiae', [["b3", "b3"]], [("MORTALITY", max_mf)])

    # Gene drive - using full dictionary here but hope to change to something more elegant in the future
    d1 = 0.95
    p_nhej = 0.5
    p_ihdr = 1e-4
    p_r_nhej = 1 / 3
    p_r_ihdr = 1 / 3
    a1_ctl = [("a0", (1 - d1) * (1 - p_nhej)), ("a1", d1 * (1 - p_ihdr)),
              ("a2", (1 - d1) * p_r_nhej * p_nhej + d1 * p_r_ihdr * p_ihdr),
              ("a3", (1 - d1) * (1 - p_r_nhej) * p_nhej + d1 * (1 - p_r_ihdr) * p_ihdr)]
    b1_ctl = [("b0", (1 - d1) * (1 - p_nhej)), ("b1", d1 * (1 - p_ihdr)),
              ("b2", (1 - d1) * p_r_nhej * p_nhej + d1 * p_r_ihdr * p_ihdr),
              ("b3", (1 - d1) * (1 - p_r_nhej) * p_nhej + d1 * (1 - p_r_ihdr) * p_ihdr)]

    vecconf.add_species_drivers(config, manifest, species='gambiae', driving_allele='a1',
                                driver_type="INTEGRAL_AUTONOMOUS", to_copy='a1', to_replace='a0',
                                likelihood_list=a1_ctl)
    vecconf.add_species_drivers(config, manifest, species='gambiae', driving_allele='a1',
                                driver_type="INTEGRAL_AUTONOMOUS", to_copy='b1', to_replace='b0',
                                likelihood_list=b1_ctl)

    # # Crete new Vector Species Params
    # config = set_vsp(config, manifest)
    return config


def update_camp_type(simulation, baseline, serialize=0, sim_duration=40 * 365):
    # simulation.task.config.parameters.Run_Number = value
    build_camp_partial = partial(build_camp, serialize=serialize, sim_duration=sim_duration, baseline=baseline)
    simulation.task.create_campaign_from_callback(build_camp_partial)

    return {"Baseline": baseline}


def build_camp(baseline=0, serialize=0, sim_duration=40 * 365):
    """
    Build a campaign input file for the DTK using emod_api.
    Right now this function creates the file and returns the filename. If calling code just needs an asset that's fine.
    """

    # This isn't desirable. Need to think about right way to provide schema (once)
    camp.set_schema(manifest.schema_file)

    # print( f"Telling emod-api to use {manifest.schema_file} as schema." )
    # camp.add(bednet.Bednet(camp, start_day=start_day_in, coverage=0.5, killing_eff=0.5, blocking_eff=0.5, usage_eff=0.5,
    #                        insecticide="pyrethroid"), first=True)
    if not serialize:
        add_treatment_seeking(camp, targets=[{"trigger": "NewClinicalCase", "coverage": 0.8, "agemin": 15,
                                              "agemax": 70, "seek": 0.4, "rate": 0.3},
                                             {"trigger": "NewSevereCase", "coverage": 0.8, "seek": 0.6, "rate": 0.5}],
                              drug=['Artemether', 'Lumefantrine'],
                              start_day=0,
                              broadcast_event_name='Received_Treatment'
                              )
        if not baseline:
            add_scheduled_mosquito_release(camp, start_day=180, released_number=1000,
                                           node_ids=find_release_nodes(pop_percent=20),
                                                    released_species='gambiae',
                                                    released_genome=[["X", "X"], ["a1", "a1"], ["b1", "b1"]])

    else:
        add_treatment_seeking(camp, targets=[{"trigger": "NewClinicalCase", "coverage": 0.8, "agemin": 15,
                                              "agemax": 70, "seek": 0.4, "rate": 0.3},
                                             {"trigger": "NewSevereCase", "coverage": 0.8, "seek": 0.6, "rate": 0.5}],
                              drug=['Artemether', 'Lumefantrine'],
                              start_day=sim_duration - 10 * 365,
                              broadcast_event_name='Received_Treatment'
                              )
    return camp


def rvg_config_builder(params):
    # params.Combine_Similar_Genomes = True
    params.Include_Vector_State_Columns = False
    params.Include_Death_By_State_Columns = False
    params.Species = 'gambiae'
    params.Stratify_By = 'ALLELE_FREQ'

    """
    E.g.,
    [
        {
            "Allele_Combination": [
                ["X", "X"],
                ["a1", "*"]
            ]
        },
        {
            "Allele_Combination": [
                ["X", "X"],
                ["a0", "a0"]
            ]
        }
    ]
    """
    return params


def rvg_config_builder_genome(params):
    params.Combine_Similar_Genomes = True
    params.Include_Vector_State_Columns = False
    params.Include_Death_By_State_Columns = False
    params.Species = 'gambiae'
    params.Stratify_By = 'SPECIFIC_GENOME'
    params.Specific_Genome_Combinations_For_Stratification = [
        {
            "Allele_Combination": [
                ["X", "X"],
                ["a0", "a0"],
                ["b0", "b0"]
            ]
        },
        {
            "Allele_Combination": [
                ["X", "X"],
                ["a0", "a0"],
                ["b1", "b1"]
            ]
        },
        {
            "Allele_Combination": [
                ["X", "X"],
                ["a1", "a1"],
                ["b1", "b1"]
            ]
        },
        {
            "Allele_Combination": [
                ["X", "X"],
                ["a0", "a1"],
                ["b1", "b1"]
            ]
        },
    ]

    return params


def rvmg_config_builder(params):
    params.Start_Day = 0
    params.End_Day = 10

    return params


def find_release_nodes(pop_percent=20):
    file = 'demographics.json'

    with open(file) as f:
        data = json.load(f)

    nodes = []
    populations = []

    for n in data['Nodes']:
        nodes += [n['NodeID']]
        populations += [n['NodeAttributes']['InitialPopulation']]

    populations = np.array(populations)
    populations = populations / np.sum(populations)

    df = pd.DataFrame({'Nodes': nodes, 'Pop': populations * 100})
    df.sort_values(['Pop'], inplace=True, ascending=False)
    df['Cumulative_Pop'] = df['Pop'].cumsum()
    df = df[df['Cumulative_Pop']<pop_percent]

    return list(df['Nodes'])
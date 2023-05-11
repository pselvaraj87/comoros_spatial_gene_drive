# import pandas as pd
# # import seaborn as sns
# # import matplotlib.pyplot as plt
# import emodpy_malaria.demographics.MalariaDemographics as Demographics
# import emod_api.migration as migration
# from functools import \
#     partial
#
# pop_file = 'data/comoros_population_com_2019-07-01.csv'
# lon_bounds = [43.21885, 43.522264]
# lat_bounds = [-11.942676, -11.362314]
#
# df = pd.read_csv(pop_file)
# df = df[['latitude', 'longitude', 'pop']]
# df = df.astype(float)
#
# df = df[(min(lat_bounds) < df['latitude']) & (df['latitude'] < max(lat_bounds))
#         & (min(lon_bounds) < df['longitude']) & (df['longitude'] < max(lon_bounds))]
#
# df.to_csv(pop_file)
#
#
# pop_file = 'data/comoros_population_com_2019-07-01.csv'
# df = pd.read_csv(pop_file)
# df['pop'] = df['pop']/10
# df.to_csv('data/comoros_population_com_temp.csv')
#
# # fig, ax = plt.subplots(figsize=(6, 12))
# # sns.scatterplot(data=df, x="longitude", y="latitude", ax=ax)
# # plt.title('Comoros population')
# # plt.savefig('figures/Comoros_population.png')
# # plt.show()
#
#
# # input_file = pop_file
# # demographics = Demographics.from_pop_csv(input_file, site='burkina')
# #
# # migration_partial = partial(migration.from_demog_and_param_gravity,
# #                             gravity_params=[7.50395776e-06, 9.65648371e-01, 9.65648371e-01, -1.10305489e+00],
# #                             id_ref='burkina', migration_type=migration.Migration.REGIONAL)
# #
# # demographics
# # migration_partial
import pandas as pd

from analyzers.SpatialOutput import SpatialOutput
import pandas as pd

filename = '/Users/pselvaraj/Downloads/SpatialReport_Population.bin'
with open(filename, mode='rb') as file:  # b is important -> binary
    fileContent = file.read()

so = SpatialOutput.from_bytes(fileContent)
x = so.to_dict()
df = pd.DataFrame(x['data'], columns=x['nodeids']).reset_index()
df.rename(columns={'index': 'time'}, inplace=True)
df = pd.melt(df,
             id_vars='time',
            value_vars=x['nodeids'], # list of days of the week
            var_name='nodes',
            value_name='prevalence')


data = {
    "id": "5740b197-6053-ed11-a9ff-b88303911bc1",
    "name": "comoros_gene_drive",
    "state": "Succeeded",
    "tags": {
        "Baseline": "true",
        "Infected_Progress": "0.1",
        "Larval_Capacity": "8.25",
        "Mortality": "1.1",
        "Run_Number": "3",
        "Serialization": "0",
        "task_type": "emodpy.emod_task.EMODTask",
        "Transmission_To_Human": "0.0"
    },
    "hpc_jobs": [
        {
            "job_id": 9372587,
            "job_state": "Finished",
            "priority": "AboveNormal",
            "working_directory": "/mnt/idm2/home/pselvaraj/output/comoros_gene_drive_20221024_055625/969/b69/a96/969b69a9-6053-ed11-a9ff-b88303911bc1",
            "output_directory_size": 180562599,
            "submit_time": "2022-10-24 05:56:55+00:00",
            "start_time": "2022-10-24 05:56:55+00:00",
            "end_time": "2022-10-24 05:58:07+00:00",
            "configuration": {
                "environment_name": "Calculon",
                "simulation_input_args": "exec Assets/EMOD_ENV_almalinux8.sif Assets/Eradication --config my_config.json --dll-path ./Assets --input-path ./Assets\\;.",
                "working_directory_root": "$COMPS_PATH(USER)/output/comoros_gene_drive_20221024_055625",
                "executable_path": "singularity",
                "node_group_name": "idm_abcd",
                "maximum_number_of_retries": 0,
                "priority": "AboveNormal",
                "min_cores": 48,
                "max_cores": 48,
                "exclusive": False,
                "asset_collection_id": "da8cfc64-1b53-ed11-a9ff-b88303911bc1"
            }
        }
    ]
}
import emodpy_malaria.demographics.MalariaDemographics as Demographics
from emod_api.demographics import DemographicsTemplates as DT
from helpers_demographics import equilibrium_age_distribution
import pandas as pd
import numpy as np

input_file = '/Users/pselvaraj/Github/testing/emodpy-vector_genetics/comoros_exploration/data/comoros_population_com_2019-07-01.csv'
dfnode = pd.read_csv(input_file)
# if 'birthrate' not in dfnode.columns:
# birthrate_file = '/Users/pselvaraj/Github/testing/emodpy-vector_genetics/comoros_exploration/data/WB_birth_rates.csv'
# birthrate_df = pd.read_csv(birthrate_file)
# birthrate = float(birthrate_df[birthrate_df['Country Name'] == 'Comoros']['2020 [YR2020]'].values[0])
# dfnode['birthrate'] = birthrate/np.array(dfnode['pop'])
# dfnode['dist_val'] = dfnode['birthrate'].apply(lambda x: equilibrium_age_distribution(x, 2.74e-06, birthrate)[1])
# dfnode['result_val'] = dfnode['birthrate'].apply(lambda x: equilibrium_age_distribution(x, 2.74e-06, birthrate)[0])
# dfnode.to_csv(input_file)

demographics = Demographics.from_pop_csv(input_file, site='comoros')
DT.MortalityRateByAge(demographics, age_bins=[0], mort_rates=[23])
DT.InitPrevUniform(demographics, low_prev=0.13, high_prev=0.15)
DT.FullMigrationHeterogeneity(demographics)





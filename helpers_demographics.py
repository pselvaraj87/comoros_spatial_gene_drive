import numpy as np
import scipy.integrate as sp


def equilibrium_age_distribution(birth_rate, mort_scale, mort_value, max_age=100):
    """
    This function generates an initial age distribution already at equilibrium,
    allowing you to run burn-ins for immunity establishment only.

    NB: You must set your config file to Birth_Rate_Dependence="FIXED_BIRTH_RATE" and
                                Age_Initialization_Distribution_Type= "DISTRIBUTION_COMPLEX"
    for this to work. If you have age-dependent birth rates, go talk to Kurt for a modified script.

    :param birth_rate: population birth rate, in units of births/node/day
    :param mort_scale: demographics["Defaults"]["IndividualAttributes"]["MortalityDistribution"]["ResultScaleFactor"]
    :param mort_value: annual deaths per 1000, set to mirror birth rates.
    :param max_age: age past which you want everyone to die. In the current implementation you will get *some* mass
                in this bin, but only the amound specified in the equilibrium distribution.
               If you have more age-specific mortality rates, you can implement them with only a slight modification to this script.
    :return: resval (a list of age bins in days from 0 to max_yr),
                and distval (a list of the cumulative proportion of the population in each age bin of resval)
    """

    # define daily mortality probability by age group
    age_vec = [365 * 0, 365 * (max_age - 0.001), 365 * max_age]
    mort_vec = [mort_scale * mort_value, mort_scale * mort_value, 1]

    max_yr = 120
    day_to_year = 365

    # add bounds around age_vec and mort_vec
    mvec_x = [-1] + age_vec + [max_yr * day_to_year + 1]
    mvec_y = [mort_vec[0]] + mort_vec + [mort_vec[-1]]

    # cumulative monthly survival probabilities
    m_x = np.arange(0, max_yr * day_to_year, 30)
    mval = (1.0 - np.interp(m_x, xp=mvec_x, fp=mvec_y)) ** 30

    # create normalized population pyramid
    popvec = birth_rate * np.cumprod(mval)
    tpop = np.trapz(popvec, x=m_x)
    npvec = popvec / tpop
    # what proportion of  people are in the simulation up to bin i
    cpvec = sp.cumtrapz(npvec, x=m_x, initial=0)

    resval = np.around(np.linspace(0, max_yr * day_to_year))
    distval = np.interp(resval, xp=m_x, fp=cpvec)

    return resval, distval


# def update_nodes_birth_mortality(demog, birth_rate):
#
#     nodes = []
#     for i, node in enumerate(demog.nodes):
#         # if res_in_degrees is custom assume node_ids are generated for a household-like setup and not based on lat/lon
#         if node.forced_id:
#             node_id = node.forced_id
#         node_attributes = node.to_dict()
#         individual_attributes = {}
#
#         node_attributes.update({'LarvalHabitatMultiplier': 1.0})
#
#
#         # if nodes are in different countries, find node-specific mortality rates
#         if "Country" in node_attributes.keys():
#
#             mod_mortality = {
#                 "NumDistributionAxes": 2,
#                 "AxisNames": ["gender", "age"],
#                 "AxisUnits": ["male=0,female=1", "years"],
#                 "AxisScaleFactors": [1, 365],
#                 "NumPopulationGroups": [2, 1],
#                 "PopulationGroups": [
#                     [0, 1],
#                     [0]
#                 ],
#                 "ResultUnits": "annual deaths per 1000 individuals",
#                 "ResultScaleFactor": mort_scale_factor,
#                 "ResultValues": [
#                     [birth_rate],
#                     [birth_rate]
#                 ]
#             }
#
#             node_attributes.pop("Country", None)
#
#             individual_attributes.update({"MortalityDistribution": mod_mortality})
#
#
#         # equilibrium age distribution
#         per_node_birth_rate = (float(node.pop) / 1000) * birth_rate / 365.0
#         node_attributes.update({'BirthRate': per_node_birth_rate})
#         resval, distval = equilibrium_age_distribution(per_node_birth_rate, mort_scale_factor,
#                                                             birth_rate)
#         mod_age = {
#             "DistributionValues":
#                 [
#                     distval.tolist()
#                 ],
#             "ResultScaleFactor": 1,
#             "ResultValues":
#                 [
#                     resval.tolist()
#                 ]
#         }
#         individual_attributes.update({"AgeDistribution": mod_age})
#
#         nodes.append({'NodeID': node_id,
#                       'NodeAttributes': node_attributes,
#                       'IndividualAttributes': individual_attributes})
#
#     return nodes
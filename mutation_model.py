from scipy import stats
import numba as nba
import numpy as np
import pandas as pd
import time
import os

BP_MUTATION_RATE = 2.91*(10**-9) # from Table 1 in Ben's paper https://www.nature.com/articles/s41467-020-14844-6/tables/1
BP_IN_EXOME = (45-18)*(10**6) # 45Mb - 18Mb of potential synonymous mutations http://www.nature.com/nature/journal/v536/n7616/full/nature19057.html?foxtrotcallback=true
CLONE_MUTATION_RATE = 1 - stats.binom.pmf(k=0, n=BP_IN_EXOME, p=BP_MUTATION_RATE) ## Or 0.1130796 from Reviewer 1, based on Ben's paper

## Probabilities of where mutations might land
GENES_IN_GENOME = 20412  # https://en.wikipedia.org/wiki/Human_genome but source is http://useast.ensembl.org/Homo_sapiens/
DRIVER_GENES = 20
BLOCKADE_GENES = 1
IMMUNOSUPPRESSIVE_GENES = 1
DELETERIOUS_GENES = GENES_IN_GENOME - DRIVER_GENES - BLOCKADE_GENES - IMMUNOSUPPRESSIVE_GENES

DELETERIOUS_MUTATION_PROB = DELETERIOUS_GENES/GENES_IN_GENOME  # Probability mutation only increases antigenicity
BLOCKADE_MUTATION_PROB = BLOCKADE_GENES / GENES_IN_GENOME
IMMUNOSUPPRESSIVE_MUTATION_PROB = IMMUNOSUPPRESSIVE_GENES/GENES_IN_GENOME
INITIAL_DRIVER_GENE_MUTATION_PROB = DRIVER_GENES/GENES_IN_GENOME  # probability decreases as unique number of driver genes decreases
assert np.isclose(DELETERIOUS_MUTATION_PROB + BLOCKADE_MUTATION_PROB + IMMUNOSUPPRESSIVE_MUTATION_PROB + INITIAL_DRIVER_GENE_MUTATION_PROB, 1)

# Related to mutation #
DRIVER_MUT_IDX = 0
BLOCKADE_MUT_IDX = 1
SUPPRESSION_MUT_IDX = 2
DELETERIOUS_MUT_IDX = 3  # Last because any left over probability in multinomial goes to last bin

# Related to phenotype #
DRIVERS_FOR_A = 2
DRIVERS_FOR_C = 4
BLOCKADE_PROTECTION = 1.0
IMMUNOSUPPRESSION = 0.2
N_ANTIGENS = 10000
E_IDX = 0
A_IDX = 1
C_IDX = 2

# Related to phenotype lookup. Values will be set by init_results_mat() #
N_BLOCKADE_STRATS = None
N_SUPPRESSION_STRATS = None
N_STRATS = None


# Related to columns in parameter matrix #

BLOCKADE_COL_IDX = 0
SUPPRESSION_COL_IDX = 1
ANT_COL_VAL_IDX = 2
DRIVER_COL_IDX = 3
SPECIES_COL_IDX = 4
ANT_ID_COL_IDX = 5
CLONE_ID_COL_IDX = 6
PARENT_COL_IDX = 7
T_CELL_INFILTRATION_IDX = 8
DIV_R_IDX = 9
IMMUNE_KILL_IDX = 10
NET_R_IDX = 11
RESULTS_MAT_ORIGIN_IDX = 12
##NOTE Reminder, if any idx above change, update INFO_COLS
INFO_COLS = ["blockade", "suppression", "antigenicity", "n_drivers", "species", "antigen_idx", "clone", "parent", "t_cell_infiltration", "div_rate", "immune_kill", "net_r", "origin"]


N_INFO_COLS = len(INFO_COLS)


PROLIF_RATIO = 1.052267  # Based on the increase in Ki67 density from CRC to CRA. See stats_ppc.R for calculation
E_DIV_RATE = 0.2  # discussed in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4007491/
A_DIV_RATE = E_DIV_RATE * PROLIF_RATIO
C_DIV_RATE = A_DIV_RATE * PROLIF_RATIO

DIV_RATES = np.array([E_DIV_RATE, A_DIV_RATE, C_DIV_RATE], dtype=np.float32)
E_K = 10 ** 7
A_K = 10 ** 8
C_K = 10 ** 9

CARRAYING_CAPACITIES = np.array([E_K, A_K, C_K], dtype=np.int64)
DEATH_RATE = 0.01
COMP_MATRIX = np.ones((3, 3), dtype=np.uint8)  # EAC vs EAC. E and A do not compete
COMP_MATRIX[0, 1] = 0
COMP_MATRIX[1, 0] = 0


# Related to ouput #
NREPS = 365*100
RECORD_FREQ = 30
RECORD_TIME = False  # Whether or not to record the simulation over time
FOUNDER_THRESH = 100  # population must reach this size for a tumor to be considered established
NO_TUMOR_STATUS = 0
EXTINCT_STATUS = 1
CONTROLLED_STATUS = 2
EXTANT_STATUS = 3
STEPS_AS_CARCINOMA_STOP_THRESH = 365
MAX_SIMS = 10
INITIAL_E_SIZE = int(E_K*(1 - DEATH_RATE))

# Info to record #
SUMMARY_SPECIES_COL_IDX = 0
ORIGIN_COL_IDX = 1
RECORD_TIME_COL_IDX = 2
W_MEAN_ANT_COL_IDX = 3
MEAN_S_STRAT_COL_IDX = 4
MEAN_B_STRAT_COL_IDX = 5
MEAN_SIZE_IDX = 6
TOTAL_SIZE_IDX = 7
TOTAL_SUPPRESSION_IDX = 8
TOTAL_BLOCKADE_IDX = 9
TOTAL_T_IDX = 10
FOUNDER_ANT_IDX = 11
FOUNDER_S_IDX = 12
FOUNDER_B_IDX = 13
N_PHENO_IDX = 14
MEAN_ANT_COL_IDX = 15

from collections import OrderedDict
summary_colname_dict = OrderedDict([(SUMMARY_SPECIES_COL_IDX, "species"),
                                    (ORIGIN_COL_IDX, "origin"),
                                    (RECORD_TIME_COL_IDX, "Time"),
                                    (W_MEAN_ANT_COL_IDX, "w_mean_ant"),
                                    (MEAN_S_STRAT_COL_IDX, "w_mean_suppression"),
                                    (MEAN_B_STRAT_COL_IDX, "w_mean_blockade"),
                                    (MEAN_SIZE_IDX, "mean_size"),
                                    (TOTAL_SIZE_IDX, "Total_Size"),
                                    (TOTAL_SUPPRESSION_IDX, "Total_Suppressive"),
                                    (TOTAL_BLOCKADE_IDX, "Total_Blockade"),
                                    (TOTAL_T_IDX, "Total_Tcell"),
                                    (FOUNDER_ANT_IDX, "founder_ant"),
                                    (FOUNDER_S_IDX, "founder_suppression"),
                                    (FOUNDER_B_IDX, "founder_protection"),
                                    (N_PHENO_IDX, "n_pheno"),
                                    (MEAN_ANT_COL_IDX, "mean_ant")
                                    ])


summary_info_cols = list(summary_colname_dict.values())

@nba.njit()
def intersect1d(ar1, ar2):

    ar1 = np.unique(ar1)
    ar2 = np.unique(ar2)

    aux = np.concatenate((ar1, ar2))

    aux.sort()

    mask = aux[1:] == aux[:-1]
    int1d = aux[:-1][mask]

    return int1d

@nba.njit()
def weighted_mean(x, w):

    unique_w = np.unique(w)
    if len(unique_w) == 1:
        return np.mean(x)
    else:
        return np.sum((x*w))/np.sum(w)

@nba.njit()
def add_founder_info(summary_mat, results_mat, founder_id, time_added):
    founder_species = int(results_mat[founder_id, SPECIES_COL_IDX])
    summary_mat[founder_species, SUMMARY_SPECIES_COL_IDX] = founder_species
    summary_mat[founder_species, FOUNDER_ANT_IDX] = results_mat[founder_id, ANT_COL_VAL_IDX]
    summary_mat[founder_species, FOUNDER_B_IDX] = results_mat[founder_id, BLOCKADE_COL_IDX]
    summary_mat[founder_species, FOUNDER_S_IDX] = results_mat[founder_id, SUPPRESSION_COL_IDX]
    summary_mat[founder_species, ORIGIN_COL_IDX] = time_added


@nba.njit()
def summarize_time_pt(summary_mat, results_mat, current_sizes, species_id, time_pt):


    species_idx = np.where(results_mat[:, SPECIES_COL_IDX] == species_id)[0]
    record_idx = intersect1d(species_idx, np.where(current_sizes > 0)[0])
    if len(record_idx) > 1 and species_id == E_IDX:
        record_idx = record_idx[record_idx != 0] ### Don't use root in calculation


    summary_mat[species_id, SUMMARY_SPECIES_COL_IDX] = species_id
    summary_mat[species_id, RECORD_TIME_COL_IDX] = time_pt
    if len(record_idx) == 0:
        ### Whole population extinct. Should only happen in Mueller's ratchet, where E goes extinct by mutating too much
        summary_mat[species_id, W_MEAN_ANT_COL_IDX] = np.nan
        summary_mat[species_id, MEAN_ANT_COL_IDX] = np.nan
        summary_mat[species_id, MEAN_S_STRAT_COL_IDX] = np.nan
        summary_mat[species_id, MEAN_B_STRAT_COL_IDX] = np.nan
        summary_mat[species_id, MEAN_SIZE_IDX] = np.nan
        summary_mat[species_id, TOTAL_SIZE_IDX] = 0
        summary_mat[species_id, TOTAL_SUPPRESSION_IDX] = 0
        summary_mat[species_id, TOTAL_BLOCKADE_IDX] = 0
        summary_mat[species_id, TOTAL_T_IDX] = 0
        summary_mat[species_id, N_PHENO_IDX] = 0

    else:
        species_sizes = current_sizes[record_idx]
        summary_mat[species_id, W_MEAN_ANT_COL_IDX] = weighted_mean(results_mat[record_idx, ANT_COL_VAL_IDX], species_sizes)
        summary_mat[species_id, MEAN_ANT_COL_IDX] = np.mean(results_mat[record_idx, ANT_COL_VAL_IDX])
        summary_mat[species_id, MEAN_S_STRAT_COL_IDX] = weighted_mean(results_mat[record_idx, SUPPRESSION_COL_IDX], species_sizes)
        summary_mat[species_id, MEAN_B_STRAT_COL_IDX] = weighted_mean(results_mat[record_idx, BLOCKADE_COL_IDX], species_sizes)
        summary_mat[species_id, MEAN_SIZE_IDX] = np.mean(species_sizes)
        summary_mat[species_id, TOTAL_SIZE_IDX] = np.sum(species_sizes)
        summary_mat[species_id, TOTAL_SUPPRESSION_IDX] = np.sum(species_sizes*results_mat[record_idx, SUPPRESSION_COL_IDX])
        summary_mat[species_id, TOTAL_BLOCKADE_IDX] = np.sum(species_sizes*results_mat[record_idx, BLOCKADE_COL_IDX])
        summary_mat[species_id, TOTAL_T_IDX] = np.sum(species_sizes*results_mat[record_idx, T_CELL_INFILTRATION_IDX])
        summary_mat[species_id, N_PHENO_IDX] = len(species_sizes)


def init_results_mat():
    '''
    No support for np.meshgrid. Cannot be jitted
    :return:
    '''

    global N_BLOCKADE_STRATS
    global N_SUPPRESSION_STRATS
    global N_STRATS

    protection_vals = np.unique(np.array([0, BLOCKADE_PROTECTION]))
    suppression_vals = np.unique(np.array([0, IMMUNOSUPPRESSION]))
    N_BLOCKADE_STRATS = len(protection_vals)
    N_SUPPRESSION_STRATS = len(suppression_vals)
    N_STRATS = N_BLOCKADE_STRATS*N_SUPPRESSION_STRATS

    ant_array = np.linspace(0, 1, N_ANTIGENS)
    driver_vals = np.arange(0, DRIVERS_FOR_C + 1)


    param_array = np.array(np.meshgrid(protection_vals, suppression_vals, ant_array, driver_vals)).T.reshape(-1, 4)

    if RECORD_TIME:
        results_mat = np.zeros((param_array.shape[0], N_INFO_COLS + int(np.ceil(NREPS / RECORD_FREQ))))
    else:
        results_mat = np.zeros((param_array.shape[0], N_INFO_COLS))

    results_mat[:, 0:4] = param_array

    results_mat[:, SPECIES_COL_IDX] = get_species_idx(results_mat[:, DRIVER_COL_IDX])
    results_mat[:, ANT_ID_COL_IDX] = get_ant_idx(np.arange(results_mat.shape[0]), results_mat[:, BLOCKADE_COL_IDX], results_mat[:, SUPPRESSION_COL_IDX], results_mat[:, DRIVER_COL_IDX], N_ANTIGENS, N_STRATS, N_SUPPRESSION_STRATS)
    results_mat[:, CLONE_ID_COL_IDX] = np.arange(0, results_mat.shape[0])
    results_mat[:, PARENT_COL_IDX] = -1

    non_zero_ant_idx = np.where(results_mat[:, ANT_COL_VAL_IDX] > 0)
    results_mat[non_zero_ant_idx, T_CELL_INFILTRATION_IDX] = 1.0 - (results_mat[non_zero_ant_idx, SUPPRESSION_COL_IDX] / results_mat[non_zero_ant_idx, ANT_COL_VAL_IDX])
    results_mat[:, T_CELL_INFILTRATION_IDX][results_mat[:, T_CELL_INFILTRATION_IDX] < 0] = 0
    results_mat[:, IMMUNE_KILL_IDX] = results_mat[:, T_CELL_INFILTRATION_IDX] * (1 - results_mat[:, BLOCKADE_COL_IDX]) * results_mat[:, ANT_COL_VAL_IDX]
    results_mat[:, DIV_R_IDX] = DIV_RATES[results_mat[:, SPECIES_COL_IDX].astype(np.int)] #[DIV_RATES[i] for i in results_mat[:, SPECIES_COL_IDX]]
    results_mat[:, NET_R_IDX] = results_mat[:, DIV_R_IDX] - results_mat[:, IMMUNE_KILL_IDX] - DEATH_RATE
    results_mat[:, RESULTS_MAT_ORIGIN_IDX] = -1

    summary_mat = np.full((3, len(summary_info_cols)), -1.0)
    add_founder_info(summary_mat, results_mat, 0, 0)

    return results_mat, summary_mat


@nba.vectorize()
def get_species_idx(n_drivers):
    if n_drivers < DRIVERS_FOR_A:
        species = E_IDX
    elif n_drivers < DRIVERS_FOR_C:
        species = A_IDX
    else:
        species = C_IDX

    return species

@nba.vectorize()
def get_pheno_idx(blockade, suppression, ant_idx, n_drivers, n_antigens, n_strats, n_suppression):
    if blockade > 0:
        blockade_strat_idx = 1
    else:
        blockade_strat_idx = 0

    if suppression > 0:
        suppression_strat_idx = 1
    else:
        suppression_strat_idx = 0

    driver_idx = n_strats * n_antigens * n_drivers
    mat_ant_idx = ant_idx * n_strats
    blockade_idx = n_suppression*blockade_strat_idx
    pheno_id = driver_idx + mat_ant_idx + blockade_idx + suppression_strat_idx

    return int(pheno_id)


@nba.vectorize()
def get_ant_idx(pheno_idx, blockade, suppression, n_drivers, n_antigens, n_strats, n_suppression):
    if blockade > 0:
        blockade_strat_idx = 1
    else:
        blockade_strat_idx = 0

    if suppression > 0:
        suppression_strat_idx = 1
    else:
        suppression_strat_idx = 0

    driver_idx = n_strats * n_antigens * n_drivers
    blockade_idx = n_suppression * blockade_strat_idx
    ant_idx = (pheno_idx - driver_idx - blockade_idx - suppression_strat_idx)/n_strats

    return int(ant_idx)



@nba.njit()
def get_phenotype_info(results_mat, pheno_id):
    blockade = results_mat[pheno_id, BLOCKADE_COL_IDX]
    suppression = results_mat[pheno_id, SUPPRESSION_COL_IDX]
    ant = results_mat[pheno_id, ANT_COL_VAL_IDX]
    n_drivers = results_mat[pheno_id, DRIVER_COL_IDX]
    species = int(results_mat[pheno_id, SPECIES_COL_IDX])
    ant_idx = int(results_mat[pheno_id, ANT_ID_COL_IDX])

    return blockade, suppression, ant, n_drivers, species, ant_idx


def check_lookup():
    results_mat, summary_mat = init_results_mat()
    for i in range(results_mat.shape[0]):
        blockade, suppression, ant, n_drivers, species, ant_idx = get_phenotype_info(results_mat, i)

        pid = get_pheno_idx(blockade, suppression, ant_idx, n_drivers, N_ANTIGENS, N_STRATS, N_SUPPRESSION_STRATS)
        if i != pid:
            raise ValueError("indices do not match. Row is", i, "But phenotype defined as", pid)

        ant_idx2 = get_ant_idx(i, blockade, suppression, n_drivers, N_ANTIGENS, N_STRATS, N_SUPPRESSION_STRATS)
        if ant_idx != ant_idx2:
            raise ValueError("Antigen indices don't match. Value should be", ant_idx, "but is determined to be ", ant_idx2)

check_lookup()
print("lookup checks out")

@nba.njit()
def calc_d(total_antigenicity, suppression, blockade_protection):
    '''
    protection reduces number killed, while suppression reduces rate at which killing occurs
    :param total_antigenicity: antigenicity
    :param suppression: immune suppression
    :param blockade_protection: protection from blockade
    :return:
    '''

    if total_antigenicity == 0:
        return 0

    t_cell_infiltration_rate = 1 - (suppression/total_antigenicity)
    if t_cell_infiltration_rate < 0:
        t_cell_infiltration_rate = 0

    d = t_cell_infiltration_rate * (1 - blockade_protection) * total_antigenicity

    if d < 0:
        d = 0

    return d

@nba.njit()
def calc_r(size, species, species_sizes, suppression, with_comp=True):
    '''
    Calculate growth rate with competition and mutualism. Also used in Tissue.test_phenotype()
    Based on equations from

    :param species: index of species
    :param species_sizes: array with sizes of all species
    :param suppression: immune suppression parameter. Can increase growth rate
    :return:
    '''

    K = CARRAYING_CAPACITIES[species]
    r = DIV_RATES[species]
    if with_comp:
        total_n = np.sum(COMP_MATRIX[species, ] * species_sizes) ### Clone is competing with all other clones, and experiences same competition as others in the species
    else:
        total_n = species_sizes[species]

    total_b = suppression * size ### benefit is only for this clone

    return r*((K - total_n + 0.5*total_b)/K) ###NOTE multiply by 0.5 so that maximum size is 2x K


@nba.njit()
def get_mutant_phenotype(mutation_idx, n_drivers, blockade, suppression):
    '''
    :return: phenotype values of new mutant.
    '''

    if mutation_idx == DELETERIOUS_MUT_IDX:
        return n_drivers, blockade, suppression

    elif mutation_idx == DRIVER_MUT_IDX:
        new_drivers = n_drivers + 1
        if new_drivers > DRIVERS_FOR_C:
            new_drivers = DRIVERS_FOR_C

        return new_drivers, blockade, suppression

    elif mutation_idx == BLOCKADE_MUT_IDX:
        return n_drivers, BLOCKADE_PROTECTION, suppression

    elif mutation_idx == SUPPRESSION_MUT_IDX:
        return n_drivers, blockade, IMMUNOSUPPRESSION


@nba.njit()
def get_mutant_types(n_mutants, n_drivers):
    '''
    :param n_mutants: how many mutants to create
    :param n_drivers: number of driver genes aleary mutated

    :return: number of each mutant type to create.
    '''
    new_driver_prob = (DRIVER_GENES - n_drivers)/GENES_IN_GENOME

    mutant_probs = np.array([new_driver_prob,
                             BLOCKADE_MUTATION_PROB,
                             IMMUNOSUPPRESSIVE_MUTATION_PROB,
                             DELETERIOUS_MUTATION_PROB + n_drivers / GENES_IN_GENOME]) ### driver gene already hit, and still gets antigenicity

    return np.random.multinomial(n_mutants, mutant_probs)

@nba.njit()
def add_mutants(current_sizes, next_sizes, next_species_sizes, pheno_id, n_mutants, results_mat, established_adenoma, established_carcinoma, time_pt):
    blockade, suppression, ant, n_drivers, species, ant_idx = get_phenotype_info(results_mat, pheno_id)
    new_mutant_types = get_mutant_types(n_mutants, n_drivers)
    # print("Adding mutants:", new_mutant_types)
    for mut_idx in range(len(new_mutant_types)):
        n_mut_of_type = new_mutant_types[mut_idx]
        # print("Adding", n_mut_of_type, "new mutants")
        if n_mut_of_type == 0:
            continue

        new_cln_n_drivers, new_cln_blockade, new_cln_suppression = get_mutant_phenotype(mut_idx, n_drivers, blockade, suppression)

        if species == E_IDX and DRIVERS_FOR_A <= new_cln_n_drivers < DRIVERS_FOR_C and established_adenoma:
            # Put new mutant back in parental epithelial population because can't add new adenoma. Without this, population would "disappear"
            next_sizes[pheno_id] += 1
            next_species_sizes[int(results_mat[pheno_id, SPECIES_COL_IDX])] += 1
            continue

        elif species == A_IDX and new_cln_n_drivers >= DRIVERS_FOR_C and established_carcinoma:
            # Put new mutant back in parental adenoma population because can't add new carcinoma. Without this, population would "disappear"
            next_sizes[pheno_id] += 1
            next_species_sizes[int(results_mat[pheno_id, SPECIES_COL_IDX])] += 1
            continue

        # Using a loop avoids creating massive array that takes up large amounts of memory.
        # If jitted, performance should not be too affected
        for mut in range(n_mut_of_type):
            min_next_ant = ant_idx + 1
            if min_next_ant >= N_ANTIGENS:
                new_ant_idx = N_ANTIGENS - 1
            else:
                new_ant_idx = np.random.randint(min_next_ant, N_ANTIGENS)

            new_pheno_idx = get_pheno_idx(new_cln_blockade, new_cln_suppression, new_ant_idx, new_cln_n_drivers, N_ANTIGENS, N_STRATS, N_SUPPRESSION_STRATS)
            if results_mat[new_pheno_idx, NET_R_IDX] < 0:
                continue

            if results_mat[new_pheno_idx, ANT_COL_VAL_IDX] > 1.0:
                print("phenotype id = ", new_pheno_idx, "new antigenicity idx =", new_ant_idx, " antigenicity=", results_mat[new_pheno_idx, ANT_COL_VAL_IDX])
                raise Warning("Antigenicity out of range")


            if results_mat[new_pheno_idx, SUPPRESSION_COL_IDX] > 1.0:
                print("phenotype id = ", new_pheno_idx, "suppression =", results_mat[new_pheno_idx, SUPPRESSION_COL_IDX])
                raise Warning("Suppression out of range")

            if results_mat[new_pheno_idx, BLOCKADE_COL_IDX] > 1.0:
                print("phenotype id = ", new_pheno_idx, "blockade =", results_mat[new_pheno_idx, BLOCKADE_COL_IDX])
                raise Warning("Blockade out of range")


            if results_mat[new_pheno_idx, DRIVER_COL_IDX] > DRIVERS_FOR_C:
                print("phenotype id = ", new_pheno_idx, "blockade =", results_mat[new_pheno_idx, DRIVERS_FOR_C])
                raise Warning("Number of drivers out of range")

            ### Verify that child ant >= parent ant
            if results_mat[new_pheno_idx, ANT_COL_VAL_IDX] < results_mat[pheno_id, ANT_COL_VAL_IDX]:
                print("child antigenicity idx =", new_ant_idx, "child antigenicity=", results_mat[new_pheno_idx, ANT_COL_VAL_IDX], "\n",
                              "parent antigenicity idx =", ant_idx, "parent antigenicity=", results_mat[pheno_id, ANT_COL_VAL_IDX]
                      )
                raise Warning("child less antigenic than parent")


            if results_mat[new_pheno_idx, BLOCKADE_COL_IDX] < results_mat[pheno_id, BLOCKADE_COL_IDX]:
                print("pheno id", new_pheno_idx, "child blockade=", results_mat[new_pheno_idx, BLOCKADE_COL_IDX], "\n",
                              "parent blockade =", results_mat[pheno_id, BLOCKADE_COL_IDX]
                      )
                raise Warning("child lost blockade")

            if results_mat[new_pheno_idx, SUPPRESSION_COL_IDX] < results_mat[pheno_id, SUPPRESSION_COL_IDX]:
                print("pheno id", new_pheno_idx, "child suppression =", results_mat[new_pheno_idx, SUPPRESSION_COL_IDX], "\n",
                              "parent suppression =", results_mat[pheno_id, SUPPRESSION_COL_IDX]
                      )
                raise Warning("child lost suppression")

            if results_mat[new_pheno_idx, DRIVER_COL_IDX] < results_mat[pheno_id, DRIVER_COL_IDX]:
                print("child drivers =", results_mat[new_pheno_idx, DRIVER_COL_IDX], "\n",
                              "parent drivers =", results_mat[pheno_id, DRIVER_COL_IDX]
                      )
                raise Warning("child lost driver mutation")

            if results_mat[new_pheno_idx, DRIVER_COL_IDX] > results_mat[pheno_id, DRIVER_COL_IDX] + 1:
                print("child drivers =", results_mat[new_pheno_idx, DRIVER_COL_IDX], "should be", new_cln_n_drivers,"\n",
                              "parent drivers =", results_mat[pheno_id, DRIVER_COL_IDX]
                      )

            if current_sizes[new_pheno_idx] == 0 and next_sizes[new_pheno_idx] == 0:
                results_mat[new_pheno_idx, PARENT_COL_IDX] = pheno_id
                results_mat[new_pheno_idx, RESULTS_MAT_ORIGIN_IDX] = time_pt + 1

            next_sizes[new_pheno_idx] += 1
            next_species_sizes[int(results_mat[new_pheno_idx, SPECIES_COL_IDX])] += 1


@nba.njit()
def grow(current_sizes, pheno_id, results_mat, species_sizes):
    '''
    Assume that clones can both divide and die within the same timestep
    :param current_sizes: current sizes of each phenotype
    :param pheno_id: row index of phenotype in results_mat
    :param results_mat: matrix containing population sizes for each possible phenotype
    :param species_sizes: current sizes of each "species" (epithelial, adenoma, carcinoma)
    :return:
    '''

    size = current_sizes[pheno_id]
    blockade, suppression, ant, n_drivers, species, ant_idx = get_phenotype_info(results_mat, pheno_id)

    ### Assume that if growth rate is <= 0, the cells are experiencing contact inhibition and can't divide.
    ### Also that invading populations consume resources, taking them away from current species, preventing divison.
    r = calc_r(size, species, species_sizes, suppression)
    loss_rate = results_mat[pheno_id, IMMUNE_KILL_IDX] + DEATH_RATE
    if r > 0:
        new_n = np.random.binomial(size, r)
        n_mutated = np.random.binomial(new_n, CLONE_MUTATION_RATE)
    else:
        new_n = 0
        n_mutated = 0
        loss_rate += abs(r)

    if loss_rate > 1:
        loss_rate = 1

    total_lost = np.random.binomial(size, loss_rate)

    next_size = size + new_n - n_mutated - total_lost
    if next_size < 0:
        next_size = 0

    return next_size, n_mutated


@nba.njit()
def erase_history(keep_idx, results_mat, remove_species, current_sizes):
    """
    Remove other species that have have formed before founder population established.
    Example: several very small adenoma populations exists, and are below `FOUNDER_THRESH`.
    Only want to track the 1st adenoma/carcinoma that forms, so when one of those goes above FOUNDER_THRESH
    will remove the other populations that belong to that species.

    :param keep_idx: row index of phenotype to keep
    :param results_mat: matrix containing population sizes of each phenotype
    :param remove_species: index of species to remove
    :param current_sizes: current sizes of each phenotype

    """
    all_species_idx = np.where(results_mat[:, SPECIES_COL_IDX] == remove_species)[0]
    for idx in all_species_idx:
        if idx != keep_idx:
            current_sizes[idx] = 0
            if RECORD_TIME:
                results_mat[idx, N_INFO_COLS:] = 0

@nba.njit()
def step(iter_order, current_sizes, next_sizes, summary_mat, results_mat, species_sizes, next_species_sizes, time_pt, established_adenoma, established_carcinoma):
    """
    Simulataes one time step in the simulation

    :param iter_order: order in which to simulate each phenotype.
    :param current_sizes: current sizes of each phenotype
    :param next_sizes: sizes of each phenotype in the next step. Initially current sizes, but gets updated in this function
    :param summary_mat: sizes of each species at all timepoints
    :param results_mat: sizes of each phenotype at all timepoints
    :param species_sizes: current size of each "species" (epithelial, adenoma, carcinoma)
    :param next_species_sizes: sizes of each species in the next step. Initially species_sizes, but gets updated in this function
    :param time_pt: current time point
    :param established_adenoma: Has an adenoma's population size gone above FOUNDER_THRESH?
    :param established_carcinoma: Has a carcinoma's population size gone above FOUNDER_THRESH?

    """

    for i in iter_order:
        if current_sizes[i] == 0 or results_mat[i, NET_R_IDX] < 0:
            continue

        cln_new_size, n_mutants = grow(current_sizes, i, results_mat, species_sizes)

        if cln_new_size == 0:
            if n_mutants == 0:
                continue

        next_sizes[i] += cln_new_size
        next_species_sizes[int(results_mat[i, SPECIES_COL_IDX])] += cln_new_size


        if n_mutants > 0:
            ### Don't let adenomas mutate until one has been recorded. Same with carcinoma. Otherwise, determining founder population isn't accurate. Earliest and most fit phenotype should still reach threshold.
            ### Only an issue with very high mutation rates (e.g. 2e-6). Shouldn't affect things too much, as growth is rapid, and so effect will mostly be delay in time
            current_species = int(results_mat[i, SPECIES_COL_IDX])
            if (current_species == A_IDX and not established_adenoma) or (current_species == C_IDX and not established_carcinoma):
                print("Had to block mutation of", i, "species=", current_species, "at time", time_pt, "size=", cln_new_size)
                next_sizes[i] += n_mutants
                next_species_sizes[current_species] += n_mutants
                cln_new_size += n_mutants
            else:
                add_mutants(current_sizes, next_sizes, next_species_sizes, i, n_mutants, results_mat, established_adenoma,
                            established_carcinoma, time_pt)

        if time_pt % RECORD_FREQ == 0 and cln_new_size > 0 and RECORD_TIME:
            results_mat[i, N_INFO_COLS + time_pt//RECORD_FREQ] = cln_new_size

        if cln_new_size >= FOUNDER_THRESH and (not established_adenoma or not established_carcinoma):
            species = results_mat[i, SPECIES_COL_IDX]
            if species == A_IDX and not established_adenoma:
                # print("Establised adenoma at time", time_pt)
                established_adenoma = True
                erase_history(i, results_mat, A_IDX, current_sizes)
                next_species_sizes[A_IDX] = cln_new_size

                add_founder_info(summary_mat, results_mat, i, time_pt + 1)
                print("Established adenoma", i,  "at time", time_pt + 1,  "size=", cln_new_size, "species=", species, "with ant=", summary_mat[int(species), FOUNDER_ANT_IDX], ", S=", summary_mat[int(species), FOUNDER_S_IDX], "B=", summary_mat[int(species), FOUNDER_B_IDX])

            if species == C_IDX and not established_carcinoma:
                established_carcinoma = True
                erase_history(i, results_mat, C_IDX, current_sizes)
                next_species_sizes[C_IDX] = cln_new_size

                ### Record what populations looked like immediately before carcinoma
                add_founder_info(summary_mat, results_mat, i, time_pt + 1)
                print("Established carinoma", i, "at time", time_pt + 1,  "size=", cln_new_size, "species=", species, "with ant=", summary_mat[int(species), FOUNDER_ANT_IDX], ", S=", summary_mat[int(species), FOUNDER_S_IDX], "B=", summary_mat[int(species), FOUNDER_B_IDX])
                summarize_time_pt(summary_mat, results_mat, current_sizes, E_IDX, time_pt)
                summarize_time_pt(summary_mat, results_mat, current_sizes, A_IDX, time_pt)

                if not established_adenoma:
                    ### Could still happen if no adenoma got above record thresh, but one created carcinoma anyway, should be rare unless mutation rate is very high
                    print("Established carcinoma, but didn't record adenoma. Current adenoma size is", species_sizes[A_IDX])
                    established_adenoma = True

                    current_p = int(results_mat[i, PARENT_COL_IDX])
                    current_p_drivers = results_mat[current_p, DRIVER_COL_IDX]
                    while current_p_drivers != DRIVERS_FOR_A:
                        # print("searching for parental CRA. current parent is", current_p, "of species", current_p_species)
                        next_p = int(results_mat[current_p, PARENT_COL_IDX])
                        next_p_drivers = results_mat[next_p, DRIVER_COL_IDX]
                        print("searching for parental CRA. current parent is", current_p, "with n drivers = ", next_p_drivers, "next p=", next_p, "with n divers=", next_p_drivers)
                        current_p = next_p
                        current_p_drivers = next_p_drivers

                        if current_p_drivers < DRIVERS_FOR_A:
                            raise Warning("Can't find CRA parent")

                        if current_p == -1:
                            raise Warning("Don't have parent recorded")

                    add_founder_info(summary_mat, results_mat, current_p, results_mat[current_p, RESULTS_MAT_ORIGIN_IDX])
                    established_adenoma = True

                    cra_s = int(results_mat[current_p, SPECIES_COL_IDX])
                    print("Setting CRA founder as", current_p, " of species=", cra_s, "with ant=",
                          summary_mat[cra_s, FOUNDER_ANT_IDX], ", S=", summary_mat[cra_s, FOUNDER_S_IDX],
                          "B=", summary_mat[cra_s, FOUNDER_B_IDX], "and originated at time", summary_mat[cra_s, ORIGIN_COL_IDX])

                    if species_sizes[A_IDX] == 0:
                        raise Warning("Created carcinoma, but no adenoma present!")

    return next_sizes, next_species_sizes, established_adenoma, established_carcinoma

@nba.njit()
def run(pop_sizes, next_pop_sizes, species_sizes, next_species_sizes, summary_mat, results_mat, iter_order):
    """
    Run the simulation

    :param pop_sizes: current sizes of each phenotype
    :param next_pop_sizes: sizes of each phenotype in the next step. Initially pop_sizes, but gets updated in this function
    :param species_sizes: current size of each "species" (epithelial, adenoma, carcinoma)
    :param next_species_sizes: sizes of each species in the next step. Initially species_sizes, but gets updated in this function
    :param summary_mat: sizes of each species at all timepoints
    :param results_mat: sizes of each phenotype at all timepoints
    :param iter_order: order in which to simulate each phenotype. Shuffled each time step

    """
    tumor_status = NO_TUMOR_STATUS
    established_adenoma = False
    established_carcinoma = False
    time_carcinoma = 0
    if RECORD_TIME:
        results_mat[0, N_INFO_COLS] = INITIAL_E_SIZE

    for time_pt in range(1, NREPS):
        # print(time_pt)

        temp_pop_sizes = pop_sizes
        pop_sizes = next_pop_sizes
        next_pop_sizes = temp_pop_sizes
        for i in range(next_pop_sizes.size):
            next_pop_sizes[i] = 0

        temp_species_sizes = species_sizes
        species_sizes = next_species_sizes
        next_species_sizes = temp_species_sizes
        for i in range(next_species_sizes.size):
            next_species_sizes[i] = 0

        if established_carcinoma:
            time_carcinoma += 1

        # Need to add mutants in random order to avoid biasing towards clones at top or bottom of results matrix
        # or favoring less antigenic clones (which have lower indices) to become adenomas or carcinomas
        np.random.shuffle(iter_order)
        next_pop_sizes, next_species_sizes, established_adenoma, established_carcinoma = step(iter_order, pop_sizes,
                                                                                              next_pop_sizes,
                                                                                              summary_mat, results_mat,
                                                                                              species_sizes, next_species_sizes,
                                                                                              time_pt, established_adenoma,
                                                                                              established_carcinoma)

        # Carcinoma existed for specidied number of steps
        if established_carcinoma and time_carcinoma >= STEPS_AS_CARCINOMA_STOP_THRESH:
            print("cacrcinoma existed for", time_carcinoma)
            tumor_status = EXTANT_STATUS
            summarize_time_pt(summary_mat, results_mat, next_pop_sizes, C_IDX, time_pt + 1)
            break

        # TUMOR WENT EXTINCT
        if established_adenoma and not established_carcinoma and species_sizes[A_IDX] == 0:
            print("adenoma went extinct at time", time_pt)
            tumor_status = EXTINCT_STATUS
            summarize_time_pt(summary_mat, results_mat, pop_sizes, E_IDX, time_pt)
            summarize_time_pt(summary_mat, results_mat, pop_sizes, A_IDX, time_pt)
            break

        if established_carcinoma and species_sizes[C_IDX] == 0:
            print("carinoma went extinct")
            summarize_time_pt(summary_mat, results_mat, pop_sizes, C_IDX, time_pt)
            tumor_status = EXTINCT_STATUS
            break

        if np.sum(species_sizes) == 0:
            print("epithelial went extinct at time", time_pt)
            summarize_time_pt(summary_mat, results_mat, pop_sizes, E_IDX, time_pt)
            tumor_status = EXTINCT_STATUS
            break

    # TUMOR WAS CONTROLLED, e.g. time ran out and stop condition not met
    if established_adenoma and time_pt + 1 >= NREPS:
        summarize_time_pt(summary_mat, results_mat, next_pop_sizes, E_IDX, time_pt + 1)
        summarize_time_pt(summary_mat, results_mat, next_pop_sizes, A_IDX, time_pt + 1)
        tumor_status = CONTROLLED_STATUS

    if established_carcinoma and time_pt + 1 >= NREPS:
        summarize_time_pt(summary_mat, results_mat, next_pop_sizes, C_IDX, time_pt + 1)
        tumor_status = CONTROLLED_STATUS

    final_size = np.sum(next_species_sizes)

    return time_pt, summary_mat, results_mat, tumor_status, established_adenoma, established_carcinoma, final_size

def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(description='Run simulation of tumor-immune eco-evolutionary dynamics')
    parser.add_argument("-dst_dir", '--dst_dir', type=str, help='directory where to save results')
    parser.add_argument("-f_prefix", '--f_prefix', type=str, help='file prefix for all results')
    parser.add_argument("-bp_mutation_rate", '--bp_mutation_rate', type=float, help='rate at which a clone mutates')
    parser.add_argument("-blockade_protection", '--blockade_protection', type=float, help='amount of protection from immune blockade')
    parser.add_argument("-immunosuppression", '--immunosuppression', type=float, help='amount of immunosuppression provided by immunosuppressive cells')
    parser.add_argument("-tumor_id", '--tumor_id', type=int, help='tumor id number')
    parser.add_argument("-record_time", '--record_time', type=str2bool, help='Whether or not time points should be recorded')

    args = vars(parser.parse_args())
    print(args)

    dst_dir = args["dst_dir"]
    f_prefix = args["f_prefix"]
    BLOCKADE_PROTECTION = args["blockade_protection"]
    IMMUNOSUPPRESSION = args["immunosuppression"]
    BP_MUTATION_RATE = args["bp_mutation_rate"]
    tumor_id = args["tumor_id"]
    RECORD_TIME = args["record_time"]


    CLONE_MUTATION_RATE = 1 - stats.binom.pmf(k=0, n=BP_IN_EXOME, p=BP_MUTATION_RATE)

    no_tumor_dir = os.path.join(dst_dir, "no_tumors") # Niether an adenoma or carcinoma evolved
    extant_carcinoma_dir = os.path.join(dst_dir, "extant_carcinomas")
    extinct_carcinoma_dir = os.path.join(dst_dir, "extinct_carcinomas")
    controlled_carcinoma_dir = os.path.join(dst_dir, "controlled_carcinomas")
    extinct_adenoma_dir = os.path.join(dst_dir, "extinct_adinomas")
    controlled_adenoma_dir = os.path.join(dst_dir, "controlled_adenomas")
    extinct_epithelial_dir = os.path.join(dst_dir, "extinct_epithelial")

    for tid in [tumor_id]:
        # tracemalloc.start()
        print("Starting sim", str(tid), "with mutation rate", BP_MUTATION_RATE, "blockade", BLOCKADE_PROTECTION, "suppression", IMMUNOSUPPRESSION)
        start = time.time()
        results_mat, summary_mat = init_results_mat()


        print("Starting matrix=", results_mat.shape)
        species_sizes = np.array([INITIAL_E_SIZE, 0, 0])
        pop_sizes = np.zeros(results_mat.shape[0])
        pop_sizes[0] = INITIAL_E_SIZE

        next_species_sizes = species_sizes.copy()
        next_pop_sizes = pop_sizes.copy()

        iter_order = np.arange(0, results_mat.shape[0])
        last_time_pt, summary_mat, results, tumor_status, established_adenoma, established_carcinoma, final_size = run(
            pop_sizes, next_pop_sizes, species_sizes, next_species_sizes, summary_mat, results_mat, iter_order)

        end = time.time()
        print("completed", last_time_pt + 1, " repetitions, in", (end - start)/60, "minutes. Result is", tumor_status)
        # print("== Print memory allocation information == If no leak, values should be the same")
        print(nba.runtime.rtsys.get_allocation_stats())

        summary_mat = summary_mat[summary_mat[:, ORIGIN_COL_IDX] >= 0, ]

        summary_df = pd.DataFrame(summary_mat)
        current_colnames = list(summary_df)
        summary_df = summary_df.rename(columns=summary_colname_dict)
        summary_df.insert(0, "blockade_param", BLOCKADE_PROTECTION)
        summary_df.insert(0, "suppression_param", IMMUNOSUPPRESSION)
        summary_df.insert(0, "bp_mutation_rate", BP_MUTATION_RATE)
        summary_df.insert(0, "tumor_id", tid)

        if tumor_status == NO_TUMOR_STATUS:
            if not os.path.exists(no_tumor_dir):
                os.makedirs(no_tumor_dir)
            print("No adenoma or carcinoma. Next")
            summary_df.insert(0, "status", "no_tumor")
            f_out = os.path.join(no_tumor_dir, str(tid) + "_no_tumor_summary.csv")
            summary_df.to_csv(f_out, index=False)

        if tumor_status == EXTANT_STATUS:
            if not os.path.exists(extant_carcinoma_dir):
                os.makedirs(extant_carcinoma_dir)
            print("Extant carcinoma. Ending")
            summary_df.insert(0, "status", "extant_carcinoma")
            f_out = os.path.join(extant_carcinoma_dir, str(tid) + "_extant_carcinoma_summary.csv")
            summary_df.to_csv(f_out, index=False)

        if tumor_status == EXTINCT_STATUS and final_size == 0:
            # Epithelial can go extinct. Happens when mutation rate so high that root population is depleted. Muller's ratchet
            if not os.path.exists(extinct_epithelial_dir):
                os.makedirs(extinct_epithelial_dir)
            print("Extinct epithelial. Next")
            summary_df.insert(0, "status", "extinct_epithelial")
            f_out = os.path.join(extinct_epithelial_dir, str(tid) + "_extinct_epithelial_summary.csv")
            summary_df.to_csv(f_out, index=False)

        if established_adenoma and not established_carcinoma:
            if tumor_status == EXTINCT_STATUS:
                if not os.path.exists(extinct_adenoma_dir):
                    os.makedirs(extinct_adenoma_dir)
                print("Extinct adenoma. Next")
                summary_df.insert(0, "status", "extinct_adenoma")
                f_out = os.path.join(extinct_adenoma_dir, str(tid) + "_extinct_adenoma_summary.csv")
                summary_df.to_csv(f_out, index=False)

            elif tumor_status == CONTROLLED_STATUS:
                if not os.path.exists(controlled_adenoma_dir):
                    os.makedirs(controlled_adenoma_dir)
                print("Controlled adenoma. Next")
                summary_df.insert(0, "status", "controlled_adenoma")
                f_out = os.path.join(controlled_adenoma_dir, str(tid) + "_controlled_adenoma_summary.csv")
                summary_df.to_csv(f_out, index=False)

        elif established_carcinoma:
            if tumor_status == EXTINCT_STATUS:
                if not os.path.exists(extinct_carcinoma_dir):
                    os.makedirs(extinct_carcinoma_dir)
                print("Extinct carcinoma. Next")
                summary_df.insert(0, "status", "extinct_carcinoma")
                f_out = os.path.join(extinct_carcinoma_dir, str(tid) + "_extinct_carcinoma_summary.csv")
                summary_df.to_csv(f_out, index=False)

            elif tumor_status == CONTROLLED_STATUS:
                if not os.path.exists(controlled_carcinoma_dir):
                    os.makedirs(controlled_carcinoma_dir)
                print("Controlled carcinoma. Next")
                summary_df.insert(0, "status", "controlled_carcinoma")
                f_out = os.path.join(controlled_carcinoma_dir, str(tid) + "_controlled_carcinoma_summary.csv")
                summary_df.to_csv(f_out, index=False)


        if RECORD_TIME:
            clone_keep_idx = np.where(np.any(results[:, N_INFO_COLS:] > 0, axis=1))[0]

            results = results[clone_keep_idx, :]
            results = results[:, :N_INFO_COLS + last_time_pt//RECORD_FREQ]

            results_df = pd.DataFrame(results)
            current_colnames = list(results_df)
            time_colnames = np.array(current_colnames[N_INFO_COLS:]) - N_INFO_COLS
            time_colnames *= RECORD_FREQ

            new_colnames = INFO_COLS.copy()
            new_colnames.extend(time_colnames)

            rename_dict = {current_colnames[j]: new_colnames[j] for j in range(len(new_colnames))}
            results_df = results_df.rename(columns=rename_dict)
            results_df.insert(0, "blockade_param", BLOCKADE_PROTECTION)
            results_df.insert(0, "suppression_param", IMMUNOSUPPRESSION)
            results_df.insert(0, "bp_mutation_rate", BP_MUTATION_RATE)
            results_df.insert(0, "tumor_id", tid)

            if tumor_status == NO_TUMOR_STATUS:
                if not os.path.exists(no_tumor_dir):
                    os.makedirs(no_tumor_dir)
                print("No adenoma or carcinoma. Next")
                results_df.insert(0, "status", "no_tumor")
                f_out = os.path.join(no_tumor_dir, str(tid) + "_no_tumor.csv")
                results_df.to_csv(f_out, index=False)

            if tumor_status == EXTANT_STATUS:
                if not os.path.exists(extant_carcinoma_dir):
                    os.makedirs(extant_carcinoma_dir)
                print("Extant carcinoma. Ending")
                results_df.insert(0, "status", "extant_carcinoma")
                f_out = os.path.join(extant_carcinoma_dir, str(tid) + "_extant_carcinoma.csv")
                results_df.to_csv(f_out, index=False)

            if tumor_status == EXTINCT_STATUS and final_size == 0:
                # Epithelial can go extinct. Happens when mutation rate so high that root population is depleted. Muller's ratchet
                if not os.path.exists(extinct_epithelial_dir):
                    os.makedirs(extinct_epithelial_dir)
                print("Extinct epithelial. Next")
                results_df.insert(0, "status", "extinct_epithelial")
                f_out = os.path.join(extinct_epithelial_dir, str(tid) + "_extinct_epithelial.csv")
                results_df.to_csv(f_out, index=False)

            if established_adenoma and not established_carcinoma:
                if tumor_status == EXTINCT_STATUS:
                    if not os.path.exists(extinct_adenoma_dir):
                        os.makedirs(extinct_adenoma_dir)
                    print("Extinct adenoma. Next")
                    results_df.insert(0, "status", "extinct_adenoma")
                    f_out = os.path.join(extinct_adenoma_dir, str(tid) + "_extinct_adenoma.csv")
                    results_df.to_csv(f_out, index=False)

                elif tumor_status == CONTROLLED_STATUS:
                    if not os.path.exists(controlled_adenoma_dir):
                        os.makedirs(controlled_adenoma_dir)
                    print("Controlled adenoma. Next")
                    results_df.insert(0, "status", "controlled_adenoma")
                    f_out = os.path.join(controlled_adenoma_dir, str(tid) + "_controlled_adenoma.csv")
                    results_df.to_csv(f_out, index=False)

            elif established_carcinoma:
                if tumor_status == EXTINCT_STATUS:
                    if not os.path.exists(extinct_carcinoma_dir):
                        os.makedirs(extinct_carcinoma_dir)
                    print("Extinct carcinoma. Next")
                    results_df.insert(0, "status", "extinct_carcinoma")
                    f_out = os.path.join(extinct_carcinoma_dir, str(tid) + "_extinct_carcinoma.csv")
                    results_df.to_csv(f_out, index=False)

                elif tumor_status == CONTROLLED_STATUS:
                    if not os.path.exists(controlled_carcinoma_dir):
                        os.makedirs(controlled_carcinoma_dir)
                    print("Controlled carcinoma. Next")

                    results_df.insert(0, "status", "controlled_carcinoma")
                    f_out = os.path.join(controlled_carcinoma_dir, str(tid) + "_controlled_carcinoma.csv")
                    results_df.to_csv(f_out, index=False)
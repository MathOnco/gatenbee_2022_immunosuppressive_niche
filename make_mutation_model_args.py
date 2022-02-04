import numpy as np
import os

MAX_SIMS = 10  # set to 100 in paper
RUN_ON_CLUSTER = False

if RUN_ON_CLUSTER:
    python_interp = ""
    arg_file_name = "mutation_model.txt"
    pyscript = ""
    cluster_pyscript = "mutation_model.py"
else:
    python_interp = "/anaconda3/bin/python3.7"

    arg_file_name = "mutation_model.sh"
    pyscript = "mutation_model.py"

argfile = open(arg_file_name, 'w')
if not RUN_ON_CLUSTER:
    argfile.writelines("#!/usr/bin/env bash\n")

# BP mutation rates from Table 1 in https://www.nature.com/articles/s41467-020-14844-6/tables/1
MED_BP_MUTATION_RATE = 2.91*(10**-9)
BP_IN_EXOME = (45-18)*(10**6)  # 45Mb - 18Mb of potential synonymous mutations http://www.nature.com/nature/journal/v536/n7616/full/nature19057.html?foxtrotcallback=true
N_VALS_FOR_EACH_STRATEGY = 25
MAX_SUPPRESSION = 1.0
MAX_PROTECTION = 1.0
RECORD_TIME = True  # Whether or not to record time points
proection_vals = np.linspace(0, MAX_PROTECTION, N_VALS_FOR_EACH_STRATEGY)
suppression_vals = np.linspace(0, MAX_SUPPRESSION, N_VALS_FOR_EACH_STRATEGY)
mutation_vals = np.array([MED_BP_MUTATION_RATE])

sim_numbers = list(range(MAX_SIMS))
param_array = np.array(np.meshgrid(proection_vals, suppression_vals, mutation_vals)).T.reshape(-1, 3) ### [get lucky ant, get smart ant, get smart protection, get_smart_suppression]
for i in range(param_array.shape[0]):
    p = param_array[i, 0]
    s = param_array[i, 1]
    m = param_array[i, 2]

    f_prefix = f"{np.around(p, 3)}p_{np.around(s, 3)}s"
    dst_dir = os.path.join(os.getcwd(), "mutation_model_results", f_prefix)
    all_completed_sims = []
    if os.path.exists(dst_dir):
        for pdir, sdir, fl in os.walk(dst_dir):
            if len(fl) > 0:
                completed_ids = [int(f.split("_")[0]) for f in fl if f.endswith("summary.csv")]
                all_completed_sims.extend(completed_ids)

    missing_sims = list(set(sim_numbers) - set(all_completed_sims))
    if len(missing_sims) == 0:
        continue

    for j in missing_sims:
        arg_line = " ".join([python_interp,
                             pyscript,
                             "-dst_dir", dst_dir,
                             "-f_prefix", f_prefix,
                             "-bp_mutation_rate", str(m),
                             "-blockade_protection", str(p),
                             "-immunosuppression", str(s),
                             "-tumor_id", str(j),
                             "-record_time", str(RECORD_TIME),
                             "\n"
                             ])
        # print(arg_line)
        argfile.writelines(arg_line)

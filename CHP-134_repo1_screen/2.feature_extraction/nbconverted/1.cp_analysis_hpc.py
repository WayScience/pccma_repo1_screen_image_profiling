#!/usr/bin/env python
# coding: utf-8

# # Perform segmentation and feature extraction using CellProfiler
# 
# NOTE: Plates with IDs that have spaces will still have spaces in the paths but the respective LoadData CSV with illum function paths will not have spaces due to `pe2loaddata` bug. This means that the LoadData from the previous IC module will be different (includes space).

# ## Import libraries

# In[1]:


import argparse
import pathlib
import pprint

import sys

sys.path.append("../../utils")
import cp_parallel

# check if in a jupyter notebook
try:
    cfg = get_ipython().config
    in_notebook = True
except NameError:
    in_notebook = False


# ## Set paths and variables

# In[2]:


#  directory where loaddata CSVs are located within the folder
loaddata_dir = pathlib.Path("./loaddata_csvs/").resolve(strict=True)

if not in_notebook:
    print("Running as script")

    parser = argparse.ArgumentParser(
        description="CellProfiler segmentation and feature extraction"
    )

    parser.add_argument(
        "--input_csv",
        type=str,
        required=True,
        help="Path to the LoadData CSV file to process images",
    )
    parser.add_argument(
        "--first_image_set",
        type=int,
        required=False,
        help="First CellProfiler image set to process",
    )
    parser.add_argument(
        "--last_image_set",
        type=int,
        required=False,
        help="Last CellProfiler image set to process",
    )
    parser.add_argument(
        "--batch_label",
        type=str,
        required=False,
        help="Label for a row batch, e.g. row_A",
    )

    args = parser.parse_args()

    loaddata_csv = pathlib.Path(args.input_csv).resolve(strict=True)
    first_image_set = args.first_image_set
    last_image_set = args.last_image_set
    batch_label = args.batch_label

else:
    print("Running in a notebook")

    loaddata_csv = pathlib.Path(
        f"{loaddata_dir}/Assay_Plate_1_3_loaddata_with_illum.csv"
    ).resolve(strict=True)
    first_image_set = None
    last_image_set = None
    batch_label = None

# set the run type for the parallelization
run_name = "analysis"

# set path for CellProfiler pipeline
path_to_pipeline = pathlib.Path("./pipeline/analysis_CHP-134.cppipe").resolve(
    strict=True
)

# set main output dir for all plates if it doesn't exist
output_dir = pathlib.Path("./sqlite_outputs")
output_dir.mkdir(exist_ok=True)


# ## Create dictionary to process data

# In[ ]:


# Extract name from LoadData CSV path (drop loaddata suffix to avoid issues getting plate name)
name = loaddata_csv.stem
if name.endswith("_loaddata_with_illum"):
    name = name[: -len("_loaddata_with_illum")]

path_to_output = output_dir / name
if batch_label:
    # if batch label is provided, add it to the output path to keep batches separate per plate
    path_to_output = path_to_output / batch_label

# create plate info dictionary with all parts of the CellProfiler CLI command to run in parallel
plate_info_dictionary = {
    name: {
        "path_to_loaddata": loaddata_csv,
        "path_to_output": path_to_output,
        "path_to_pipeline": path_to_pipeline,
    }
}

if first_image_set is not None:
    plate_info_dictionary[name]["first_image_set"] = first_image_set
if last_image_set is not None:
    plate_info_dictionary[name]["last_image_set"] = last_image_set

# view the dictionary to assess that all info is added correctly
pprint.pprint(plate_info_dictionary, indent=4)


# ## Perform segmentation and feature extraction (analysis)
# 
# Note: This code cell was not ran as we prefer to perform CellProfiler processing tasks via `sh` file (bash script) which is more stable.

# In[4]:


cp_parallel.run_cellprofiler_parallel(
    plate_info_dictionary=plate_info_dictionary, run_name=run_name
)


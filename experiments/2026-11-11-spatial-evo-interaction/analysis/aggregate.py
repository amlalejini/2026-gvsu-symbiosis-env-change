'''
Aggregate data.
This script generates the following output files:
- A summary file with one line per-replicate.
'''

import argparse
import os
import sys
import pathlib
from scipy.stats import entropy
import statistics

# Add scripts directory to path, import utilities from scripts directory.
sys.path.append(
    os.path.join(
        pathlib.Path(os.path.dirname(os.path.abspath(__file__))).parents[2],
        "scripts"
    )
)
import utilities as utils

run_identifier = "RUN_"

# Run configuration fields to keep as fields in summary output file.
run_cfg_fields_summary = {
    "CYCLES_PER_UPDATE",
    "ENABLE_HEALTH",
    "ENABLE_NUTRIENT",
    "ENABLE_STRESS",
    "EVENTS_CFG_PATH",
    "FIND_NEIGHBOR_HOST_ATTEMPTS",
    "HEALTH_INTERACTION_CHANCE",
    "HEALTH_TYPE",
    "HORIZ_TRANS",
    "HORIZONTAL_TRANSMISSION_COMPATIBILITY_MODE",
    "HOST_AGE_MAX",
    "HOST_INT",
    "HOST_MIN_CYCLES_BEFORE_REPRO",
    "HOST_ONLY_FIRST_TASK_CREDIT",
    "HOST_PROGRAM_PATH",
    "HOST_REPRO_RES",
    "INIT_POP_SIZE",
    "MUTATION_RATE",
    "MUTATION_SIZE",
    "MUTUALIST_CYCLE_DONATE_MULTIPLIER",
    "MUTUALIST_CYCLE_GAIN_PROP",
    "OUSTING",
    "PARASITE_BASE_CYCLE_PROP",
    "PARASITE_CYCLE_LOSS_PROP",
    "PARASITE_CYCLE_STEAL_MULTIPLIER",
    "SEED",
    "SGP_MUT_PER_BIT_RATE",
    "SPATIAL_STRUCT_CFG_PATH",
    "SPATIAL_STRUCT_LOAD_MODE",
    "SPATIAL_STRUCT_MODE",
    "START_MOI",
    "SYM_AGE_MAX",
    "SYM_HORIZ_TRANS_RES",
    "SYM_INT",
    "SYM_MIN_CYCLES_BEFORE_REPRO",
    "SYM_ONLY_FIRST_TASK_CREDIT",
    "SYM_PROGRAM_PATH",
    "SYM_VERT_TRANS_RES",
    "TASK_ENV_CFG_PATH",
    "TASK_IO_BANK_SIZE",
    "TASK_IO_UNIQUE_OUTPUT",
    "TASK_PROFILE_COMPATIBILITY_MODE",
    "TASK_PROFILE_MODE",
    "UPDATES",
    "VERTICAL_TRANSMISSION",
    "VT_TASK_MATCH",
    "WORLD_HEIGHT",
    "WORLD_WIDTH",
    "SPATIAL_STRUCT_CFG_PATH"
}

sym_int_vals_fields_summary = {
    "mean_intval"
}

run_cfg_fields_time_series = {
    "SEED",
    "VERTICAL_TRANSMISSION",
    "SPATIAL_STRUCT_MODE",
    "SPATIAL_STRUCT_CFG_PATH"
}

org_counts_fields_time_series = {
    "host_count",
    "hosted_sym_count"
}

cur_update_info_fields_time_series = {
    "NOT_in_host_profile_counts",
    "NAND_in_host_profile_counts",
    "OR_NOT_in_host_profile_counts",
    "AND_in_host_profile_counts",
    "OR_in_host_profile_counts",
    "AND_NOT_in_host_profile_counts",
    "NOR_in_host_profile_counts",
    "XOR_in_host_profile_counts",
    "EQU_in_host_profile_counts",
    "NOT_in_sym_profile_counts",
    "NAND_in_sym_profile_counts",
    "OR_NOT_in_sym_profile_counts",
    "AND_in_sym_profile_counts",
    "OR_in_sym_profile_counts",
    "AND_NOT_in_sym_profile_counts",
    "NOR_in_sym_profile_counts",
    "XOR_in_sym_profile_counts",
    "EQU_in_sym_profile_counts",
    "NOT_in_host_parent_org_counts",
    "NAND_in_host_parent_org_counts",
    "OR_NOT_in_host_parent_org_counts",
    "AND_in_host_parent_org_counts",
    "OR_in_host_parent_org_counts",
    "AND_NOT_in_host_parent_org_counts",
    "NOR_in_host_parent_org_counts",
    "XOR_in_host_parent_org_counts",
    "EQU_in_host_parent_org_counts",
    "NOT_in_sym_parent_org_counts",
    "NAND_in_sym_parent_org_counts",
    "OR_NOT_in_sym_parent_org_counts",
    "AND_in_sym_parent_org_counts",
    "OR_in_sym_parent_org_counts",
    "AND_NOT_in_sym_parent_org_counts",
    "NOR_in_sym_parent_org_counts",
    "XOR_in_sym_parent_org_counts",
    "EQU_in_sym_parent_org_counts",
    "host_parent_entropy_task_sets",
    "host_current_entropy_task_sets",
    "host_parent_num_task_sets",
    "host_current_num_task_sets",
    "sym_parent_entropy_task_sets",
    "sym_current_entropy_task_sets",
    "sym_current_num_task_sets",
    "sym_parent_num_task_sets",
    "host_sym_perfect_matches_total",
    "host_sym_any_matches_total",
    "CurUpdate_sym_mean_generations",
    "CurUpdate_host_mean_generations"
}

tasks_file_fields_time_series = {
    "host_task_NOT",
    "host_task_NAND",
    "host_task_OR_NOT",
    "host_task_AND",
    "host_task_OR",
    "host_task_AND_NOT",
    "host_task_NOR",
    "host_task_XOR",
    "host_task_EQU",
    "sym_task_NOT",
    "sym_task_NAND",
    "sym_task_OR_NOT",
    "sym_task_AND",
    "sym_task_OR",
    "sym_task_AND_NOT",
    "sym_task_NOR",
    "sym_task_XOR",
    "sym_task_EQU"
}

transmission_rates_fields_time_series = {
    "attempts_horiztrans",
    "successes_horiztrans",
    "attempts_verttrans",
    "successes_verttrans"
}

def extract_summary_data(data, target_update, fields, prefix=None):
        info = {}

        # Grab the data line that matches the target update for this run
        summary_data = [
            line
            for line in data
            if (target_update is None) or (int(line["update"]) == target_update)
        ][-1]

        # Add specified fields to run summary data
        for field in summary_data:
            if field in fields:
                if prefix is None:
                    info[field] = summary_data[field]
                else:
                    info[f"{prefix}_{field}"] = summary_data[field]

        return info

def add_time_series_info(
    time_series_data,
    run_data,
    fields,
    prefix = None
):
    # For each relevant line in run data, add relevant fields to time_series_data
    for line in run_data:
        # Skip over updates we don't want to sample
        line_update = int(line["update"])
        if not line_update in time_series_data:
            continue
        for field in line:
            if field in fields:
                if prefix is None:
                    time_series_data[line_update][field] = line[field]
                else:
                    time_series_data[line_update][f"{prefix}_{field}"] = line[field]


def main():
    parser = argparse.ArgumentParser(description = "Run submission script.")
    parser.add_argument("--data_dir", type=str, help="Where is the base output directory for each run?")
    parser.add_argument("--dump_dir", type=str, help="Where to dump this?", default=".")
    parser.add_argument("--summary_update", type=int, help="Update to pull summary data for?")
    parser.add_argument("--time_series_units", type=str, default="interval", choices=["interval", "total"], help="Unit for resolution of time series")
    parser.add_argument("--time_series_resolution", type=int, default=1, help="What resolution should we collect time series data at?")
    # -- Spatial data aggregation options --
    parser.add_argument("--spatial_agg", type=bool, default=False, help="Also aggregate spatial data?")
    parser.add_argument("--spatial_struct_dir", type=str, default=None, help="If aggregating spatial data, where to find spatial structure node info files?")

    args = parser.parse_args()
    data_dir = args.data_dir
    dump_dir = args.dump_dir
    target_update = args.summary_update
    time_series_units = args.time_series_units
    time_series_resolution = args.time_series_resolution
    aggregate_spatial_data = args.spatial_agg
    spatial_struct_dir = args.spatial_struct_dir

    if not os.path.exists(data_dir):
        print("Unable to find data directory.")
        exit(-1)

    # Verify time series resolution >= 1
    if time_series_resolution < 1:
        print("Time series resolution must be >= 1")
        exit(-1)

    if aggregate_spatial_data:
        if (spatial_struct_dir is None) or (not os.path.exists(spatial_struct_dir)):
            print(f"Unable to find spatial_struct_dir: {spatial_struct_dir}")
            exit(-1)

    utils.mkdir_p(dump_dir)

    # Aggregate run directories.
    run_dirs = [run_dir for run_dir in os.listdir(data_dir) if run_identifier in run_dir]
    print(f"Found {len(run_dirs)} run directories.")

    # Lazy-read spatial structure node info files.
    # Format {structure id: {loc_id: {region: ..., connection_point: ...}, }, ..}
    spatial_struct_node_info = {}
    spatial_regions_by_struct = {}

    # Create file to hold time series data
    time_series_content = []    # This will hold all the lines to write out for a single run; written out for each run.
    time_series_header = None   # Holds the time series file header (verified for consistency across runs)
    time_series_fpath = os.path.join(dump_dir, f"time_series.csv")

    with open(time_series_fpath, "w") as fp:
        fp.write("")

    # spatial_data_content holds per-location data
    spatial_data_out_path = os.path.join(dump_dir, "spatial_data.csv")
    spatial_data_header = None
    spatial_data_content = []
    # spatial_region_data_content holds per-region data
    spatial_region_data_out_path = os.path.join(dump_dir, "spatial_region_data.csv")
    spatial_region_data_header = None
    spatial_region_content = []
    if aggregate_spatial_data:
        # Create file to hold spatial data
        with open(spatial_data_out_path, "w") as fp:
            fp.write("")
        # Create file to hold spatial region data
        with open(spatial_region_data_out_path, "w") as fp:
            fp.write("")

    task_discovery_out_path = os.path.join(dump_dir, "task_discovery.csv")
    task_discovery_header = None
    task_discovery_content = []
    with open(task_discovery_out_path, "w") as fp:
        fp.write("")

    # For each run directory...
    # summary_header = None
    summary_content_lines = []
    sym_interaction_values_content_lines = []
    incomplete_runs = []
    for run_dir_i in range(len(run_dirs)):
        run_dir = run_dirs[run_dir_i]
        print(f"...({run_dir_i + 1}/{len(run_dirs)}) aggregating from {run_dir}")
        run_path = os.path.join(data_dir, run_dir)

        run_summary_info = {} # Hold summary information about this run.
        sym_int_vals_info = {}
        time_series_info = {} # Hold time series information. Indexed by update.
        spatial_out_info = {} # Hold spatial data output info. Indexed by location.
        task_discovery_info = {} # Hold task discovery information. Indexed by task.
        region_out_info = {} # Hold information by region. Indexed by region name.
        ########################################
        # Extract run parameters
        ########################################
        run_cfg_path = os.path.join(run_path, "output", "run_config.csv")
        if not os.path.isfile(run_cfg_path):
            print("Run did not finish, skipping")
            incomplete_runs.append(run_dir)
            continue

        run_cfg_data = utils.read_csv(run_cfg_path)
        run_params = {}
        for line in run_cfg_data:
            param = line["parameter"]
            value = line["value"]
            run_params[param] = value
            # Add a subset of parameters to summary information for this run.
            if param in run_cfg_fields_summary:
                run_summary_info[param] = value
                sym_int_vals_info[param] = value

        max_pop_size = 0

        spatial_struct_mode = run_params["SPATIAL_STRUCT_MODE"]
        if spatial_struct_mode in {"well-mixed", "grid"}:
            max_pop_size = int(run_params["WORLD_WIDTH"]) * int(run_params["WORLD_HEIGHT"])
        elif spatial_struct_mode == "load":
            spatial_struct_fname = run_params["SPATIAL_STRUCT_CFG_PATH"]
            spatial_struct_path = os.path.join(run_path, spatial_struct_fname)
            with open(spatial_struct_path, "r") as fp:
                max_pop_size = len(fp.readlines())
        else:
            print("Unknown spatial structure mode!")
            exit()

        run_summary_info["max_pop_size"] = max_pop_size
        sym_int_vals_info["max_pop_size"] = max_pop_size

        ########################################
        # Extract data from OrganismCounts.csv
        ########################################
        org_counts_path = os.path.join(run_path, "output", "OrganismCounts_data.csv")
        org_counts_data = utils.read_csv(org_counts_path)

        # --- Analyze updates represented, setup time series info --
        # Grab list of updates represented in data
        updates = [int(row["update"]) for row in org_counts_data]
        if len(updates) == 0:
            continue

        # Did run finish with respect to target update?
        run_finished_target = target_update in updates

        # Extract time series updates (only if run reached target)
        time_series_updates = utils.filter_time_points(
            updates,
            method = time_series_units,
            resolution = time_series_resolution
        ) if run_finished_target else []
        time_series_updates = set(time_series_updates)
        # Add run cfg information to time_series info
        time_series_info = {
            update:{field:run_params[field] for field in run_cfg_fields_time_series}
            for update in time_series_updates
        }
        for update in time_series_updates:
            time_series_info[update]["update"] = update

        run_target_update = utils.nearest(target_update, updates)
        run_summary_info["update"] = run_target_update
        run_summary_info["reached_target_update"] = run_finished_target
        sym_int_vals_info["update"] = run_target_update
        sym_int_vals_info["reached_target_update"] = run_finished_target
        # ---

        # Extract summary info
        org_counts_fields = set(org_counts_data[0].keys())
        org_counts_fields.remove("update")
        run_summary_info.update(
            extract_summary_data(
                data = org_counts_data,
                target_update = run_target_update,
                fields = org_counts_fields,
                prefix = "OrgCounts"
            )
        )

        # Extract time series info
        if run_finished_target:
            add_time_series_info(
                time_series_data = time_series_info,
                run_data = org_counts_data,
                fields = org_counts_fields_time_series,
                prefix = "OrgCounts"
            )

        del org_counts_data
        ########################################
        # Extract data from CurrentUpdateInfo.csv
        ########################################
        cur_update_info_path = os.path.join(run_path, "output", "CurrentUpdateInfo_data.csv")
        cur_update_info_data = utils.read_csv(cur_update_info_path)

        # Extract summary info
        cur_update_info_fields = set(cur_update_info_data[0].keys())
        cur_update_info_fields.remove("update")
        run_summary_info.update(
            extract_summary_data(
                data = cur_update_info_data,
                target_update = run_target_update,
                fields = cur_update_info_fields,
                prefix = "CurUpdate"
            )
        )

        # Extract time series info
        if run_finished_target:
            add_time_series_info(
                time_series_data = time_series_info,
                run_data = cur_update_info_data,
                fields = cur_update_info_fields_time_series,
                prefix = "CurUpdate"
            )

        del cur_update_info_data
        ########################################
        # Extract data from Tasks.csv
        ########################################
        tasks_path = os.path.join(run_path, "output", "Tasks_data.csv")
        tasks_data = utils.read_csv(tasks_path)

        # Extract summary info
        tasks_fields = set(tasks_data[0].keys())
        tasks_fields.remove("update")
        run_summary_info.update(
            extract_summary_data(
                data = tasks_data,
                target_update = run_target_update,
                fields = tasks_fields,
                prefix = "Tasks"
            )
        )

        # Extract time series info
        if run_finished_target:
            add_time_series_info(
                time_series_data = time_series_info,
                run_data = tasks_data,
                fields = tasks_file_fields_time_series,
                prefix = "Tasks"
            )

        del tasks_data
        ########################################
        # Extract data from TransmissionRates.csv
        ########################################
        transmission_rates_path = os.path.join(run_path, "output", "TransmissionRates_data.csv")
        transmission_rates_data = utils.read_csv(transmission_rates_path)

        # Extract summary info
        transmission_rates_fields = set(transmission_rates_data[0].keys())
        transmission_rates_fields.remove("update")
        run_summary_info.update(
            extract_summary_data(
                data = transmission_rates_data,
                target_update = run_target_update,
                fields = transmission_rates_fields,
                prefix = "TransmissionRates"
            )
        )

        # Extract time series info
        if run_finished_target:
            add_time_series_info(
                time_series_data = time_series_info,
                run_data = transmission_rates_data,
                fields = transmission_rates_fields_time_series,
                prefix = "TransmissionRates"
            )

        del transmission_rates_data
        ########################################
        # Extract data from SymbiontInteractionValues.csv
        ########################################
        # update,mean_intval,count,Hist_-1,Hist_-0.9,Hist_-0.8,Hist_-0.7,Hist_-0.6,Hist_-0.5,Hist_-0.4,Hist_-0.3,Hist_-0.2,Hist_-0.1,Hist_0.0,Hist_0.1,Hist_0.2,Hist_0.3,Hist_0.4,Hist_0.5,Hist_0.6,Hist_0.7,Hist_0.8,Hist_0.9
        sym_int_vals_path = os.path.join(run_path, "output", "SymbiontInteractionValues_data.csv")
        sym_int_vals_data = utils.read_csv(sym_int_vals_path)

        # Update run summary info
        run_summary_info.update(
            extract_summary_data(
                data = sym_int_vals_data,
                target_update = run_target_update,
                fields = sym_int_vals_fields_summary,
                prefix = "sym_int_vals"
            )
        )

        # Update interaction value file
        sym_int_fields = set(sym_int_vals_data[0].keys())
        sym_int_fields.remove("update")
        sym_int_vals_info.update(
            extract_summary_data(
                data = sym_int_vals_data,
                target_update = run_target_update,
                fields = sym_int_fields
            )
        )
        del sym_int_vals_data
        ########################################
        # Task discovery data
        ########################################
        task_discovery_path = os.path.join(run_path, "output", "TaskDiscovery_data.csv")
        task_discovery_data = utils.read_csv(task_discovery_path)
        # Task discovery file is output at end of run (no need to filter lines or keep update).
        # Also, because of quirks in output timing, update is one off.
        task_discovery_info = {
            line["task_name"]: {
                c: line[c] for c in line if c != "output_update"
            } for line in task_discovery_data
        }
        # Add run config + task complexity info to each task discovery line
        for task in task_discovery_info:
            task_discovery_info[task].update(
                {field:run_params[field] for field in run_cfg_fields_time_series}
            )
            task_discovery_info[task]["task_complexity"] = str(utils.get_task_complexity(task))

        del task_discovery_data
        ########################################
        # Spatial structure aggregation?
        ########################################
        if aggregate_spatial_data:
            # Get spatial structure file name
            spatial_struct_path = pathlib.Path(run_params["SPATIAL_STRUCT_CFG_PATH"])
            spatial_struct_mode = run_params["SPATIAL_STRUCT_MODE"]
            spatial_struct_name = None
            # Need different behavior based on spatial structure mode
            if spatial_struct_mode == "well-mixed":
                # Use well-mixed
                spatial_struct_name = "well-mixed"
            elif spatial_struct_mode == "grid":
                spatial_struct_name = "grid"
            elif spatial_struct_mode == "load":
                # Get just the file name
                spatial_struct_name = spatial_struct_path.stem
            else:
                print("Unrecognized SPATIAL_STRUCT_MODE parameter:", run_params["SPATIAL_STRUCT_MODE"])

            if not (spatial_struct_name in spatial_struct_node_info):
                # If in load mode, read info file.
                # Otherwise, need to generate info.
                if run_params["SPATIAL_STRUCT_MODE"] == "load":
                    info_fname = f"node_info_{spatial_struct_name}.csv"
                    info_path = os.path.join(spatial_struct_dir, info_fname)
                    if not os.path.exists(info_path):
                        print(f"Failed to find node info file for this spatial structure: {info_path}")
                        exit(-1)
                    node_info_data = utils.read_csv(info_path)
                    # Key by loc id, add to node info dictionary + group locations by region
                    spatial_struct_node_info[spatial_struct_name] = {}
                    spatial_regions_by_struct[spatial_struct_name] = {}
                    for line in node_info_data:
                        loc_id = line["loc_id"]
                        spatial_struct_node_info[spatial_struct_name][loc_id] = {p: line[p] for p in line}
                        # If region column, note region. Otherwise, use structure name as name.
                        region = line["region"] if "region" in line else spatial_struct_name
                        # Add any required info in case not in info file
                        if not "region" in spatial_struct_node_info[spatial_struct_name][loc_id]:
                            spatial_struct_node_info[spatial_struct_name][loc_id]["region"] = spatial_struct_name
                        if not "connection_point" in spatial_struct_node_info[spatial_struct_name][loc_id]:
                            spatial_struct_node_info[spatial_struct_name][loc_id]["connection_point"] = "0"
                        if not region in spatial_regions_by_struct[spatial_struct_name]:
                            spatial_regions_by_struct[spatial_struct_name][region] = set()
                        spatial_regions_by_struct[spatial_struct_name][region].add(loc_id)
                    # Initialize spatial regions mapping
                    del node_info_data
                elif run_params["SPATIAL_STRUCT_MODE"] == "well-mixed":
                    # Get world size
                    world_size = int(run_params["WORLD_HEIGHT"]) * int(run_params["WORLD_WIDTH"])
                    # Fill out node info manually
                    spatial_struct_node_info[spatial_struct_name] = {
                        str(loc_id): {
                            "region": "clique",
                            "connection_point": "False",
                            "loc_id": loc_id
                        } for loc_id in range(world_size)
                    }
                    # Add clique region to spatial regions by structure
                    spatial_regions_by_struct[spatial_struct_name]["clique"] = {str(loc_id) for loc_id in range(world_size)}
                elif run_params["SPATIAL_STRUCT_MODE"] == "grid":
                    # Get world size
                    world_size = int(run_params["WORLD_HEIGHT"]) * int(run_params["WORLD_WIDTH"])
                    # Fill out node info manually
                    spatial_struct_node_info[spatial_struct_name] = {
                        str(loc_id): {
                            "region": "grid",
                            "connection_point": "False",
                            "loc_id": loc_id
                        } for loc_id in range(world_size)
                    }
                    # Add grid region to spatial regions by structure
                    spatial_regions_by_struct[spatial_struct_name]["grid"] = {str(loc_id) for loc_id in range(world_size)}

            # At this point, we know spatial_struct_node_info is filled out
            node_info = spatial_struct_node_info[spatial_struct_name]
            spatial_regions = spatial_regions_by_struct[spatial_struct_name]
            # Read spatial data file
            spatial_data_path = os.path.join(run_path, "output", "Spatial_data.csv")
            spatial_data = utils.read_csv(spatial_data_path)
            # Filter spatial data down to target update
            spatial_data = [line for line in spatial_data if line["update"] == str(run_target_update)]
            # Initialize spatial output info with config info (use same as for time series)
            spatial_out_info = {
                loc_id: {field:run_params[field] for field in run_cfg_fields_time_series} for loc_id in node_info
            }
            # Add node info information + run data
            for loc_id in spatial_out_info:
                loc_region = node_info[loc_id]["region"]
                spatial_out_info[loc_id].update(
                    {
                        "update": run_target_update,
                        "location_id": loc_id,
                        "connection_point": node_info[loc_id]["connection_point"],
                        "region": loc_region
                    }
                )

            # NOTE: adjust accordingly if ever using mult-sym
            def get_endo(s):
                endo_vals = s.strip("[] ").split(" ")
                return endo_vals[-1] if len(endo_vals) > 0 else "NONE"

            def categorize_int_val(val):
                if val < -0.2:
                    return "antagonistic"
                elif val >= -0.2 and val <= 0.2:
                    return "nearly_neutral"
                elif val > 0.2:
                    return "mutualistic"
            # Initialize regional statistics
            region_out_info = {
                region: {
                    "host_antagonistic_count": 0,
                    "endosym_antagonistic_count": 0,
                    "host_nearly_neutral_count": 0,
                    "endosym_nearly_neutral_count": 0,
                    "host_mutualistic_count": 0,
                    "endosym_mutualistic_count": 0,
                    "host_present_count": 0,
                    "endosym_present_count": 0,
                    "host_sym_task_profile_compatible_count": 0,
                    "host_sym_matching_parent_task_profile_count": 0,
                    "avg_host_lineage_length": [] # Will get converted to float later
                } for region in spatial_regions
            }
            # Add configuration information to regional output information
            for region in region_out_info:
                region_out_info[region].update(
                    {field:run_params[field] for field in run_cfg_fields_time_series}
                )
                region_out_info[region]["region_size"] = len(spatial_regions[region])
                region_out_info[region]["region"] = region

            for line in spatial_data:
                loc_id = line["location_id"]
                region = node_info[loc_id]["region"]
                spatial_out_info[loc_id].update(
                    {
                        "host_present": line["host_present"],
                        "endosymbionts_present": get_endo(line["endosymbionts_present"]),
                        "host_interaction_value": line["host_interaction_value"],
                        "endosym_interaction_value": get_endo(line["endosym_interaction_values"]),
                        "host_lineage_length": line["host_lineage_length"],
                        "endosym_lineage_length": get_endo(line["endosym_lineage_length"])
                    }
                )
                host_present = line["host_present"] == "1"
                endosymbionts_present =  get_endo(line["endosymbionts_present"]) == "1"
                host_sym_compatibility = get_endo(line["host_sym_task_profile_compatibility"])
                region_out_info[region]["host_sym_task_profile_compatible_count"] += host_sym_compatibility == "1"
                region_out_info[region]["avg_host_lineage_length"].append(int(line["host_lineage_length"]))
                host_parent_task_profile = line["host_task_profile"]
                endosym_parent_task_profile = get_endo(line["endosym_task_profiles"])
                parent_task_profile_match = host_parent_task_profile == endosym_parent_task_profile
                region_out_info[region]["host_sym_matching_parent_task_profile_count"] += parent_task_profile_match
                # Categorize interaction values
                if host_present:
                    region_out_info[region]["host_present_count"] += 1
                    host_val = float(line["host_interaction_value"])
                    cat = categorize_int_val(host_val)
                    region_out_info[region][f"host_{cat}_count"] += 1
                if endosymbionts_present:
                    region_out_info[region]["endosym_present_count"] += 1
                    endo_val = float(get_endo(line["endosym_interaction_values"]))
                    cat = categorize_int_val(endo_val)
                    region_out_info[region][f"endosym_{cat}_count"] += 1

            # Calculate regional proportions
            for region in region_out_info:
                # Calculate proportion for each count type
                region_size = region_out_info[region]["region_size"]
                num_hosts = region_out_info[region]["host_present_count"]
                num_endosyms = region_out_info[region]["endosym_present_count"]
                # Host regional proportions
                region_out_info[region]["host_occupancy_regional_prop"] = region_out_info[region]["host_present_count"] / region_size
                region_out_info[region]["host_nearly_neutral_regional_prop"] = region_out_info[region]["host_nearly_neutral_count"] / num_hosts if num_hosts > 0 else 0
                region_out_info[region]["host_mutualistic_regional_prop"] = region_out_info[region]["host_mutualistic_count"] / num_hosts if num_hosts > 0 else 0
                region_out_info[region]["host_antagonistic_regional_prop"] = region_out_info[region]["host_antagonistic_count"] / num_hosts if num_hosts > 0 else 0
                # Endosymbiont regional proportions
                region_out_info[region]["endosym_occupany_regional_prop"] = region_out_info[region]["endosym_present_count"] / region_size
                region_out_info[region]["endosym_nearly_neutral_regional_prop"] = region_out_info[region]["endosym_nearly_neutral_count"] / num_endosyms if num_endosyms > 0 else 0
                region_out_info[region]["endosym_mutualistic_regional_prop"] = region_out_info[region]["endosym_mutualistic_count"] / num_endosyms if num_endosyms > 0 else 0
                region_out_info[region]["endosym_antagonistic_regional_prop"] = region_out_info[region]["endosym_antagonistic_count"] / num_endosyms if num_endosyms > 0 else 0
                # Host-endosymbiont profile compatibility proportion in this region
                region_out_info[region]["host_sym_task_profile_compatible_regional_prop"] = region_out_info[region]["host_sym_task_profile_compatible_count"] / num_endosyms if num_endosyms > 0 else 0
                region_out_info[region]["host_sym_matching_parent_task_profile_regional_prop"] = region_out_info[region]["host_sym_matching_parent_task_profile_count"] / num_endosyms if num_endosyms > 0 else 0
                # Average host lineage lineage length
                region_out_info[region]["avg_host_lineage_length"] = statistics.mean(region_out_info[region]["avg_host_lineage_length"])

            # Annotate task discovery info with region information
            discovery_types = ["host", "host_parent_offspring", "sym", "sym_parent_offspring"]
            # Initialize regional complex task discovery counts to 0
            for disc_type in discovery_types:
                for region in region_out_info:
                    region_out_info[region][f"{disc_type}_complex_task_discoveries"] = 0
                    region_out_info[region][f"{disc_type}_complexity_gt4_task_discoveries"] = 0

            for task in task_discovery_info:
                task_complexity = utils.get_task_complexity(task)
                for disc_type in discovery_types:
                    discovered = task_discovery_info[task][f"{disc_type}_discovered"] == "1"
                    location = task_discovery_info[task][f"{disc_type}_discovered_location"]
                    region = "NONE" if not discovered else spatial_out_info[location]["region"]
                    task_discovery_info[task][f"{disc_type}_discovered_region"] = region
                    # Count task discoveries with complexity greater than 1
                    if discovered:
                        region_out_info[region][f"{disc_type}_complex_task_discoveries"] += task_complexity > 1
                        region_out_info[region][f"{disc_type}_complexity_gt4_task_discoveries"] += task_complexity >= 4

            # ---------------------------
            # Go ahead and write out spatial out info lines
            locations = [int(loc_id) for loc_id in spatial_out_info]
            locations.sort()
            spatial_out_fields = list(spatial_out_info[str(locations[0])].keys())
            spatial_out_fields.sort()
            write_header = False
            if spatial_data_header is None:
                write_header = True
                spatial_data_header = ",".join(spatial_out_fields)
            elif spatial_data_header != ",".join(spatial_out_fields):
                print("Spatial data header mismatch!")
                exit(-1)
            # Write spatial content line-by-line
            spatial_data_content = []
            for loc in locations:
                loc_id = str(loc)
                spatial_data_content.append(",".join([str(spatial_out_info[loc_id][field]) for field in spatial_out_fields]))
            with open(spatial_data_out_path, "a") as fp:
                if write_header:
                    fp.write(spatial_data_header)
                fp.write("\n")
                fp.write("\n".join(spatial_data_content))
            del spatial_data_content
            del spatial_data
            del spatial_out_info

            # Write out region info lines
            regions = list(region_out_info.keys())
            regions.sort()
            region_fields = list(region_out_info[regions[0]].keys())
            region_fields.sort()
            write_header = False
            region_header = ",".join(region_fields)
            if spatial_region_data_header is None:
                write_header = True
                spatial_region_data_header = region_header
            elif spatial_region_data_header != region_header:
                print("Spatial region data header mismatch!")
                exit(-1)
            # Write out content line-by-line
            spatial_region_content = []
            for region in regions:
                spatial_region_content.append(
                    ",".join([str(region_out_info[region][field]) for field in region_fields])
                )
            with open(spatial_region_data_out_path, "a") as fp:
                if write_header:
                    fp.write(spatial_region_data_header)
                fp.write("\n")
                fp.write("\n".join(spatial_region_content))
            del spatial_region_content
            del region_out_info
            # ---------------------------

        ########################################
        # Add summary info to summary content lines
        summary_content_lines.append(run_summary_info)
        sym_interaction_values_content_lines.append(sym_int_vals_info)

        ############################################################
        # Output time series data and task discovery data for this run
        if run_finished_target:
            # ---- Write time series information ----
            # Order the updates
            time_series_update_order = list(time_series_updates)
            time_series_update_order.sort()
            # Order the fields
            time_series_fields = list(time_series_info[time_series_update_order[0]].keys())
            time_series_fields.sort()
            # If we haven't written the header, write it.
            write_header = False
            if time_series_header is None:
                write_header = True
                time_series_header = ",".join(time_series_fields)
            elif time_series_header != ",".join(time_series_fields):
                print("Time series header mismatch!")
                exit(-1)

            # Write time series content line-by-line
            time_series_content = []
            for u in time_series_update_order:
                time_series_content.append(",".join([str(time_series_info[u][field]) for field in time_series_fields]))
            with open(time_series_fpath, "a") as fp:
                if write_header:
                    fp.write(time_series_header)
                fp.write("\n")
                fp.write("\n".join(time_series_content))
            time_series_content = []
            # ---- ----

            # ----- Write task discovery information ----
            task_disc_order = list(task_discovery_info.keys())
            task_disc_order.sort()
            # Order the fields
            task_disc_fields = list(task_discovery_info[task_disc_order[0]].keys())
            task_disc_fields.sort()
            # If we haven't written the header, write it.
            write_header = False
            if task_discovery_header is None:
                write_header = True
                task_discovery_header = ",".join(task_disc_fields)
            elif task_discovery_header != ",".join(task_disc_fields):
                print("Task discovery header mismatch!")
                exit(-1)

            task_discovery_content = []
            for task in task_disc_order:
                task_discovery_content.append(",".join([str(task_discovery_info[task][field]) for field in task_disc_fields]))
            with open(task_discovery_out_path, "a") as fp:
                if write_header:
                    fp.write(task_discovery_header)
                fp.write("\n")
                fp.write("\n".join(task_discovery_content))
            # ---- ----

    ############################################################
    # Write summary info out
    summary_path = os.path.join(dump_dir, "summary.csv")
    utils.write_csv(summary_path, summary_content_lines)
    sym_int_path = os.path.join(dump_dir, "symbiont_interaction_values.csv")
    utils.write_csv(sym_int_path, sym_interaction_values_content_lines)

    # print incomplete runs
    print("Incomplete runs:")
    print("\n".join(incomplete_runs))

if __name__ == "__main__":
    main()
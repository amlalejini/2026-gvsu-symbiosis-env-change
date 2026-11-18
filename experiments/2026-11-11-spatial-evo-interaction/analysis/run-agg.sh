#!/usr/bin/env bash

EXP_NAME=2026-11-11-spatial-evo-interaction
DATA_DIR=/mnt/scratch/lalejini/data/2026-gvsu-symbiosis-env-change/${EXP_NAME}/
REPO_DIR=/mnt/home/lalejini/devo_ws/2026-gvsu-symbiosis-env-change/
EXP_DIR=${REPO_DIR}/experiments/${EXP_NAME}/
SPATIAL_STRUCT_DIR=${EXP_DIR}/hpc/spatial-structs/

DUMP_DIR=./dump
FINAL_UPDATE=200000
TS_UNITS=interval
TS_RES=10

python3 aggregate.py \
  --data_dir ${DATA_DIR} \
  --dump_dir ${DUMP_DIR} \
  --summary_update ${FINAL_UPDATE} \
  --time_series_units ${TS_UNITS} \
  --time_series_resolution ${TS_RES} \
  --spatial_agg \
  --spatial_struct_dir ${SPATIAL_STRUCT_DIR}
#!/usr/bin/env bash
cd /home/js741/all_pdb_models/progress_logs/
current_date=$(date '+%m_%d_%Y')
make -f /home/js741/all_pdb_models/progress_logs/Makefile clean
make -j4 -f /home/js741/all_pdb_models/progress_logs/Makefile all -k
bash /home/js741/all_pdb_models/progress_logs/report_stats.sh >"/home/js741/all_pdb_models/progress_logs/updates/$current_date"
make -f /home/js741/all_pdb_models/progress_logs/Makefile clean

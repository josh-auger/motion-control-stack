# POST-RUN ANALYSIS SCRIPT
# =====================================
# Bash script to analyze results of motion-control-stack run.
#
# Example run command:
#   sh analyze_results.sh

# DATA=/home/jauger/Radiology_Research/Scan_data/20260924_Phantom_scan/fromPythonFireServer_20260924/savedData_20260924T142618_func-bold_task-rest480_run-01_SLIMMON 
DATA=/home/jauger/GitHubRepos/motion-control-stack/data/savedData_20260921_fetal1_SMS1_TR2500_LIFO_persistent_withProfiling
OUT=analysis_baseline_2

python3 shared/analysis_tools/analyze_queue_profile.py "$DATA" --output-dir "$OUT"
python3 shared/analysis_tools/analyze_motion_dashboard_profile.py "$DATA" --output-dir "$OUT"
python3 shared/analysis_tools/analyze_fire_moco_profile.py "$DATA" --output-dir "$OUT"
python3 shared/analysis_tools/analyze_registration_completeness.py "$DATA" --output-dir "$OUT"
python3 shared/analysis_tools/analyze_registration_calls.py "$DATA" --output-dir "$OUT"
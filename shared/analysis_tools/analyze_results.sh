
DATA=/home/jauger/GitHubRepos/motion-control-stack/data/savedData_20260930T155711_func-bold_task-rest480_run-01_SLIMMON 
OUT=analysis_test_3

python3 shared/analysis_tools/analyze_queue_profile.py "$DATA" --output-dir "$OUT"
python3 shared/analysis_tools/analyze_motion_dashboard_profile.py "$DATA" --output-dir "$OUT"
python3 shared/analysis_tools/analyze_fire_moco_profile.py "$DATA" --output-dir "$OUT"
python3 shared/analysis_tools/analyze_registration_completeness.py "$DATA" --output-dir "$OUT"
python3 shared/analysis_tools/analyze_registration_calls.py "$DATA" --output-dir "$OUT"
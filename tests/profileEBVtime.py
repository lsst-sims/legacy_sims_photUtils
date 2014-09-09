import pstats

p=pstats.Stats('control_profile.sav')
p.strip_dirs().sort_stats('filename').print_stats('EBV.py')

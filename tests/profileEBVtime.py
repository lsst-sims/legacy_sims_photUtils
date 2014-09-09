import pstats

p=pstats.Stats('junk.sav')
p.strip_dirs().sort_stats('filename').print_stats('EBV.py')

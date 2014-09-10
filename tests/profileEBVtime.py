import pstats

print "control\n"
p=pstats.Stats('control_profile.sav')
p.strip_dirs().sort_stats('filename').print_stats('EBV.py')

print "test\n"
pp=pstats.Stats('test_profile.sav')
pp.strip_dirs().sort_stats('filename').print_stats('EBV.py')

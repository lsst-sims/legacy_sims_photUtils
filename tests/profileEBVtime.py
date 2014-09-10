import pstats

print "control\n"
for i in range(10):

    name = 'control_profile_'+str(i)+'.sav'
    p=pstats.Stats(name)
    p.strip_dirs().sort_stats('name').print_stats('get_EBV')

print "test\n"
for i in range(10):
    name='test_profile_'+str(i)+'.sav'
    pp=pstats.Stats(name)
    pp.strip_dirs().sort_stats('name').print_stats('get_EBV')

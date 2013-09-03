# profile the code
import multipath as mp
import hotshot, hotshot.stats 

prof = hotshot.Profile("multipath.prof")
prof.run( 'mp.multipath_moore(0.3, 0.10, "LEIAT504        NONE", 2,1)')

prof.close()

# print the results
stats = hotshot.stats.load("multipath.prof")
stats.strip_dirs()
stats.sort_stats('time','calls')
stats.print_stats(20)


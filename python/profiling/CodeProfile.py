import pstats

p = pstats.Stats('Stats.out')
p.strip_dirs() # remove extraneous path from all the module names

## To see the profile by cumulative time in a function
# p.sort_stats('cumulative') # sort profile by cumulative time in a function

## To see what functions were looping a lot, and taking a lot of time:
p.sort_stats('time')

## Print 15 lines of the output based on the above criteria
p.print_stats(10)

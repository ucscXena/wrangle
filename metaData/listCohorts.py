import string, sys
import xena_query as xena

if len(sys.argv[:]) != 3:
    print "python listCohort.py huburl output_txt"
    sys.exit()

hub = sys.argv[1]
output = sys.argv[2]

allCohorts = xena.all_cohorts(hub)
fout = open(output,'w')
fout.write (hub + '\n')
fout.write('\n')
for cohort in allCohorts:
    fout.write(cohort+'\n')
fout.close()

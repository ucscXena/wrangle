import string, xenaPython, sys

if len(sys.argv[:])!=2:
	print "python hubSampleN.py hub"
	sys.exit()

hub = sys.argv[1]

cohorts = xenaPython.all_cohorts(hub)
total =0
for cohort in cohorts:
	n = len(xenaPython.cohort_samples(hub, cohort, None))
	total = total +n
	print cohort, n
print len(cohorts), total
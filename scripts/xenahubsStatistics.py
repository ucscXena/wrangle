import xena_query as xena

def uniqSamples(huburl):
  allSamples ={}
  allCohorts = xena.all_cohorts(huburl)

  for cohort in allCohorts:
    samples = xena.all_samples(huburl, cohort)
    for sample in samples:
      if sample not in allSamples:
        allSamples[sample]=""
  return allSamples

def totalSamples(hubs):
  totalSamples={}
  for hub in hubs:
    samples = uniqSamples(hub)
    print hub, "samples", len(samples)
    for sample in samples:
      if sample not in totalSamples:
        totalSamples[sample]=""
  print "total samples:", len(totalSamples)


def totalDatasets(hubs):
  totalDatasets =[]
  for hub in hubs:
    totalDatasets.extend(xena.datasets_list(hub))
  print "total datasets:", len(totalDatasets)
  print xena.datasets_list(hub)

hubs = [
#  "https://tcga.xenahubs.net",
#  "https://icgc.xenahubs.net",
#  "https://toil.xenahubs.net",
  "https://ucscpublic.xenahubs.net",
#  "https://icgc.xenahubs.net",
#  "https://gdc.xenahubs.net",
#  "https://treehouse.xenahubs.net"
]
totalSamples(hubs)

totalDatasets(hubs)

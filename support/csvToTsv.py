import pandas
import sys

if len(sys.argv[:])!=3:
    print ("python csvToTsv.py input_csv output_tsv")
    sys.exit()

input_csv = sys.argv[1]
ouptupt_tsv = sys.argv[2]

df = pandas.read_csv(input_csv, index_col=0)
print(df)
df.to_csv(ouptupt_tsv, sep="\t")


import h5py
import string, sys

def print_attrs(name, obj):
    print name, len(obj), obj
    if name[-6:] == "/shape":
        print obj[0],obj[1]

def get_h5_info (h5_file):
    hF = h5py.File(h5_file)
    print hF.keys()
    hF.visititems(print_attrs)
    print "data :5", hF["mm10/data"][:5]
    print "indices :5", hF["mm10/indices"][:5]
    print "indptr :3", hF["mm10/indptr"][0:3]

if __name__ == "__main__" and len(sys.argv[:])!=2:
    print "pyton h5_xena_T.py h5file"
    sys.exit()

h5file = sys.argv[1]
get_h5_info (h5file)


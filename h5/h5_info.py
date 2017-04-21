import h5py
import string, sys

def print_attrs(name, obj):
    print name, len(obj),obj

def get_h5_info (h5_file):
    hF = h5py.File(h5_file)
    print hF.keys()
    hF.visititems(print_attrs)

if __name__ == "__main__" and len(sys.argv[:])!=2:
    print "pyton h5_xena_T.py h5file"
    sys.exit()

h5file = sys.argv[1]
get_h5_info (h5file)

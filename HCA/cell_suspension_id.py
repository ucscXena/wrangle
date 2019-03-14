import sys

def cellSuspensionID (cell_suspension_file):
	# get the cell suspension id
	fin = open(cell_suspension_file, 'U')
	J = json.loads(fin.read())
	fin.close()
	cell_suspension_id = J["provenance"]["document_id"]

	return cell_suspension_id

if __name__ == "__main__" and len(sys.argv[:]) != 2:
	print "pyton cell_suspension_id.py cell_suspension_json_file\n"
	sys.exit()

cell_suspension_json_file = sys.argv[1]
cell_suspension_id = cellSuspensionID(cell_suspension_json_file)
print cell_suspension_id

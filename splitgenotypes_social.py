#!/usr/bin/python
import glob,os,re,shutil

commongenos = ["het","hetandwt","hom","wt","wtwt","homhom","hethet","hethom","homhet","hetwt","wthet","wtandhet","homwt","wthom","mut","drug","dmso"]

wellsections = {'1':['A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'A7', 'A8', 'A9', 'A10', 'A11', 'A12', 'B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8'], '2':['B9', 'B10', 'B11', 'B12', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9', 'C10', 'C11', 'C12', 'D1', 'D2', 'D3', 'D4'], '3':['D5', 'D6', 'D7', 'D8', 'D9', 'D10', 'D11', 'D12', 'E1', 'E2', 'E3', 'E4', 'E5', 'E6', 'E7', 'E8', 'E9', 'E10', 'E11', 'E12'], '4':['F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8', 'F9', 'F10', 'F11', 'F12', 'G1', 'G2', 'G3', 'G4', 'G5', 'G6', 'G7', 'G8'], '5':['G9', 'G10', 'G11', 'G12', 'H1', 'H2', 'H3', 'H4', 'H5', 'H6', 'H7', 'H8', 'H9', 'H10', 'H11', 'H12','I1','I2','I3','I4'], '6':['I5','I6','I7','I8','I9','I10','I11','I12','J1','J2','J3','J4','J5','J6','J7','J8','J9','J10','J11','J12'], '7':['K1','K2','K3','K4','K5','K6','K7','K8','K9','K10','K11','K12','L1','L2','L3','L4','L5','L6','L7','L8']}

#welltonum = {'A1':20,'A2':19,'A3':18,'A4':17,'A5':16,'A6':15,'A7':14,'A8':13,'A9':12,'A10':11,'A11':10,'A12':9,'B1':8,'B2':7,'B3':6,'B4':5,'B5':4,'B6':3,'B7':2,'B8':1,
#		'B9':20,'B10':19,'B11':18,'B12':17,'B'}
#welltonum = {}
#rows = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
#for row in rows:
#    for col in range(1, 13):
#        well = f"{row}{col}"
#        if row == 'B' and col <= 8:
#            value = 20 - col
#        else:
#            value = 20 - (col - 8)
#        welltonum[well] = value

welltonum = {
    'A1': 20, 'A2': 19, 'A3': 18, 'A4': 17, 'A5': 16, 'A6': 15, 'A7': 14, 'A8': 13, 'A9': 12, 'A10': 11, 'A11': 10, 'A12': 9,
    'B1': 8, 'B2': 7, 'B3': 6, 'B4': 5, 'B5': 4, 'B6': 3, 'B7': 2, 'B8': 1, 'B9': 20, 'B10': 19, 'B11': 18, 'B12': 17,
    'C1': 16, 'C2': 15, 'C3': 14, 'C4': 13, 'C5': 12, 'C6': 11, 'C7': 10, 'C8': 9, 'C9': 8, 'C10': 7, 'C11': 6, 'C12': 5,
    'D1': 4, 'D2': 3, 'D3': 2, 'D4': 1, 'D5': 20, 'D6': 19, 'D7': 18, 'D8': 17, 'D9': 16, 'D10': 15, 'D11': 14, 'D12': 13,
    'E1': 12, 'E2': 11, 'E3': 10, 'E4': 9, 'E5': 8, 'E6': 7, 'E7': 6, 'E8': 5, 'E9': 4, 'E10': 3, 'E11': 2, 'E12': 1,
    'F1': 20, 'F2': 19, 'F3': 18, 'F4': 17, 'F5': 16, 'F6': 15, 'F7': 14, 'F8': 13, 'F9': 12, 'F10': 11, 'F11': 10, 'F12': 9,
    'G1': 8, 'G2': 7, 'G3': 6, 'G4': 5, 'G5': 4, 'G6': 3, 'G7': 2, 'G8': 1, 'G9': 20, 'G10': 19, 'G11': 18, 'G12': 17,
    'H1': 16, 'H2': 15, 'H3': 14, 'H4': 13, 'H5': 12, 'H6': 11, 'H7': 10, 'H8': 9, 'H9': 8, 'H10': 7, 'H11': 6, 'H12': 5,
    'I1': 4, 'I2': 3, 'I3': 2, 'I4': 1, 'I5': 20, 'I6': 19, 'I7': 18, 'I8': 17, 'I9': 16, 'I10': 15, 'I11': 14, 'I12': 13,
    'J1': 12, 'J2':11, 'J3':10, 'J4':9, 'J5':8, 'J6':7, 'J7':6, 'J8':5, 'J9':4, 'J10':3, 'J11':2, 'J12':1,
    'K1': 20, 'K2': 19, 'K3': 18, 'K4': 17, 'K5': 16, 'K6': 15, 'K7': 14, 'K8': 13, 'K9': 12, 'K10': 11, 'K11': 10, 'K12': 9,
    'L1': 8, 'L2': 7, 'L3': 6, 'L4': 5, 'L5': 4, 'L6': 3, 'L7': 2, 'L8': 1, 'L9': 20, 'L10': 19, 'L11': 18, 'L12': 17
}


alphatonum = {'A': 0, 'B': 12, 'C': 24, 'D': 36, 'E': 48, 'F': 60, 'G': 72, 'H': 84, 'I':96, 'J':108, 'K':120, 'L':132, 'M': 144, 'N':156, 'O':168, 'P':180, 'Q':192}
alphatonum_col = {'A': 0, 'B': 1, 'C': 2, 'D': 3, 'E': 4, 'F': 5, 'G': 6, 'H': 7, 'I':8, 'J':9, 'K':10, 'L':11, 'M': 12, 'N':13, 'O':14, 'P':15, 'Q':16}
alphatonum_colflip = {'A': 7, 'B': 6, 'C': 5, 'D': 4, 'E': 3, 'F': 2, 'G': 1, 'H': 0}

def Reverse(lst):
	return [ele for ele in reversed(lst)]

def zipper(a,b):
	list = [a[i] + b[i] for i in range(len(a))]
	return list

current_directory = os.getcwd()
directory_names = [item for item in os.listdir(current_directory) if os.path.isdir(os.path.join(current_directory, item))and ("clear" in item or "white" in item)]
numberlist = []
for name in directory_names:
	#print(name)
	match = re.search(r'\d+(?=(clear|white))', name)
	if match:
		number = match.group(0)
		#print(number)
		numberlist.append(str(number))
numberlist = list(set(numberlist))
for n in numberlist:
	gfile = open("genotyping", 'w')
	for file in glob.glob('*matrix'):
		f = open(file, 'r')
		line1 = f.readline()
		lines = f.readlines()
		wells = {}
		for line in lines:
			if len(line.split()) != 13:
				print( "Either row not correct length or duplicated gene set, length = ", len(line.split()))
				if len(line.split()) != 26:
					print( "REAL ERROR REAL ERROR REAL ERROR, length = ", len(line.split()))
				wells[line.split()[0]] = zipper(line.split()[1:13], line.split()[14:26])
				# IF YOUR CAMERA IS FLIPPED AROUND, UNCOMMENT BELOW LINE AND COMMENT ABOVE LINE
				#wells[line.split()[0]] = Reverse(zipper(line.split()[1:13], line.split()[14:26]))
			else:
				wells[line.split()[0]] = line.split()[1:13]
				# IF YOUR CAMERA IS FLIPPED AROUND, UNCOMMENT BELOW LINE AND COMMENT ABOVE LINE
				#wells[line.split()[0]] = Reverse(line.split()[1:13])
		valuelist = []
		for v in wells.values():
			valuelist = valuelist + v
		valueset = set(valuelist)
		finaldata = {}
		for v2 in valueset:
			finaldata[v2] = []
		#print(wells)
		for a in wells.keys():
			for v3 in range(0, len(wells[a])):
				#num = ((v3)*8)+1 + alphatonum_col[a]
				# IF YOUR CAMERA IS FLIPPED AROUND, UNCOMMENT BELOW LINE AND COMMENT ABOVE LINE
				well = str(a+str(v3+1))
				#print(a, v3, v3+1, well, welltonum[well])
				#num = ((v3)*8)+1 + alphatonum_colflip[a]
				if well in wellsections[n]:
					num = welltonum[well]
					finaldata[wells[a][v3]].append(num)
		for k in finaldata.keys():
			finaldata[k].sort()
		if "het" in finaldata and "wt" in finaldata:
			if finaldata["het"] and finaldata["wt"]:
				finaldata["hetandwt"] = finaldata["het"] + finaldata["wt"]
		#print(file)
		#print(finaldata)
		gfile.write(file)
		gfile.write('\n')
		for key in finaldata.keys():
			title = file.split('_')[0] + "_" + key
			if key in commongenos:
				title = "*" + title
			gfile.write(title + ": ")
			gfile.write(str(finaldata[key]).strip().strip("[").strip("]").replace(" ", ""))
			gfile.write('\n')
		gfile.close()
		for name in directory_names:
			match = re.search(r'\d+(?=(clear|white))', name)
			if match:
				number = int(match.group(0))
				if number == int(n):
					#print(name)
					destination_dir = os.path.join(name, "genotyping")
					#print(n,destination_dir)
					shutil.copy("genotyping", destination_dir)

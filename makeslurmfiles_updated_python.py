#!/usr/bin/python3

import os,sys,glob,argparse,shutil

# Input arguments
parser = argparse.ArgumentParser(description='options for slurm file generation')
parser.add_argument('-module', type=str, action="store", dest="module", default="Anaconda3")
parser.add_argument('-partition', type=str, action="store", dest="partition", default="short")
parser.add_argument('-statsfile', type=str, action="store", dest="statsfile", default="linearmodel") # prefix for this file
parser.add_argument('-gfile', type=str, action="store", dest="gfile", default="genotyping")
parser.add_argument('-efile', type=str, action="store", dest="efile", default="../scheduled-events")
parser.add_argument('-prefix', type=str, action="store", dest="prefix", default="../pre-processed.csv") # slow-speed data prefix
parser.add_argument('-hprefix', type=str, action="store", dest="hprefix", default="../hsmovie") # high-speed data prefix
parser.add_argument('-pfile', type=str, action="store", dest="pfile", default="../PlotParameters")
parser.add_argument('-esfile', type=str, action="store", dest="esfile", default="../sectionsfile")
parser.add_argument('-rois', type=str, action="store", dest="rois", default="../zebrafish-tracker.cells")
parser.add_argument('-scriptpath', type=str, action="store", dest="scriptpath", default="/pi/summer.thyme-umw/behaviordata/behavior-analysis-scripts/") # path for the analysis suite
parser.add_argument('-other', type=str, action="store", dest="other", default="") # args to add at the end. YOU CANNOT INCLUDE THE INITIAL DASH FOR A NEW ARG WITH ARGPARSER, WILL BE ADDED BY THIS SCRIPT"

args = parser.parse_args()
module = args.module
partition = args.partition
statsfile = args.statsfile
gfile = args.gfile
efile = args.efile
prefix = args.prefix
hprefix = args.hprefix
pfile = args.pfile
esfile = args.esfile
rois = args.rois
scriptpath = args.scriptpath
other = args.other

# This logic determins what is a control and what is the "mutant" comparison. Add identifiers as needed.
logic = {				 "wt":0,
					"wtwt":0,
					"hetwt":2,
					"wthet":1,
					"wtandhet":1,
					"hetandwt":1,
					"hethet":3,
					"het":3,
					"dmso-het":3,
					"dmso-hetandwt":2,
					"ptz-hetandwt":5,
					"drug-hetandwt":5,
					"dmso-hom":8,
					"homwt":5,
					"wthom":4,
					"hethom":6,
					"homhet":7,
					"hom":7,
					"mut":7,
					"homhom":8,
					"transgenic":8,
					"drug":1,
					"dmso":0,
					"dmso-wt":1,
					"ptz-hom":10,
					"drug-hom":10,
					"ptz-wt":4,
					"drug-wt":4,
					"drug-het":6,
					"ptz-het":6
				}


def make_slurm_file(fd1, fd2, hfile):
	ffile = open('submission_script_' + fd1 + "_vs_" + fd2 + '.slurm', 'w')
	ffile.write("#!/bin/bash\n")
	ffile.write("#BSUB -q short # Partition to submit to\n")
	ffile.write("#BSUB -n 1 # Number of cores requested\n")
	ffile.write("#BSUB -R \"span[hosts=1]\"\n")
	ffile.write("#BSUB -W 400 # Runtime in minutes\n")
	ffile.write("#BSUB -R \"rusage[mem=64000]\"\n")
	ffile.write("#BSUB -e hostname_%J_%I.err # Standard err goes to this filehostname\n")
	ffile.write("#BSUB -o " + statsfile + "_" + fd1 + "_vs_" + fd2 + "_%J_%I.out # Standard out goes to this file\n")
	#ffile.write("module load " + module + "\n")
	ffile.write("cd outputfulldata_" + fd1 + "_vs_" + fd2 + '\n')
	if os.path.exists("outputfulldata_" + fd1 + "_vs_" + fd2):
		print("Output directories already generated.")
	else:
		os.mkdir("outputfulldata_" + fd1 + "_vs_" + fd2)
	ffile.write("~/miniconda3/bin/python3 " + scriptpath + "processmotiondata.py -e ")
	ffile.write("\"" + efile + "\" -c \"" + prefix)
	ffile.write("\" -m \"../")
	ffile.write("\" -g \"../" + hfile)
	ffile.write("\" -s \"" + esfile)
	ffile.write("\" -j \"" + pfile)
	ffile.write("\" -r \"" + rois)
	ffile.write("\" -python)
	ffile.write("\" -output_folder \"." + "\"")
	if other != "":
		ffile.write(" -" + other)
	return 'submission_script_' + fd1 + "_vs_" + fd2 + '.slurm'

#~/miniconda3/bin/python3 /pi/summer.thyme-umw/behaviordata/behavior-analysis-scripts/processmotiondata.py " -e "../scheduled-events" -c "../pre-processed.csv" -m "../" -g "../python_wt_vs_python_hom_scripted_inputgenotypeids" -s "../sectionsfile" -j "../PlotParameters"" -r "../zebrafish-tracker.cells"" -output_folder "."

#~/miniconda3/bin/python3 /pi/summer.thyme-umw/behaviordata/behavior-analysis-scripts/processmotiondata.py -t "../testlog.timestamp1.Tue, Dec 23, 2025" -e "../fulltestrun_final_01_27_2020" -c "../testlog.centroid1.Tue, Dec 23, 2025" -d "../testlog.motion1.Tue, Dec 23, 2025" -m "../" -g "../labview_wt_vs_labview_hom_scripted_inputgenotypeids" -r "../rois_string" -s "../sectionsfile" -j "../PlotParameters" -output_folder "."

genofile = open(gfile, 'r')
lines = genofile.readlines()
dirnames = []
iddict = {}
for line in lines:
	if line.startswith('*'):
		destdir = line.strip().split(':')[0][1:]
		itype = destdir.split('_')[1]
		dirnames.append(destdir)
		ids = line.strip().split(':')[1].strip()
		iddict[itype] = ids
sfile = open("jobsubmission.sh", 'w')
sfile.write("#!/bin/bash\n")
for fd1 in dirnames:
	for fd2 in dirnames:
		if fd1.split('_')[0] == fd2.split('_')[0]: # Check if same gene
			if fd1 != fd2:
				testfile = 'submission_script_' + fd1 + "_vs_" + fd2 + '.slurm'
				testfileo = 'submission_script_' + fd2 + "_vs_" + fd1 + '.slurm'
				if os.path.exists(testfile) or os.path.exists(testfileo):
					continue
				type1 = fd1.split('_')[1]
				type2 = fd2.split('_')[1]
				if (type1 == "het") and (type2 == "hetandwt"):
					continue
				if (type1 == "wt") and (type2 == "hetandwt"):
					continue
				if (type1 == "hetandwt") and (type2 == "het"):
					continue
				if (type1 == "hetandwt") and (type2 == "wt"):
					continue
				if logic[type1] > logic[type2]:
					hfile = fd2 + "_vs_" + fd1 + "_scripted_inputgenotypeids"
					ifile = open(hfile, 'w')
					ifile.write("controlgroup_" + type2 + ":" + iddict[type2])
					ifile.write("\ntestgroup_" + type1 + ":" + iddict[type1] + "\n")
					fname2r = make_slurm_file(fd2, fd1, hfile)
				else:
					hfile = fd1 + "_vs_" + fd2 + "_scripted_inputgenotypeids"
					ifile = open(hfile, 'w')
					ifile.write("controlgroup_" + type1 + ":" + iddict[type1])
					ifile.write("\ntestgroup_" + type2 + ":" + iddict[type2] + "\n")
					fname2r = make_slurm_file(fd1, fd2, hfile)
				sfile.write("bsub <")
				sfile.write(fname2r)
				sfile.write("\nsleep 1\n")

os.system("chmod +x jobsubmission.sh")

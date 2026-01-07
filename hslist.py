import subprocess
import os
import sys

def find_all_videos_for_tracking(path=None, ext="avi"):
    """Finds all avi files in current working directory, if path given finds all
    avi files in path that aren't output of the algorithm."""
    if path is None:
        path = os.getcwd()

    files_to_read = []
    for ex in os.walk(path):
        for filename in ex[2]:
            if filename[-3:] == ext and "long" not in filename:
                files_to_read.append(ex[0] + "/" + filename)

    return files_to_read

if __name__ == "__main__":
    path = sys.argv[1]
    files = find_all_videos_for_tracking(path)
    files.sort(key=lambda name: int(name.split('/')[-1].split('.')[0].split('-')[1]))
    print(files)

    for file in files:
#   subprocess.run(["python", "highspeedmovieanalysis.py", "-r", "../../labview_comp_06_03_25/zebrafish-tracker-6-3-25.cells", f"-m {file}", "-e", "../../labview_comp_06_03_25/scheduled-events", "-fd", "660, 992", "-f", "274"])
#   print(["python", "highspeedmovieanalysis.py", "-r", "../../labview_comp_06_03_25/zebrafish-tracker-6-3-25.cells", f"-m {file}", "-e", "../../labview_comp_06_03_25/scheduled-events", "-fd", "660, 992", "-f", "274"])
        try:
            filenum = int(file.split('/')[-1].split('.')[0].split('-')[1])

            cmd = f'python highspeedmovieanalysis.py -r "{path}/zebrafish-tracker.cells" -m "{file}" -e "{path}/scheduled-events" -fd "660, 992" -f 285 -filenumber {filenum}'
            print("COMMAND:", cmd)
            os.system(cmd)
        except Exception as e:
            print(e)


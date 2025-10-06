import subprocess
import os

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

files = find_all_videos_for_tracking("../../labview_comp_06_24_25")

for file in files:
#   subprocess.run(["python", "highspeedmovieanalysis.py", "-r", "../../labview_comp_06_03_25/zebrafish-tracker-6-3-25.cells", f"-m {file}", "-e", "../../labview_comp_06_03_25/scheduled-events", "-fd", "660, 992", "-f", "274"])
#   print(["python", "highspeedmovieanalysis.py", "-r", "../../labview_comp_06_03_25/zebrafish-tracker-6-3-25.cells", f"-m {file}", "-e", "../../labview_comp_06_03_25/scheduled-events", "-fd", "660, 992", "-f", "274"])
    try:
        cmd = f'python highspeedmovieanalysis.py -r "../../labview_comp_06_24_25/zebrafish-tracker-6-24-25-realrun.cells" -m "{file}" -e "../../labview_comp_06_24_25/scheduled-events" -fd "660, 992" -f 285'
        print("COMMAND:", cmd)
        os.system(cmd)
    except Exception as e:
        print(e)

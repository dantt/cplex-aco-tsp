import os
import subprocess


for filename in os.listdir('instances'):
    print filename
    subprocess.call("./cplex/main" + " ./instances/" + filename, shell=True)
    subprocess.call("./heur/main" + " ./instances/" + filename, shell=True)
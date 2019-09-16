import subprocess
import sys
import os

thetas = [0.5,1.5,2.5,3.5,4.5,5.5,6.5]
for i in range(7):
    subprocess.call(['python', 'runplotSED.py','7','1', '100', '30', str(thetas[i])])
    print(i);


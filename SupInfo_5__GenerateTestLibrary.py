import random
import math
### Generates a list of example Ka values based upon values
### reported by Clark etal

def readFile(name):
    '''takes name and reads the file row by row to list.  Returns list''' 
    content = 0
    with open(name) as f:
        content = f.readlines()
    f.close()
    return content


input_kvalues = readFile('Kvalues_FromPaper.txt')
log_values = []
for x in input_kvalues:
    log_values.append(math.log10(float(x)))
    


count = 0
while count < 707788-len(input_kvalues):
    value = abs(random.gauss(0,2))
    if value < 5.9:
        log_values.append(value)
    count += 1
log_values = sorted(log_values, reverse = True)

file1 = open('example_kavalues.txt', 'w')

for x in log_values:
    file1.write(str(x) + "\n")

    
file1.close()    


    
##5.789

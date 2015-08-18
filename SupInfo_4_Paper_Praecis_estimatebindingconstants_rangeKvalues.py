import random
import math

####################################################################################################
###  Sampling of an eqilibrium association constant calculated from the number of times                                  
###  read by highthroughput seqeuncing.  Random variables include synthetic yeild of the library,
###  total amount of recovered DNA, and concentration of active protein.
###################################################################################################

read_count = 205                    ## sequence reads for a particular library member(ligand)
total_reads = 190363                ## total number of usable sequence reads
diversity = 7077888.0               ## total number of distinct molecules in the library


def getConc(content, total_reads, del_conc_final):
    '''takes list of sequence counts and returns list of
       final concentrations for each library member'''
    list1 = []
    for n in content:
        counts = float(n)
        concitem = (counts/total_reads) * del_conc_final
        list1.append(concitem)
    return list1

def getKa(B, A, list_conc_final):
    '''takes floats A, B, and list of library member final concentrations
       and returns list of calculated association equilibrium constants'''
    list1 = []
    for i in list_conc_final:
        a = (B**3)*i - (B**3)*A
        b = 3*(B**2)*i
        c = 3*B*i
        d = i
        delta = 18*a*b*c*d - 4*(b**3)*d + (b**2)*c**2 - 4*a*c**3 - 27*(a**2)*(d**2)
        delta0 = b*b - 3*a*c
        delta1 = 2*(b**3) - 9*a*b*c + 27*(a**2)*d
        delta2 = -27*(a**2)*delta
        u1 = 1
        u2 = (-1 + 1.0j*(3**(1.0/2)))/2
        u3 = (-1 - 1.0j*(3**(1.0/2)))/2
        C = ((delta1 + (delta2**(1.0/2)))/2)**(1.0/3)
        K1 = (-1/(3*a))*(b + u1*C + delta0/(u1*C))
        K2 = (-1/(3*a))*(b + u2*C + delta0/(u2*C))
        K3 = (-1/(3*a))*(b + u3*C + delta0/(u3*C))
        list1.append(K1)
    return list1

def outputFile(list1, name):
    '''takes list and str filename and writes list to file'''
    file1 = open(name, 'w')
    for x in list1:
        file1.write(str(x) + "\n")
    file1.close()

def readFile(name):
    '''takes name and reads the file row by row to list.  Returns list''' 
    content = 0
    with open(name) as f:
        content = f.readlines()
    f.close()
    return content
    

list_Ka_range = []
count = 0
while count < 100000:

    del_conc_initial = random.gauss(1.7 * 10**(-5), 0.5*10**(-5))       ## initial concentration of DEL.                                                                      ## for the final product
    protein_conc = random.gauss(0.16*10**(-6), 0.05*10**(-6))           ## concentration of active protein
    del_conc_final = random.gauss(33.0* (10**(-12)), 30.0*(10**(-12)))  ## see section S3
    if protein_conc > 0 and del_conc_initial > 0 and del_conc_final > 0:
        conc_ligand_initial = del_conc_initial/diversity                ## starting concentration of each ligand
        list_conc_final = getConc([read_count], total_reads, del_conc_final)
        list_Ka = getKa(protein_conc, conc_ligand_initial, list_conc_final)
        if list_Ka[0] > 0:
            list_Ka_range.append(list_Ka[0])
            count += 1

outputFile(list_Ka_range, 'sampling_Kvalues.txt')
    

    


    
    
    
    


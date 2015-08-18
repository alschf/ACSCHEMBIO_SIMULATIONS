import random
import math

####################################################################################################
###  This script takes a list of 'counts' from sequencing data and calculates theoretical       ###
###  equilibrium associaton constants for each library memeber.  Assumes 3 rds of Selection.    ###
###################################################################################################




### required values  
total_reads = 190363                ## total number of usable sequence reads
diversity = 7077888.0               ## total number of distinct molecules in the library
del_conc_initial = 1.7 * 10**(-5)   ## initial concentration of DEL.  Here we estimate an average yeild of 33%
                                    ## for the final product. In 100ul

protein_conc = float(.16*10**(-6))  ## We assume 80% of the protein is active 
del_conc_final = 24* (10**(-12))    ## Estimated assuming a DNA concentration of 440 nM after 20 rounds of PCR
                                    ## and that 1/2 of each sample was set aside following each round of selection.

conc_ligand_initial = del_conc_initial/diversity ## concentration of each library member prior to selection   


#########################################

def getConc(content, total_reads, del_conc_final):
    '''takes list of sequence counts and returns list of
       final concentrations for each library member'''
    print type(content)
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
        delta = 18*a*b*c*d - 4*(b**3)*d + (b**2)*c*c - 4*a*c**3 - 27*(a**2)*(d**2)
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
    

content = readFile('listofrawcounts.txt')

list_conc_final = getConc(content, total_reads, del_conc_final)

list_Ka = getKa(protein_conc, conc_ligand_initial, list_conc_final)

outputFile(list_Ka, 'Kvalues.txt')
    

    


    
    
    
    


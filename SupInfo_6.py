import random
import math

##############################################################
## The values below can be adjusted to mimic different selection conditions
volume = 100.0*10**(-6)         ## total volume of the selection 
diversity = 707788.0            ## total # of unique sequences.  
del_conc_initial = 50*10**(-6)   ## conc of DEL during incubation with protein
protein_conc = 500*10**(-9)     ## conc of protein used during a selection 
reads = 190363                  ## sequence reads
mean_yield = 45                 ## avg yield individual library member
yield_sd  = 0                  ## standard deviation of mean_yield
                               
##########################################################
ligand_conc = del_conc_initial/diversity                 
umolesbegin = float((del_conc_initial*volume)*1000000)  


def readFile(name):
    '''takes name and reads the file row by row to list.  Returns list''' 
    content = 0
    with open(name) as f:
        content = f.readlines()
    f.close()
    return content

def assignStartConc(listK, mean, sd, ligand_conc):
    '''Takes a list of Ka values, and 3 floats (mean, sd, ligand_conc).
       Generates a random 'yield' form a norml distribution
       and assigns a starting ligand concentration for each library member.
       Returns a list of concentrations (float).
       Also returns total micromoles of DNA before selection.''' 
    total = 0.0
    list_concentrations = []
    for x in range(0, len(list_gen_Kavalues)):
        concAB = random.gauss(mean, sd) * 0.01 * ligand_conc
        if concAB < 0:
            concAB = 0
        list_concentrations.append(concAB)  
        total = total + concAB
    umolesend = (total*volume)*1000000.0
    return list_concentrations, umolesend

def simulate_n_RoundsSelection(list_concentrations, list_gen_Kavalues, n):
    '''takes two lists and an integer.  Simulates n rounds of selections
       Returns list of concentrations for each library member (floats).
       Also returns total micromoles of DNA after selection is complete.''' 
    for x in range(0, n):
        total = 0.0
        for i in range(0, len(list_gen_Kavalues)):
            ka = 10.0**(float(list_gen_Kavalues[i]))
            list_concentrations[i] = (ka*protein_conc*list_concentrations[i]/(1+ka*protein_conc))
            total = total + list_concentrations[i]*0.5  
        umolesend = (total*volume)*1000000.0
    return list_concentrations, umolesend

def calcEnrichment(list_concentrations,ligand_conc, del_conc_initial):
    '''takes list of floats and returns list of floats.  Converts concentrations
       to enrichment ratio.''' 
    list_enrichment = []
    for x in range(0, len(list_concentrations)):
        list_enrichment.append((list_concentrations[x]/ligand_conc) /(ligand_conc/del_conc_initial))
    return list_enrichment

def convertEnrichmentToInteger(list_enrichment, mult_factor):
    '''takes list of floats and a float.  Enrichment values are adjusted to prevent
       int boxes_total from exceeding 10**8.  Returns list and int''' 
    boxes_total = 10**9
    list_boxsize = []
    lessthan1 = 0.0

    while boxes_total > 10**8:
        list_boxsize = []
        lessthan1 = 0.0
        mult_factor *= 2
        boxes_total = 0
        for x in list_enrichment:
            x = float(x/mult_factor) 
            if x >= 1:
                boxes_total += int(x)
                list_boxsize.append(int(x))  
                lessthan1 += x - int(x)
            if x < 1:
                list_boxsize.append(0)
                lessthan1 += x

        list_boxsize.append(int(lessthan1))
        boxes_total += int(lessthan1)
    return list_boxsize, boxes_total



def randomChoiceBoxes(boxes_total):
    '''takes int.  Returns list with length equal to boxes_total.
       Indices of the list are populated randomly to simulate
       sequencing experiment.''' 
    list_pickboxes = []
    for x in range(0, boxes_total):
        list_pickboxes.append(0)

    for x in range(0, reads):
        r = random.randrange(0, boxes_total)
        list_pickboxes[r] = list_pickboxes[r] + 1
    return list_pickboxes

def assignRandomChoicesToLigands(list_boxsize, list_pickboxes, mult_factor):
    '''takes two lists and an float.  Indices populated in list_pickboxes are
       assigned to particular library members according to the information stored
       in list_boxsize.  Counts and enrichment for each library member is stored
       in listofcount and list_enrich.  Both lists are returned.  '''
    total = 0
    item = 0
    listofcounts = []
    list_enrich = []
    for x in list_boxsize:
        total = total + x
        listofcounts.append(0)
        for i in range (total-x, total):
            listofcounts[item] = listofcounts[item] + list_pickboxes[i]
        list_enrich.append(x* mult_factor )     ## enrichment
        item +=1
    return list_enrich, listofcounts

def printToFile(numberlines, list_enrich, listofcounts, list_gen_Kavalues):
    '''takes an int and 3 lists.  Writes to a file information from the 3 lists.
       Writes the first number of lines dictated by the int numberline'''
    if numberlines > len(list_enrich)-1:
        numberlines = len(list_enrich)-1
    file1 = open('counts_ryields_july15_withend.txt', 'w')    
    for x in range(0,numberlines-1):
        file1.write(str(list_enrich[x]) + " " + str(listofcounts[x]) + " " + str(x+1) + " " + str(list_gen_Kavalues[x]) )
    ## the line below may be used to observe the counts due to all the <1 enrichments
    ##file1.write(str(list_enrich[-1]) + " " + str(listofcounts[-1]))
    file1.close()


list_gen_Kavalues = readFile("example_kavalues_new.txt")
list_concentrations, umolesend = assignStartConc(list_gen_Kavalues, mean_yield, yield_sd, ligand_conc)
list_concentrations, umolesend = simulate_n_RoundsSelection(list_concentrations, list_gen_Kavalues, 3)
list_enrichment_adjusted = calcEnrichment(list_concentrations,ligand_conc, del_conc_initial)

##Sequencing
mult_factor = 0.5
list_boxsize, boxes_total = convertEnrichmentToInteger(list_enrichment_adjusted, mult_factor)
list_pickboxes = randomChoiceBoxes(boxes_total)
list_enrich, listofcounts = assignRandomChoicesToLigands(list_boxsize, list_pickboxes, mult_factor)
printToFile(40000, list_enrich, listofcounts, list_gen_Kavalues)

        
        

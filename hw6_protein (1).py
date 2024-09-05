"""
15-110 Hw6 - Protein Sequencing Project
Name: Kotnana Venkatesh Naidu
AndrewID:2023501096
"""

import hw6_protein_tests as test

project = "Protein" # don't edit this

### WEEK 1 ###

'''
readFile(filename)
#1 [Check6-1]
Parameters: str
Returns: str
'''
def readFile(filename):
    with open(filename,'r') as file:
        content=file.read()
        contents=content.replace("\n","")
    return contents



'''
dnaToRna(dna, startIndex)
#2 [Check6-1]
Parameters: str ; int
Returns: list of strs
'''
def dnaToRna(dna, startIndex):
    list = []
    for i in range(startIndex,len(dna),3):
        codon = dna[i:i+3]
        if "T" in codon:
            codon = codon.replace("T","U")
        list.append(str(codon))
        if codon in ["UAA","UAG","UGA"] :
            break
    return list
'''
makeCodonDictionary(filename)
#3 [Check6-1]
Parameters: str
Returns: dict mapping strs to strs
'''
def makeCodonDictionary(filename):
    import json
    with open(filename,'r') as file:
        calling_json = json.load(file)
        load_dict = {}
        for amino,codons in calling_json.items():
            for codon in codons:
                if 'T' in codon:
                    codon= codon.replace('T','U')
                load_dict[codon] = amino
    #print(load_dict["AAA"])
    return load_dict


'''
generateProtein(codons, codonD)
#4 [Check6-1]
Parameters: list of strs ; dict mapping strs to strs
Returns: list of strs
'''
def generateProtein(codons, codonD):
    list = []
    for i in range(len(codons)):
        if i == 0 and codons[i] == 'AUG':
            list.append('Start')
        elif codons[i] in ["UAA","UAG","UGA"]:
            list.append('Stop')
            break
        else:
            list.append(codonD[codons[i]])

    return list


'''
synthesizeProteins(dnaFilename, codonFilename)
#5 [Check6-1]
Parameters: str ; str
Returns: 2D list of strs
'''
def synthesizeProteins(dnaFilename, codonFilename):
        dna = readFile(dnaFilename)
        codonD = makeCodonDictionary(codonFilename)
        proteins = []
        i = 0
        count = 0
        while i < len(dna):
            if dna[i:i+3] == 'ATG':
                rna = dnaToRna(dna,i)
                protein = generateProtein(rna,codonD)
                proteins.append(protein)
                i +=3 * len(rna)
            else:
                i += 1
                count +=1
        print(f"the total number of bases : {len(dna)}")
        print(f"the number of unused base : {count}")
        print(f'the total number of proteins synthesized : {len(proteins)}')
        return proteins
        


def runWeek1():
    print("Human DNA")
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    print("Elephant DNA")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")


### WEEK 2 ###

'''
commonProteins(proteinList1, proteinList2)
#1 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs
Returns: 2D list of strs
'''
def commonProteins(proteinList1, proteinList2):
    protein_list = []
    for i in proteinList1:
        for j in proteinList2:
            if i == j and i not in protein_list:
                protein_list.append(i)

    return protein_list


'''
combineProteins(proteinList)
#2 [Check6-2]
Parameters: 2D list of strs
Returns: list of strs
'''
def combineProteins(proteinList):
    combine_list = []
    for i in proteinList:
        for j in i:
            combine_list.append(j)
 
    return combine_list


'''
aminoAcidDictionary(aaList)
#3 [Check6-2]
Parameters: list of strs
Returns: dict mapping strs to ints
'''
def aminoAcidDictionary(aaList):
    dic  = {}
    for i in aaList:
            if i in dic:
                dic[i] += 1
            else:
                dic[i] = 1

    return dic


'''
findAminoAcidDifferences(proteinList1, proteinList2, cutoff)
#4 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs ; float
Returns: 2D list of values
'''
def findAminoAcidDifferences(proteinList1, proteinList2, cutoff):
    list = []
    combine_list1 = combineProteins(proteinList1)
    combine_list2 = combineProteins(proteinList2)
    protein_list1 = aminoAcidDictionary(combine_list1)
    protein_list2 = aminoAcidDictionary(combine_list2)
    for i in protein_list1:
        if i not in protein_list2:
            protein_list2[i] = 0

    for j in protein_list2:
        if j not in protein_list1:
            protein_list1[j] = 0

    for z in protein_list1:
        if z not in ['Start','Stop']:
            probs1 = protein_list1[z]/len(combine_list1)
            probs2 = protein_list2[z]/len(combine_list2)

            if abs(probs1-probs2) > cutoff:
                element = [z,probs1,probs2]
                list.append(element)
                element = []

    return list


'''
displayTextResults(commonalities, differences)
#5 [Check6-2]
Parameters: 2D list of strs ; 2D list of values
Returns: None
'''
def displayTextResults(commonalities, differences):
    print("most common proteins :")
    for i in commonalities:
        for j in i:
            if j not in ['Start','Stop']:
                print(j)
    print("most-different amino acids :")
    for j in differences:
        print(f"{j[0]} : {round(j[1] * 100,2)}% in probs1, {round(j[2] * 100,2)}% in probs2")
        
    return
def runWeek2():
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")

    commonalities = commonProteins(humanProteins, elephantProteins)
    differences = findAminoAcidDifferences(humanProteins, elephantProteins, 0.005)
    displayTextResults(commonalities, differences)


## # WEEK 3 ###

'''
makeAminoAcidLabels(proteinList1, proteinList2)
#2 [Hw6]
Parameters: 2D list of strs ; 2D list of strs
Returns: list of strs
'''
def makeAminoAcidLabels(proteinList1, proteinList2):
    combine_list1 = combineProteins(proteinList1)
    combine_list2 = combineProteins(proteinList2)
    combine_list3 = combine_list1 + combine_list2
    amino_acid = aminoAcidDictionary(combine_list3)
    result = list(amino_acid.keys())
    result.sort()
    return result


'''
setupChartData(labels, proteinList)
#3 [Hw6]
Parameters: list of strs ; 2D list of strs
Returns: list of floats
'''
def setupChartData(labels, proteinList):
    freq = []
    combine_list = combineProteins(proteinList)
    amino_acid = aminoAcidDictionary(combine_list)
    for i in labels:
        if i in amino_acid:
            freq.append(amino_acid[i]/len(combine_list))
        else:
            freq.append(0)
    return freq


'''
createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None)
#4 [Hw6] & #5 [Hw6]
Parameters: list of strs ; list of floats ; str ; list of floats ; str ; [optional] list of strs
Returns: None
'''
def createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None):
    import matplotlib.pyplot as plt
    import numpy as np
    bar_width = 0.35

    # Create the x values for the bar charts
    x_values1 = np.arange(len(xLabels)) - bar_width/2
    x_values2 = np.arange(len(xLabels)) + bar_width/2

    # Plot the first set of bars
    plt.bar(x_values1, freqList1, width=bar_width, label=label1,edgecolor = edgeList)

    # Plot the second set of bars
    plt.bar(x_values2, freqList2, width=bar_width, label=label2,edgecolor=edgeList)

    # Set x-axis ticks and labels
    plt.xticks(np.arange(len(xLabels)), xLabels)

    # Set labels and title
    plt.xlabel('Amino Acid Labels')
    plt.ylabel('Frequency')
    plt.title('Gene Frequency Comparison')

    # Add legend
    plt.legend()

    # Show the plot
    plt.show()

    return


'''
makeEdgeList(labels, biggestDiffs)
#5 [Hw6]
Parameters: list of strs ; 2D list of values
Returns: list of strs
'''
def makeEdgeList(labels, biggestDiffs):
    edgeList = [] 
    temp =[]
    for i in biggestDiffs:
            temp.append(i[0])
    for i in labels:
            if i in temp:
                edgeList.append('black')
            else:
                edgeList.append('white')

    return edgeList


'''
runFullProgram()
#6 [Hw6]
Parameters: no parameters
Returns: None
'''
import json
import matplotlib.pyplot as plt
import numpy as np

def runFullProgram():
    #1 Synthesize Proteins
    human_proteins = synthesizeProteins('data/human_p53.txt', 'data/codon_table.json')
    elephant_proteins = synthesizeProteins('data/elephant_p53.txt', 'data/codon_table.json')
    #  2: Text Report
    common_proteins = commonProteins(human_proteins, elephant_proteins)
    differences = findAminoAcidDifferences(human_proteins, elephant_proteins, 0.005)
    displayTextResults(common_proteins, differences)
    # 3: Generate Bar Chart
    x_labels = makeAminoAcidLabels(human_proteins, elephant_proteins)
    freq_list_human = setupChartData(x_labels, human_proteins)
    freq_list_elephant = setupChartData(x_labels, elephant_proteins)
    edge_list = makeEdgeList(x_labels, differences)
    createChart(x_labels, freq_list_human, 'Human p53', freq_list_elephant, 'Elephant p53', edge_list)
    plt.show()
    return


### RUN CODE ###

# This code runs the test cases to check your work
if __name__ == "__main__":

     # print("\n" + "#"*15 + " WEEK 1 TESTS " +  "#" * 16 + "\n")
     # test.week1Tests()
     # print("\n" + "#"*15 + " WEEK 1 OUTPUT " + "#" * 15 + "\n")
     # runWeek1()

    ## Uncomment these for Week 2 ##
    
    # print("\n" + "#"*15 + " WEEK 2 TESTS " +  "#" * 16 + "\n")
    # test.week2Tests()
    # print("\n" + "#"*15 + " WEEK 2 OUTPUT " + "#" * 15 + "\n")
    # runWeek2()
    

    # ## Uncomment these for Week 3 ##
    
    # print("\n" + "#"*15 + " WEEK 3 TESTS " +  "#" * 16 + "\n")
    # test.week3Tests()
    # print("\n" + "#"*15 + " WEEK 3 OUTPUT " + "#" * 15 + "\n")
    runFullProgram()


#NAME:Sharma
#------------------------------------------------------------------------------------------
#PYTHON LAB 6
#------------------------------------------------------------------------------------------
#PART 1 -

#FUNCTION NAME: rna2protein
#PARAMETERS: 1 (An RNA sequence)
#PURPOSE: The function should:
#           (1) Divide the RNA sequence into a list of 3-base codons.
#HINT: You could use CodonList inside this function from previous lecture.
#           (2) Create a new protein string.
#           (3) Use the "standard_code" dictionary to find the amino acid for each codon in the list.
#           (4) Return the new protein string.
#RETURN VALUES: A protein sequence. (A string)

#Hint: When building the protein you are trying to return, add the dictionary values below (the letters
#      of the amino acids to an empty string (e.g., prot="").

standard_code = {
     "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L", "UCU": "S",
     "UCC": "S", "UCA": "S", "UCG": "S", "UAU": "Y", "UAC": "Y",
     "UAA": "*", "UAG": "*", "UGA": "*", "UGU": "C", "UGC": "C",
     "UGG": "W", "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
     "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P", "CAU": "H",
     "CAC": "H", "CAA": "Q", "CAG": "Q", "CGU": "R", "CGC": "R",
     "CGA": "R", "CGG": "R", "AUU": "I", "AUC": "I", "AUA": "I",
     "AUG": "M", "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
     "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K", "AGU": "S",
     "AGC": "S", "AGA": "R", "AGG": "R", "GUU": "V", "GUC": "V",
     "GUA": "V", "GUG": "V", "GCU": "A", "GCC": "A", "GCA": "A",
     "GCG": "A", "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
     "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G"}

def rna2protein(rna_seq):
     codon_list = []
     for number in range(0,len(rna_seq),3):
          codon = rna_seq[number:number + 3]
          if len(codon) == 3:
              codon_list.append(codon)
     new_protein = ""
     for codon in codon_list:
          if codon in standard_code:
               new_protein = new_protein + standard_code[codon]
     return new_protein

#EXAMPLE:
#protein=rna2protein("GCGAGGGUCUGA")
#print (protein)

#This should print:
#ARV*

#------------------------------------------------------------------------------------------
#PART 2 -

#FUNCTION NAME: dna2protein
#PARAMETERS: 1 (A DNA sequence)
#PURPOSE: The function should:
#           (1) Clean the DNA sequence and convert it to RNA. Just change T's to U's (don't reverse compliment)
#           (2) Divide the RNA sequence into a list of 3-base codons.
#           (3) Create a new protein string.
#           (4) Use the "standard_code" dictionary to find the amino acid for each codon.
#                   [NOTE: Stop translating after "Stop" codons.]
#           (5) Return the new protein string.
#RETURN VALUES: A protein sequence. (A string)


#== FUNCTION 2 ==
def dna2protein(dna_seq):
    clean_dna = dna_seq.strip().upper().replace(" ", '').replace("?","N")
    rna_seq = clean_dna.replace("T", "U")
    codon_list = []
    for number in range(0, len(rna_seq), 3):
        codon = rna_seq[number : number + 3]
        if len(codon) == 3:
            codon_list.append(codon)
    new_protein = ""
    for codon in codon_list:
        try :
            codon in standard_code
            new_protein = new_protein + standard_code[codon]
            if standard_code[codon] == "*":
                return new_protein
        except:
            new_protein = new_protein + "?"
    return new_protein


#EXAMPLE:
#protein=dna2protein("\n   ATGCaaaGAGacTGAgCC  \n\t\n")
#print (protein)

#The function would return:
#   MQRD*

#############Extra PRACTICE#############
# Modify your dna2protein function
# Use a try/except loop with your dictionary so that if the key is NOT in the dictionary,
#  It adds a "?" instead to the protein sequence.

#EXAMPLE:
#protein2=dna2protein("\n   ATGC-aaGAGacTGAgCC  \n\t\n")
#print (protein2)
###
#Should print:
#   M?RD*


#------------------------------------------------------------------------------------------

#== FUNCTION 3 ==

#FUNCTION NAME: EC_translate
#PARAMETERS: 2 (A sequence, and a default integer parameter set equal to "0")  
#PURPOSE: The function should:
#           (1) Clean the sequence.
#           (2) Check to see if the sequence is RNA or DNA. If it is DNA, convert it to RNA.
#           (3) Divide the RNA sequence into a list of 3-base codons.
#           (4) Create a new protein string.
#           (5) Based on the second parameter, use the appropriate dictionary 
#                   (if 0 use "standard_code", else use "mitochondrial_code")

#           (6) Return the new protein string.

#RETURN VALUES: A protein sequence. (A string)


mitochondrial_code = {
     "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L", "UCU": "S",
     "UCC": "S", "UCA": "S", "UCG": "S", "UAU": "Y", "UAC": "Y",
     "UAA": "*", "UAG": "*", "UGU": "C", "UGC": "C", "UGA": "W",
     "UGG": "W", "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
     "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P", "CAU": "H",
     "CAC": "H", "CAA": "Q", "CAG": "Q", "CGU": "R", "CGC": "R",
     "CGA": "R", "CGG": "R", "AUU": "I", "AUC": "I", "AUA": "M",
     "AUG": "M", "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
     "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K", "AGU": "S",
     "AGC": "S", "AGA": "*", "AGG": "*", "GUU": "V", "GUC": "V",
     "GUA": "V", "GUG": "V", "GCU": "A", "GCC": "A", "GCA": "A",
     "GCG": "A", "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
     "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G"}


#== FUNCTION 3 ==
def EC_translate(fseq, fcode=0):   #fcode is set to a default of 0 (see above)
    if fcode == 0:
        new_protein = dna2protein(fseq)
    else:
        clean_seq = fseq.strip().upper().replace(" ", '').replace("?", "N")
        rna_seq = clean_seq.replace("T", "U")
        codon_list = []
        for number in range(0, len(rna_seq), 3):
            codon = rna_seq[number: number + 3]
            if len(codon) == 3:
                codon_list.append(codon)
        new_protein = ""
        for codon in codon_list:
            try:
                codon in mitochondrial_code
                new_protein = new_protein + mitochondrial_code[codon]
                if mitochondrial_code[codon] == "*":
                    return new_protein
            except:
                new_protein = new_protein + "?"


    return new_protein


#HINT(1):
#   If NO second parameter is passed to the function (the default is not changed),
#       the function should use the "standard_code".

#EXAMPLE(1):
#protein1=EC_translate("\n  AUGccaaGAGActGAgCC \t\n")
#print (protein1)

#The function would return:
#   MPRD*
 
#HINT(2):
#   If a second parameter is passed to the function (the default is changed),
#       the function should use the "mitochondrial_code".

#EXAMPLE(2):
#protein2=EC_translate("\n  AUGccaaGAGActGAgCC \t\n", 1)
#print (protein2)

#The function would return:
#   MP*

#Test the function with the following sequences using both "standard" and "mitochondrial" codes:
#rna1=" \n \tAUGcaaGCAGuuACAUGAGagguAGGCAAGCACGCAGGAAC   \n\t"
#dna1=" \n atGTTCAtagTCATTATagTTacagTATTATtCTGa \n\t"

#print(EC_translate(rna1, 0))
#print(EC_translate(rna1, 1))
#print(EC_translate(dna1, 0))
#print(EC_translate(dna1, 1))

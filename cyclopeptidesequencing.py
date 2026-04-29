#This implementation uses the **Leaderboard Cyclopeptide Sequencing algorithm**, which:

 #Starts with an empty peptide  
 #Iteratively expands candidate peptides using possible amino acid masses  
 #Filters candidates based on total mass  
 #Scores candidates against the experimental spectrum  
 #Keeps only the top `N` best scoring candidates (leaderboard)  
 #Returns the best-scoring cyclic peptide(s) 
with open (r"Tyrocidine_B1_Spectrum_10.txt", "r") as f:
    
    N=1000
    spectrum=list(map(int,f.readline().strip().split()))

from collections import Counter
def cyclic_spectrum(peptide):
    prefix_mass=[0] 
    for m in peptide:
        prefix_mass.append(prefix_mass[-1]+m)
    peptide_mass=prefix_mass[-1]
    spectrum=[0]
    for i in range(len(peptide)):
        for j in range(i+1, len(peptide)+1):
             spectrum.append(prefix_mass[j]-prefix_mass[i])
             if i>0 and j<len(peptide):
                spectrum.append(peptide_mass - (prefix_mass[j] - prefix_mass[i]))
                                                                             
    return sorted(spectrum)
AMINO_ACID_MASSES = [

        57, 71, 87, 97, 99, 101, 103, 113, 114,

        115, 128, 129, 131, 137, 147, 156, 163, 186

    ]
def linear_spectrum(peptide):
    prefix_mass=[0]
    for m in peptide:
        prefix_mass.append(prefix_mass[-1]+m)
    spectrum=[0]
    for i in range(len(peptide)):
        for j in range(i+1, len(peptide)+1):
            spectrum.append(prefix_mass[j]-prefix_mass[i])
    return sorted(spectrum)
def expand(peptides):
    return [peptide + [mass] for peptide in peptides for mass in AMINO_ACID_MASSES]

def scoring(amino,table, spectrum):
    score=0
    amino_spectrum=linear_spectrum(amino)
    
    amino_spectrum_count=Counter(amino_spectrum)
    spectrum_count=Counter(spectrum)
    
    
    
    for i in amino_spectrum_count.keys():
        if i in spectrum_count:
            a=min(amino_spectrum_count[i],spectrum_count[i])
            score+=a
    return score

def cyclic_score(amino,table,spectrum):
    cyclic_score=0
    cyclic_spect=cyclic_spectrum(amino)
    cyclic_spect_count=Counter(cyclic_spect)
    spectrum_count=Counter(spectrum)
    for i in cyclic_spect_count.keys():
        if i in spectrum_count.keys():
            b=min(cyclic_spect_count[i],spectrum_count[i])
            cyclic_score+=b
    return cyclic_score
    
    

def leaderboardCyclopeptideSequencing(spectrum,N):
    Leaderboard=[[]]
    leaderpeptide=[]
    n_top=[]
    parent_mass=spectrum[-1]
    while Leaderboard:
        Leaderboard=expand(Leaderboard)
        Leaderboard=[peptide for peptide in Leaderboard if sum(peptide)<=parent_mass]
        for peptide in Leaderboard:
            if sum(peptide)==parent_mass:
                if cyclic_score(peptide,table,spectrum)>cyclic_score(leaderpeptide,table,spectrum):
                    leaderpeptide=peptide
                    n_top.append((leaderpeptide,cyclic_score(leaderpeptide,table,spectrum)))
        scored_peptides=[(peptide, scoring(peptide,table,spectrum)) for peptide in Leaderboard]
        scored_peptides.sort(key=lambda x:x[1], reverse=True)
        if len(scored_peptides)>N:
            cutoff_score=scored_peptides[N-1][1]
            Leaderboard=[p for p, s in scored_peptides if s>=cutoff_score]
        else:
            Leaderboard = [p for p, s in scored_peptides]
            
        if not Leaderboard:
            break
    n_top.sort(key=lambda x:x[1], reverse=True)
    return n_top
    
table = {
    "G": 57,
    "A": 71,
    "S": 87,
    "P": 97,
    "V": 99,
    "T": 101,
    "C": 103,
    "I": 113,
    "L": 113,
    "N": 114,
    "D": 115,
    "K": 128,
    "Q": 128,
    "E": 129,
    "M": 131,
    "H": 137,
    "F": 147,
    "R": 156,
    "Y": 163,
    "W": 186}              
print(leaderboardCyclopeptideSequencing(spectrum, N))    

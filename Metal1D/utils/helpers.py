import sys
import numpy as np
import pandas as pd
import py3Dmol
from biopandas.pdb import PandasPdb
from scipy.spatial import distance_matrix


def ExtractProbMap(resultsCOORDfile, coord_toexclude = 1):
  """ Extract probability map from coordination file
      obtained with LINK analysis (output of LINKanalysis.py )
    
    Parameters
    ----------
    resultsCOORDfile : str
      path of the probability map file of interest

    coord_toexclude : int
      coordinations equal or lower than COORD_THRESHOLD are excluded
        
    Returns
    -------
    ProbMap : list
      Probability map containing  [0] coordination (1-letter codes for AA)
                                  [1] probability 
        
    """
  CoordFile = open(resultsCOORDfile, 'r')
  COORD_THRESHOLD = coord_toexclude
  ln = 0
  ProbMap = [[],[]]
  for line in CoordFile:
      fileline = line.split(';')[:-1]
      if(ln==0):
          for c in range(1,len(fileline)):
              ProbMap[0].append(fileline[c])
      else:
          for c in range(1,len(fileline)):
              if(ln==1):
                  ProbMap[1].append(int(fileline[c]))
              else:
                  ProbMap[1][c-1]+=int(fileline[c])
      ln+=1 
  toremove = []
  for i in range(0,len(ProbMap[0])):    
      if(len(ProbMap[0][i])<=COORD_THRESHOLD):
          toremove.append(i)

  for i in range(0,len(toremove)):
      index = toremove[i]-i
      ProbMap[0].pop(index)
      ProbMap[1].pop(index)

  CoordFile.close()
  #PrbMap now contains coordinations observed with corresponding occurrency
  Tot = sum(ProbMap[1])
  for i in range(0, len(ProbMap[1])):
      ProbMap[1][i] /= Tot
  #PrbMap now contains coordinations observed with corresponding probability

  return ProbMap

def ProteinRead(pdb_file, Include_dAA = True, IncludeWATER = False):
  """ Read pdb file using biopandas library       
      excluding non-standard amino acids
    
    Parameters
    ----------
    pdb_file : str
      name of pdb file to analyze

    Include_dAA : bool
      Flag for D-Amino acids:
        if True, considered with probabilities of corresponding L-ones
        if False, excluded (considered as non-standard) 

    IncludeWATER : bool
      Flag for water molecules present in pdb structure 
      probability map from LINK analysis contains waters (when resolved)
        if True water molecules considered and scored as done with other AA 
        if False water molecules in the structure excluded
      Note: False by default due to low reliability of water location in structures
        used to build the probability map
        and to difficulty of placing water molecules (expecially buried ones) 
        in new structures, expecially when also metal ions are present   
        
    Returns
    -------
    ppdb_ATOM : pandas DataFrame
      biopandas dataset with amino acids (and water) of the input structure considered 

    Chains : list
      chains in structure (used for scoring part)

    """
  # structure from input file or fetched if not present
  if(pdb_file[-4:] == '.pdb' or pdb_file[-3:] == '.gz'):
      ppdb = PandasPdb().read_pdb(pdb_file)
  else:
      ppdb = PandasPdb().fetch_pdb(pdb_file)
  
  # lists for standard and d-AA used to save structure to dataset 
  standardAA = ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']
  d_AA = ['DAL','DAR','DSG','DAS','DCY','DGN','DGL','GLY','DHI','DIL','DLE','DLY','MED','DPN','DPR','DSN','DTH','DTR','DTY','DVA']#scan takes into account only standard amino acids

  for aa in standardAA: #ATOM entries, excluding water molecules 
      if(aa==standardAA[0]):
          ppdb_ATOM = ppdb.df['ATOM'][ppdb.df['ATOM']['residue_name'] == aa] 
      else:
          ppdb_ATOM = pd.concat([ppdb_ATOM, ppdb.df['ATOM'][ppdb.df['ATOM']['residue_name'] == aa]], ignore_index=True) 

  if(Include_dAA):
      for i in range(0,len(d_AA)): 
          if(d_AA[i]!='GLY'):
              ppdb_d_AA = pd.concat([ppdb.df['ATOM'][ppdb.df['ATOM']['residue_name'] == d_AA[i]],ppdb.df['HETATM'][ppdb.df['HETATM']['residue_name'] == d_AA[i]]], ignore_index=True)
              pd.options.mode.chained_assignment = None 
              ppdb_d_AA['residue_name'].iloc[:] = standardAA[i] #dAA considered as standard one for scan 
              ppdb_ATOM = pd.concat([ppdb_ATOM, ppdb_d_AA], ignore_index=True) 

  ppdb_PROTEIN = ppdb_ATOM #protein atoms saved here 
  ppdb_WATER = pd.concat([ppdb.df['HETATM'][ppdb.df['HETATM']['residue_name'] == 'HOH'],ppdb.df['ATOM'][ppdb.df['ATOM']['residue_name'] == 'HOH'],ppdb.df['HETATM'][ppdb.df['HETATM']['residue_name'] == 'WAT'],ppdb.df['ATOM'][ppdb.df['ATOM']['residue_name'] == 'WAT']], ignore_index=True) #oxygen atoms of water molecules
                                                                                                                                                                                                                                                                                              #can be both HETATM (standard pdb file) or ATOM (vmd output)
  if(len(ppdb_WATER)>0 and IncludeWATER):
      pd.options.mode.chained_assignment = None 
      ppdb_WATER['residue_name'].iloc[:] = 'HOH'
      ppdb_WATER['chain_id'].iloc[:] = 'water'
      ppdb_ATOM = pd.concat([ppdb_ATOM, ppdb_WATER], ignore_index=True)

  Chains = []
  for i in range(0,len(ppdb_ATOM)):
      if(ppdb_ATOM['chain_id'].iloc[i] in Chains):
          continue
      else:
          Chains.append(ppdb_ATOM['chain_id'].iloc[i])  
  return ppdb_ATOM, Chains

def RefAtom(Residue):
  """ Assign reference atom(s) from which to search for neighbouring residues

    Parameters
    ----------
    Residue : str
      3-letter code for AA  
  
    Returns
    -------
    REF : list
      reference atom(s) for given residue
    """

  RES =  ['HOH','ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']
  REF = [['O'],['O'], ['NH1','NH2'], ['OD1'], ['OD1', 'OD2'], ['SG'], ['OE1'], ['OE1', 'OE2'], ['O'], ['ND1', 'NE2'], ['O'], ['O'], ['NZ'], ['SD'], ['O'], ['O'], ['OG'], ['OG1'], ['O'], ['OH'], ['O']]
  return REF[RES.index(Residue)][:]

def res_1Letter(RES):
  """ Convert 3letter code AA to 1letter code AA

    Parameters
    ----------
    RES : str
      3-letter code AA  
  
    Returns
    -------
    R_LetterCode : str
      1letter code AA
    """
  RES_LetterCode = ['HOH', 'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']
  R_LetterCode = ['O', 'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']	
  return R_LetterCode[RES_LetterCode.index(RES)]

def ScoreLoc(Loc, ProbMap):
  """ Score a location based on the probability map
      score computed as sum of probabilities in ProbMap
      of coordinations compatible with the observed one
    
    Parameters
    ----------
    Loc : str
      local coordination in the structure

    ProbMap : list
      Probability map containing  [0] coordination (1-letter codes for AA)
                                  [1] probability 

    Returns
    -------
    Score : float
      score [0,1] of the given Loc

    """
  Score = 0.
  for i in range(0,len(ProbMap[0])):
    PotentialCoord = ProbMap[0][i]
    LocTest = Loc
    Flag = False
    for c in PotentialCoord:
      if c in LocTest:
        LocTest = LocTest.replace(c, '', 1)
        Flag = True
      else:
        Flag = False
        break
    if (Flag):
      Score += ProbMap[1][i]
  return Score

 # Chemical element specified in output file
  # Site prediction for each group
  # Weighted average reported in output file (.xyz file)
  # If 3 arguments are given last = chemical element
  ChemicalElement = 'H' # Default = H    
  if(len(sys.argv)==4):
      ChemicalElement = sys.argv[-1].upper()

  console_output = sys.stdout     
  if(sys.argv[1][-4:] == '.pdb'):
      OutFile = open(sys.argv[1][:-4].replace('../','')+'_PredictedSites.xyz', 'w') #overwrite file if already present
  elif(sys.argv[1][-3:] == '.gz'):
      OutFile = open(sys.argv[1][:-7].replace('../','')+'_PredictedSites.xyz', 'w') #overwrite file if already present
  else:
      OutFile = open(sys.argv[1].replace('../','')+'_PredictedSites.xyz', 'w') #overwrite file if already present
  OutFile.close()


def SitesPredict(Chain, pdb_file, ppdb_ATOM, ProbMap, ScoreThreshold = 0.75, SearchRadius = 5.5, ChemicalElement = 'H'):
  """ Scan all residues in each chain
      (scoring based on nearby residue in the same and other chains)
      and site determination as a weighted average between the residues with higher score
      where weights are given from scores
      Each predictions is written in a temporary output file 
    
    Parameters
    ----------
    Chain: str
      code for chain to scan

    pdb_file : str
      name of pdb file to analyze

    ProbMap : list
      Probability map containing  [0] coordination (1-letter codes for AA)
                                  [1] probability 
    ScoreThreshold : float
      prediction done considering residues within ScoreThreshold% 
      of the highest-scored one
      Default (0.75) resulted to be the best compromise 
      between sites found and false positives for ZN testset             
     
    SearchRadius : float
      Radius used to perform search around each amino acid
      For a given metal with d_{M-L} = D --> SearchRadius = 2*D + eps
      where eps accounts for possible rearrangements/structure relaxation
      Default (5.5) from average LINK distance of 2.2+-0.2 for ZN structures

    ChemicalElement : str
      Chemical element in output file
      Default H, changing it does not affect prediction
      only for visualization (proper vdw radius)

    """

  ResStart = min(ppdb_ATOM['residue_number'][ppdb_ATOM['chain_id']==Chain]) #First residue number
  ResStop = max(ppdb_ATOM['residue_number'][ppdb_ATOM['chain_id']==Chain])  #Last residue number
  LOCAL_coord = []
  Score = [[],[]]

  console_output = sys.stdout     
  if(pdb_file[-4:] == '.pdb'):
      OutFile = open(pdb_file[:-4].replace('../','')+'_PredictedSites.xyz', 'a') #overwrite file if already present
  elif(pdb_file[-3:] == '.gz'):
      OutFile = open(pdb_file[:-7].replace('../','')+'_PredictedSites.xyz', 'a') #overwrite file if already present
  else:
      OutFile = open(pdb_file.replace('../','')+'_PredictedSites.xyz', 'a') #overwrite file if already present
  OutFile.close()

  #Scan all protein structure assigning scores to each residue
  for i in range(ResStart, ResStop+1):
    # Scan one residue per time    
    # Missing residues present in the pdb structure are skipped
    try:  
      ppdb_RES = ppdb_ATOM[ppdb_ATOM['residue_number']==i] # Atoms of residue considered
      ppdb_RES = ppdb_RES[ppdb_RES['chain_id']==Chain]
      RES = ppdb_RES['residue_name'].iloc[0]
      REF = RefAtom(RES) # reference atoms(s) for residue considered   
      reference_coord = [0.0, 0.0, 0.0]
      add_SearchRadius = 0.0  # For residues with more reference points
                                    # SearchRadius augmented, adding distance from midpoint
      if(len(REF)==1):        
        reference_coord = [ppdb_RES['x_coord'][ppdb_RES['atom_name']==str(REF[0])].iloc[0], ppdb_RES['y_coord'][ppdb_RES['atom_name']==str(REF[0])].iloc[0], ppdb_RES['z_coord'][ppdb_RES['atom_name']==str(REF[0])].iloc[0]]            
      else:
      # if more than one atom is used as reference, the mean position is calculated as reference point
        for r in range(0,len(REF)):
          reference_coord[0] += (ppdb_RES['x_coord'][ppdb_RES['atom_name']==str(REF[r])].iloc[0])
          reference_coord[1] += (ppdb_RES['y_coord'][ppdb_RES['atom_name']==str(REF[r])].iloc[0])
          reference_coord[2] += (ppdb_RES['z_coord'][ppdb_RES['atom_name']==str(REF[r])].iloc[0]) 
        reference_coord[0] /= len(REF)
        reference_coord[1] /= len(REF)
        reference_coord[2] /= len(REF)
        dist = 0.0
        for r in range(0,len(REF)):
          dist += (reference_coord[0]-(ppdb_RES['x_coord'][ppdb_RES['atom_name']==str(REF[r])].iloc[0]))**2+(reference_coord[1]-(ppdb_RES['y_coord'][ppdb_RES['atom_name']==str(REF[r])].iloc[0]))**2+(reference_coord[2]-(ppdb_RES['z_coord'][ppdb_RES['atom_name']==str(REF[r])].iloc[0]))**2     
          dist = np.sqrt(dist)
          add_SearchRadius += dist
        add_SearchRadius /= len(REF)
          
      # Distance from reference point to any other atom of the structure
      distances = PandasPdb.distance_df(ppdb_ATOM, xyz=(reference_coord[0],reference_coord[1],reference_coord[2]))
      ppdb_LOCAL = ppdb_ATOM[distances <= (SearchRadius + add_SearchRadius)]
      ppdb_LOCAL = ppdb_LOCAL[ppdb_LOCAL['residue_number']!=i]	
      
      # Save local environment for residue considered
      # Other residues in LOCAL if their reference point(s) distance from the reference_coord
      # within SearchRadius+add_SearchRadius threshold
      LOCAL = [[RES],[i]]
      for l in ppdb_LOCAL['atom_number']:
          l_index = int(ppdb_LOCAL[ppdb_LOCAL['atom_number']==l].index.to_numpy())
          if(ppdb_LOCAL['atom_name'].loc[l_index] in RefAtom(ppdb_LOCAL['residue_name'].loc[l_index])):
              LOCAL[0].append(ppdb_LOCAL['residue_name'].loc[l_index])
              LOCAL[1].append(ppdb_LOCAL['residue_number'].loc[l_index])

      # Conversion to 1-letter format, alphabetically sorted
      for l in range(0,len(LOCAL[0])):
          LOCAL[0][l] = res_1Letter(LOCAL[0][l])
      LOCAL_coord.append(''.join(sorted(LOCAL[0])))
      # Score vector with score for each residue considered
      # assigned from ProbMap
      Score[0].append(i)
      Score[1].append(ScoreLoc(LOCAL_coord[-1], ProbMap))
    except:
      continue
  sortedlist = sorted(zip(Score[1][:],Score[0][:]), reverse=True)
 
  coeffScore = 1. - ScoreThreshold 
  highestscoredRES = [] # residue number, score, ref. point coord    
  for i in range(0,len(Score[0])):
    if(sortedlist[[i][0]][0]>=coeffScore*sortedlist[[0][0]][0]):
      x, y, z = 0.0, 0.0, 0.0
      ppdb_RES = ppdb_ATOM[ppdb_ATOM['residue_number']==sortedlist[[i][0]][1]]
      ppdb_RES = ppdb_RES[ppdb_RES['chain_id']==Chain]
      for l in range(0,len(RefAtom(ppdb_RES['residue_name'].iloc[0]))):
        # last iloc[0] to consider only first in case of alt_loc
        x += float(ppdb_RES['x_coord'][ppdb_RES['atom_name']==RefAtom(ppdb_RES['residue_name'].iloc[0])[l]].iloc[0])
        y += float(ppdb_RES['y_coord'][ppdb_RES['atom_name']==RefAtom(ppdb_RES['residue_name'].iloc[0])[l]].iloc[0])
        z += float(ppdb_RES['z_coord'][ppdb_RES['atom_name']==RefAtom(ppdb_RES['residue_name'].iloc[0])[l]].iloc[0])
      x /= len(RefAtom(ppdb_RES['residue_name'].iloc[0]))
      y /= len(RefAtom(ppdb_RES['residue_name'].iloc[0]))
      z /= len(RefAtom(ppdb_RES['residue_name'].iloc[0]))
      highestscoredRES.append([sortedlist[[i][0]][1],sortedlist[[i][0]][0], x, y, z])
  highestscoredRES_df = pd.DataFrame([row[2:5] for row in highestscoredRES], columns=['x','y','z'], index=[row[0] for row in highestscoredRES])
  highestscoredRES_dmatrix = pd.DataFrame(distance_matrix(highestscoredRES_df.values,highestscoredRES_df.values), index=highestscoredRES_df.index, columns=highestscoredRES_df.index)

  # nearby (closer than 2*Searchradius) residues grouped in clusters
  # producing 1 predicted site for each group
  Clusters = []
  possibleCluster = []
  for rw in range(0,len(highestscoredRES_dmatrix.index.values)):
    possibleCluster = []
    for cl in range(0,len(highestscoredRES_dmatrix.index.values)):
      if(any(highestscoredRES_dmatrix.index.values[cl] in row for row in Clusters)):
        continue 
      else:
        if(highestscoredRES_dmatrix.values[rw][cl]<2*SearchRadius):                    
          if(any(highestscoredRES_dmatrix.index.values[rw] in row for row in Clusters)):
            # residue added to group already existing
            Clusters[[Clusters.index(row) for row in Clusters if (highestscoredRES_dmatrix.index.values[rw]) in row][0]].append(highestscoredRES_dmatrix.index.values[cl])
          else:
            possibleCluster.append(highestscoredRES_dmatrix.index.values[cl])
    if(len(possibleCluster)>0):
      Clusters.append(possibleCluster)
  # highestscoredRES contains RES number, score, x_pos, y_pos, z_pos for each residue
  # in case of residues with more than one reference atom, midpoint is calculated
  for i in range(0,len(Clusters)):
        if(len(Clusters[i][:])==1):
            # append closest highest scored within threshold
            for j in range(0,len(Score[0])): 
                if(sortedlist[[j][0]][1] != Clusters[i][0]):
                    x, y, z = 0.0, 0.0, 0.0
                    ppdb_RES = ppdb_ATOM[ppdb_ATOM['residue_number']==sortedlist[[j][0]][1]]
                    ppdb_RES = ppdb_RES[ppdb_RES['chain_id']==Chain]
                    for l in range(0,len(RefAtom(ppdb_RES['residue_name'].iloc[0]))):
                        #last iloc[0] to consider only first in case of alt_loc
                        x += float(ppdb_RES['x_coord'][ppdb_RES['atom_name']==RefAtom(ppdb_RES['residue_name'].iloc[0])[l]].iloc[0])
                        y += float(ppdb_RES['y_coord'][ppdb_RES['atom_name']==RefAtom(ppdb_RES['residue_name'].iloc[0])[l]].iloc[0])
                        z += float(ppdb_RES['z_coord'][ppdb_RES['atom_name']==RefAtom(ppdb_RES['residue_name'].iloc[0])[l]].iloc[0])
                    x /= len(RefAtom(ppdb_RES['residue_name'].iloc[0]))
                    y /= len(RefAtom(ppdb_RES['residue_name'].iloc[0]))
                    z /= len(RefAtom(ppdb_RES['residue_name'].iloc[0]))
                    x_c = ([row[2] for row in highestscoredRES])[([row[0] for row in highestscoredRES]).index(Clusters[i][0])]
                    y_c = ([row[3] for row in highestscoredRES])[([row[0] for row in highestscoredRES]).index(Clusters[i][0])]
                    z_c = ([row[4] for row in highestscoredRES])[([row[0] for row in highestscoredRES]).index(Clusters[i][0])]
                    if(np.sqrt((x-x_c)**2+(y-y_c)**2+(z-z_c)**2)<2*(SearchRadius)):
                        Clusters[i].append(sortedlist[[j][0]][1])
                        if(sortedlist[[j][0]][0]>=coeffScore*([row[1] for row in highestscoredRES])[([row[0] for row in highestscoredRES]).index(Clusters[i][0])]):
                            highestscoredRES.append([sortedlist[[j][0]][1],sortedlist[[j][0]][0], x, y, z])
                        else:
                            highestscoredRES.append([sortedlist[[j][0]][1],coeffScore*([row[1] for row in highestscoredRES])[([row[0] for row in highestscoredRES]).index(Clusters[i][0])], x, y, z])
                        break 
                else:
                    continue
    # For each cluster xyz coord for site predicted as weighted average
    # between positions of RefAtom for each residue
    # weight depending on the score of the residue
    # Each predicted location is scored based on surrounding residues
    # within 0.6*SearchRadius 
    # if the site is farther than SearchRadius from any protein atom the site is discarded
    # (possible for sites predicted in bulk water)
  for j in range(0,len(Clusters)):
        x_avg, y_avg, z_avg = 0.0, 0.0, 0.0
        sumWeights = 0.0
        for k in range(0,len(Clusters[j])):
            resindex = [row[0] for row in highestscoredRES].index(Clusters[j][k])
            x_avg += highestscoredRES[resindex][1]*highestscoredRES[resindex][2]
            y_avg += highestscoredRES[resindex][1]*highestscoredRES[resindex][3]
            z_avg += highestscoredRES[resindex][1]*highestscoredRES[resindex][4]
            sumWeights += highestscoredRES[resindex][1]
        x_avg /= sumWeights 
        y_avg /= sumWeights 
        z_avg /= sumWeights
	# site discarded if farther than SearchRadius from any protein atom the site 
        distances = PandasPdb.distance_df(ppdb_ATOM, xyz=(x_avg,y_avg,z_avg))
        if(min(distances)>SearchRadius):
            break            	
        # during the loop coordinations temporary written in the output file
        # edited at the end of the loop according to xyz files formatting
        # and sorting predicted sites based on score
        if(pdb_file[-4:] == '.pdb'):
            sys.stdout = OutFile = open(pdb_file[:-4].replace('../','')+'_PredictedSites.xyz', 'a')
        elif(pdb_file[-3:] == '.gz'):
            sys.stdout = OutFile = open(pdb_file[:-7].replace('../','')+'_PredictedSites.xyz', 'a')
        else:
            sys.stdout = OutFile = open(pdb_file.replace('../','')+'_PredictedSites.xyz', 'a')

        sys.stdout.write('\n'+ChemicalElement+'\t'+str(x_avg)+'\t'+str(y_avg)+'\t'+str(z_avg))
        # Score the predicted location based on surrounding residues
        distances = PandasPdb.distance_df(ppdb_ATOM, xyz=(x_avg,y_avg,z_avg))
        ppdb_LOCAL = ppdb_ATOM[distances <= 0.6*SearchRadius]
        SITE_coord = ''
        for l in ppdb_LOCAL.index.tolist():
            if(ppdb_LOCAL['atom_name'].loc[l] in RefAtom(ppdb_LOCAL['residue_name'].loc[l])):
                SITE_coord += str(res_1Letter(ppdb_LOCAL['residue_name'].loc[l]))        
        # site score added last column in the temporary output file
        sys.stdout.write('\t'+str(ScoreLoc(''.join(sorted(SITE_coord)), ProbMap)))
        OutFile.close()
        sys.stdout = console_output



def SortPredictions(pdb_file, ScoreThreshold = 0.75):
  """ Sorts predictions in the temporary output file of SitesPredict function
    xyz file formatting
    sites sorted based on score (final score for each site as #comment in xyz file)
    
    Parameters
    ----------
    pdb_file : str
      name of pdb file to analyze
    
    ScoreThreshold : float
      final re-scoring of sites excludes sites with score lower than
      ScoreThreshold% of the highest-scored one
      Default (0.75) resulted to be the best compromise 
      between sites found and false positives for ZN testset             

     """
  console_output = sys.stdout     
  #read temporary OutFile    
  if(pdb_file[-4:] == '.pdb'):
      OutFile = open(pdb_file[:-4].replace('../','')+'_PredictedSites.xyz') 
  elif(pdb_file[-3:] == '.gz'):
      OutFile = open(pdb_file[:-7].replace('../','')+'_PredictedSites.xyz') 
  else:
      OutFile = open(pdb_file.replace('../','')+'_PredictedSites.xyz')
  PredictedSites = pd.read_table(OutFile, names = ['Element', 'x_coord', 'y_coord', 'z_coord', 'Score'])
  OutFile.close()
  #Open again for writing predicted sites location + score
  #according to xyz files formatting
  #sites sorted according to score, in descending order
  #score of the site added as comment after the coordinates (EL x_site y_site z_site #site_score)
  PredictedSites = PredictedSites.sort_values(by=['Score'], ascending=False)

  PredictedSites = PredictedSites[PredictedSites['Score']>(1-ScoreThreshold)*max(PredictedSites['Score'])] #site with score lower than ScoreThreshold% of higest one excluded
  Num_PredSites = len(PredictedSites)
  if(pdb_file[-4:] == '.pdb'):
      OutFile = open(pdb_file[:-4].replace('../','')+'_PredictedSites.xyz', 'w') 
  elif(pdb_file[-3:] == '.gz'):
      OutFile = open(pdb_file[:-7].replace('../','')+'_PredictedSites.xyz', 'w') 
  else:
      OutFile = open(pdb_file.replace('../','')+'_PredictedSites.xyz', 'w')
  OutFile.write(str(Num_PredSites)+'\n\n')
  #OutFile.write(str(Num_PredSites)+'\n'+'Copy and paste the following line in VMD representation to select highest-scored residues used to compute predicted sites:\n'+VMD_output+'\n')
  for i in PredictedSites.index:
      OutFile.write(str(PredictedSites.loc[i]['Element'])+'\t'+str(PredictedSites.loc[i]['x_coord'])+'\t'+str(PredictedSites.loc[i]['y_coord'])+'\t'+str(PredictedSites.loc[i]['z_coord'])+'\t#'+str(PredictedSites.loc[i]['Score'])+'\n')
  OutFile.close()

  sys.stdout = console_output
  print('----------')
  print('SCAN COMPLETED')
  print('\tPredicted sites can be found in:')
  if(pdb_file[-4:] == '.pdb'):
      print('\t'+pdb_file[:-4].replace('../','')+'_PredictedSites.xyz')
  elif(pdb_file[-3:] == '.gz'):
      print('\t'+pdb_file[:-7].replace('../','')+'_PredictedSites.xyz')
  else:
      print('\t'+pdb_file.replace('../','')+'_PredictedSites.xyz')
  print('----------')



def CreateOutFile(pdb_file):
  """ Output file initialization
    ovewrite any pre-existing file
    
    Parameters
    ----------
    pdb_file : str
      name of pdb file to analyze
     
    """   
  if(pdb_file[-4:] == '.pdb'):
    OutFile = open(pdb_file[:-4].replace('../','')+'_PredictedSites.xyz', 'w') #overwrite file if already present
  elif(pdb_file[-3:] == '.gz'):
    OutFile = open(pdb_file[:-7].replace('../','')+'_PredictedSites.xyz', 'w') #overwrite file if already present
  else:
    OutFile = open(pdb_file.replace('../','')+'_PredictedSites.xyz', 'w') #overwrite file if already present



def show_map(pdb,show_sticks_all=False, show_sticks_metalbinding=True, show_probes=True, show_pdb_metals=True):
    """ Show a protein using py3Dmol and the predicted metal sites
    
    Parameters
    ----------
    pdb : str
        Name of the pdb file
    show_sticks_all : bool
        If True, show all the residues as sticks in the protein
    show_sticks_metalbinding : bool
        If True, show the residues that are metal binding as sticks in the protein
    show_probes : bool
        If True, show predicted sites as filled spheres
        coloured according to probability 
    show_pdb_metals : bool
        If True, show metal ions in the structure (if present) as transparent spheres 

    Returns
    -------
    view : py3Dmol.view
        The view of the protein and the probes
    """
    view=py3Dmol.view(width=1000, height=800)

    view.addModel(open(pdb+'.pdb', 'r').read(),'pdb')
    if show_probes:
        view.addModel(open(pdb+'_PredictedSites.xyz', 'r').read(),'xyz')
        probes = open(pdb+'_PredictedSites.xyz', 'r').readlines()
        if(int(probes[0])!=0):
            probabilities = [p.replace('#','').split()[-1] for p in probes[2:]] # read p from comment in xyz file
            colors = {}
            # use different colors for the probabilities
            for i,x in enumerate(probabilities):
                colors[i] = '#%02x%02x%02x' % (0, 128, int(float(x)/float(probabilities[0])*255))
        else: #no predicted site
            colors = [] 
            view.addLabel("No probe predicted", {'position': {'x':0, 'y':0, 'z':0}, 'backgroundColor': '#0080FF', 'fontColor': 'white'});
   
    view.zoomTo()
    view.setBackgroundColor('white')
    view.setStyle({},{'cartoon': {'color':'gray'}})
    if show_sticks_all:
        view.setStyle({}, {'stick':{},'cartoon': {'color':'gray'}})
    if show_pdb_metals:
        view.getModel(0).setStyle({'resn':"ZN"},{'sphere': {'opacity':.75}})
        view.getModel(0).setStyle({'resn':"CA"},{'sphere': {'opacity':.75}})
        view.getModel(0).setStyle({'resn':"CU"},{'sphere': {'opacity':.75}})
        view.getModel(0).setStyle({'resn':"HG"},{'sphere': {'opacity':.75}})
        view.getModel(0).setStyle({'resn':"MG"},{'sphere': {'opacity':.75}})
        view.getModel(0).setStyle({'resn':"FE"},{'sphere': {'opacity':.75}})
        view.getModel(0).setStyle({'resn':"MN"},{'sphere': {'opacity':.75}})
        view.getModel(0).setStyle({'resn':"NI"},{'sphere': {'opacity':.75}})
        view.getModel(0).setStyle({'resn':"MB"},{'sphere': {'opacity':.75}})
    
    if show_probes:
        view.getModel(1).setStyle({},{'sphere': {'colorscheme':{'prop':'index', 'map':colors}}})
    
    # add hoverable labels for the residues and the predicted metals
    # two callbacks are needed, one for the residues and one for the metals
    # the metal one displays the probability
    view.getModel(0).setHoverable({},True,'''function(atom,viewer,event,container) {
                    if(!atom.label) {
                        atom.label = viewer.addLabel(atom.resn+atom.resi+":"+atom.atom,{position: atom, backgroundColor: 'mintcream', fontColor:'black'});
                    }}''',
                '''function(atom,viewer) { 
                    if(atom.label) {
                        viewer.removeLabel(atom.label);
                        delete atom.label;
                    }
                    }''')
    view.getModel(1).setHoverable({},True,'''function(atom,viewer,event,container) {
                    if(!atom.label) {
                        atom.label = viewer.addLabel(atom.atom+" ["+atom.serial+"]",{position: atom, backgroundColor: 'mintcream', fontColor:'black'});
                    }}''',
                '''function(atom,viewer) { 
                    if(atom.label) {
                        viewer.removeLabel(atom.label);
                        delete atom.label;
                    }
                    }''')
    if show_sticks_metalbinding:
        view.setStyle({'resn':"HIS"},{'stick': {}, 'cartoon': {'color':'gray'}})
        view.setStyle({'resn':"ASP"},{'stick': {}, 'cartoon': {'color':'gray'}})
        view.setStyle({'resn':"GLU"},{'stick': {}, 'cartoon': {'color':'gray'}})
        view.setStyle({'resn':"CYS"},{'stick': {}, 'cartoon': {'color':'gray'}})

    return view.show()
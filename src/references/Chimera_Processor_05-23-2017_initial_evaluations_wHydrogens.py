# Call this script from cygwin as: ./Chimera\ 1.10.1/bin/chimera.exe --start "Attribute Calculator" --start "Reply Log" --script "./Chimera_Processor_05-23-2017_initial_evaluations_wHydrogens.py"
# Or call it on Mac as: /Applications/Chimera.app/Contents/Resources/bin/chimera --start "Reply Log" --script "/Users/matthewrobinson/Documents/Radiodurans_Project/Chimera/Chimera_Processor_8-2-2016_Drad_Evaluation.py"
# Ubuntu/Debian: Add chimera to bash profile. Then you can just run: chimera --start "Attribute Calculator" --start "Reply Log" --script "./script.py"

import os
import chimera
import re
from chimera import runCommand as rc # use 'rc' as shorthand for runCommand
from chimera import replyobj # for emitting status messages
from chimera import selection
from chimera import dialogs
from DockPrep import prep
import time
import AddH
from DockPrep.prefs import prefs, defaults, INCOMPLETE_SC, MEMORIZED_SETTINGS

# Specify radius for Bubble feature computation in Angstroms.
radius = 12
outputDir = '/home/julian/chimera_test/original_output'

# Gather the names of .pdb files in the specified folder	
directoryName = "/home/julian/chimera_test"	
	
# To get all pdb files in the specified directory:	
#file_names = [fn for fn in os.listdir(directoryName) if fn.endswith(".pdb")][0]	
file_names = ["/home/julian/chimera_test/2zjr_C.pdb"]

# Make a logfile name
logfile = "log_07162019_initial_eval"


### EXAMPLE OF DEFINING A DICTIONARY	
# Dictionaries:	
# Amino acid 3-to-1 letter code dictionary.	
aaDict = dict({'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V', 'PTR': 'Y', 'CSX': 'C', 'CSO': 'C', 'CSD': 'C', 'KCX': 'K', 'UNK': 'X', 'MLY': 'K', 'PVH': 'H'})
# Residue number to 1 letter code dictionary.
resDict = dict({})	
MC_resDict = dict({})	

# Index 0: AAindPolarity
# Index 1: AAindSS
# Index 2: AAindMolVol
# Index 3: AAindCodonDiv
# Index 4: AAindElecCharge

AAindex_mat = [	
[-0.591, -1.302, -0.733, 1.57, -0.146],	
[-1.343, 0.465, -0.862, -1.020, -0.255],	
[1.05, 0.302, -3.656, -0.259, -3.242],	
[1.357, -1.453, 1.477, 0.113, -0.837],	
[-1.006, -0.590, 1.891, -0.397, 0.412],	
[-0.384, 1.652, 1.33, 1.045, 2.064],	
[0.336, -0.417, -1.673, -1.474, -0.078],	
[-1.239, -0.547, 2.131, 0.393, 0.816],	
[1.831, -0.561, 0.533, -0.277, 1.648],	
[-1.019, -0.987, -1.505, 1.266, -0.912],	
[-0.663, -1.524, 2.219, -1.005, 1.212],	
[0.945, 0.828, 1.299, -0.169, 0.933],	
[0.189, 2.081, -1.628, 0.421, -1.392],	
[0.931, -0.179, -3.005, -0.503, -1.853],	
[1.538, -0.055, 1.502, 0.44, 2.897],	
[-0.228, 1.399, -4.760, 0.67, -2.647],	
[-0.032, 0.326, 2.213, 0.908, 1.313],	
[-1.337, -0.279, -0.544, 1.242, -1.262],	
[-0.595, 0.009, 0.672, -2.128, -0.184],	
[0.26, 0.83, 3.097, -0.838, 1.512]	
]	

AAindex_dict = ({'A': AAindex_mat[0], 'C': AAindex_mat[1], 'D': AAindex_mat[2], 'E': AAindex_mat[3],	
 'F': AAindex_mat[4], 'G': AAindex_mat[5], 'H': AAindex_mat[6], 'I': AAindex_mat[7], 'K': AAindex_mat[8],	
 'L': AAindex_mat[9], 'M': AAindex_mat[10], 'N': AAindex_mat[11], 'P': AAindex_mat[12], 'Q': AAindex_mat[13],	
 'R': AAindex_mat[14], 'S': AAindex_mat[15], 'T': AAindex_mat[16], 'V': AAindex_mat[17], 'W': AAindex_mat[18],	
 'Y': AAindex_mat[19],})	
	
charge_index_dict = ({'A': 0, 'C': 0, 'D': -1, 'E': -1, 'F': 0, 'G': 0, 'H': 0, 'I': 0, 'K': 1,	
 'L': 0, 'M': 0, 'N': 0, 'P': 0, 'Q': 0, 'R': 1, 'S': 0, 'T': 0, 'V': 0, 'W': 0, 'Y': 0,})

OHrxnConst_index_dict = ({'A': 77000000, 'C': 34000000000, 'D': 75000000, 'E': 230000000, 'F': 6500000000, 'G': 17000000, 'H': 13000000000, 'I': 1800000000, 'K': 340000000,	
 'L': 1700000000, 'M': 8300000000, 'N': 49000000, 'P': 480000000, 'Q': 540000000, 'R': 3500000000, 'S': 320000000, 'T': 510000000, 'V': 760000000, 'W': 13000000000, 'Y': 13000000000,})

reactivity_index_dict = ({'A': 0, 'C': 0, 'D': 0, 'E': 0, 'F': 0, 'G': 0, 'H': 0, 'I': 0, 'K': 1.009211433,	
 'L': 0, 'M': 0, 'N': 0, 'P': 1.296184165, 'Q': 0, 'R': 0.900412937, 'S': 0, 'T': 0.830171969, 'V': 0, 'W': 0, 'Y': 0,})
		
# Set the number of times that we should retry when things fail
maxAttempts = 5

# Loop through the files, opening, processing, and closing each in turn.	
for fn in file_names:	
                print("Working with file: " + fn)

                # Add file name to the log file
                #log = open("log_Chimera_Processor_05_23", 'a')
                #log.write("Attempt " + attempt + ": " + fn + "\n")
                #log.close()

                # Open file in Chimera.	
                replyobj.status("Processing " + fn) # show what file we're working on	
                #rc("open " + directoryName + "/" + fn)	
                rc("open " + fn)	

                
        	
                # Parse input file name:	
                fn = fn.split("/")[0]
                #fn = re.sub(r'./', '', fn)	
                #fn = re.sub(r'.pdb', '', fn)
        	#fn = re.sub(r'\./allStructThatWork/', '', fn)
        	#fn = re.sub(r'allStructThatWork/', '', fn)
        	#fn = re.sub(r'allStructThatWork', '', fn)
                #fn = re.sub(r'allStruct', '', fn)
                #fn = re.sub(r'./Drad_GEM-PRO/','',fn)
                #fn = re.sub(r'./newStructVersions/','',fn)
        		
                # Output filenames:	
                areaSASFile = outputDir + '/chimeraVOiD_areaSAS_' + fn + '.txt'	
                f_areaSAS = open(str(areaSASFile), 'w')	
                depthFile = outputDir + '/chimeraVOiD_depth_' + fn + '.txt'	
                f_depth = open(str(depthFile), 'w')	
                distMatrixFile =  outputDir + '/chimeraVOiD_distMatrix_' + fn + '.txt'	
                f_distMatrix = open(distMatrixFile, 'w')	
                rlbcoorFile =  outputDir + '/' + fn + '.rlbcoor'	
                f_rlbcoor = open(rlbcoorFile, 'w')	
                ssFile =  outputDir + '/chimeraVOiD_SS_' + fn + '.txt'	
                f_SS = open(str(ssFile), 'w')	
                hydroFile =  outputDir + '/chimeraVOiD_hydro_' + fn + '.txt'	
                f_hydro = open(str(hydroFile), 'w')	
                surfaceMC_file =  outputDir + '/chimeraVOiD_surfaceMC_' + fn + '.txt'	
                f_surfaceMC = open(str(surfaceMC_file),'w')	
                cum_hydro_file =  outputDir + '/chimeraVOiD_cum_hydro_' + fn + '.txt'	
                f_cum_hydro = open(str(cum_hydro_file),'w')	
                num_rkpt_contacts_file =  outputDir + '/chimeraVOiD_num_rkpt_contacts_' + fn + '.txt'	
                f_num_rkpt_contacts = open(str(num_rkpt_contacts_file),'w')	
                num_contacts_file =  outputDir + '/chimeraVOiD_num_contacts_' + fn + '.txt'	
                f_num_contacts = open(str(num_contacts_file),'w')	
                AA_factors_file =  outputDir + '/chimeraVOiD_AA_factors_' + fn + '.txt'	
                f_AA_factors = open(str(AA_factors_file),'w')	
                charge_file =  outputDir + '/chimeraVOiD_net_charge_' + fn + '.txt'	
                f_charge = open(str(charge_file),'w')	
                OHrxnConst_file =  outputDir + '/chimeraVOiD_OH_rxn_const_' + fn + '.txt'	
                f_OHrxnConst = open(str(OHrxnConst_file),'w')	
                reactivity_file =  outputDir + '/chimeraVOiD_reactivity_' + fn + '.txt'	
                f_reactivity = open(str(reactivity_file),'w')	
                surface_MC_contacts_file =  outputDir + '/chimeraVOiD_surface_MCs_contacts_' + fn + '.txt'	
                f_MC_contacts = open(str(surface_MC_contacts_file),'w')	
        #        error_file = '.' + outputDir + '/chimeraVOiD_errors_' + fn + '.txt'	
        #        f_error_file = open(str(error_file),'w')	

        #	# FOR DEBUGGING, OTHERWISE COMMENT OUT
        #	debug_file = '.' + outputDir + '/chimeraVOiD_debug_' + fn + '.txt'
        #	f_debug = open(str(debug_file),'w')

        	
                # Clean out all non-peptide atoms and ligands.
                rc("select protein")
                rc("select invert")
                rc("delete selected")
        	rc("select ligand")
        	rc("delete selected")

                # Delete all alternative locations for residues.
                rc("select #0:@.B")
                rc("delete selected")

        # ADDING HYDROGENS AT THIS STEP LEADS TO TOO MANY FAILED COMPUTATIONS. INSTEAD, LET'S ADD THEM AFTER COMPUTING THE SURFACE.
        #	# Add hydrogen atoms/charges, mutate MSE and other non-standard AAs to standard, and perform other basic DockPrep operations.
        #	models = chimera.openModels.list(modelTypes=[chimera.Molecule])
        #	prep(models, addHFunc=AddH.hbondAddHydrogens, hisScheme=None, mutateMSE=True, mutate5BU=True, mutateUMS=True, mutateCSL=True, delSolvent=True, delIons=False, delLigands=False, delAltLocs=True, incompleteSideChains="rotamers", nogui=True, rotamerLib=defaults[INCOMPLETE_SC], rotamerPreserve=True, memorize=False, memorizeName=None)

                # Initialize areaSAS for each atom of protein, otherwise throws error for any atom without areaSAS attribute defined.
                rc("select protein")
                atoms = selection.currentAtoms()
                for a in atoms:
                        a.areaSAS = 0

                # Generate surface to protein.
                rc("select protein")
                rc("split")
        	r = dialogs.find('reply')
        	r.Clear()
                rc("surface allComponents false") # excludes bubbles from surface computation
        #        rc("surface probeRadius 2.8")

        	# Check the reply log for surface computation failure.
        	surfText = r.text.get('1.0', 'end')
        	surfTest = re.sub(r'connected surface components', '', surfText)
        	r.Clear()
        	# Unless there was a surface computed, skip to next structure file now.
        	if surfTest == surfText:
        		f_areaSAS.close()
        		os.remove('.' + outputDir + '/chimeraVOiD_areaSAS_' + fn + '.txt')
        		f_depth.close()
        		os.remove('.' + outputDir + '/chimeraVOiD_depth_' + fn + '.txt')
        		f_distMatrix.close()
        		os.remove('.' + outputDir + '/chimeraVOiD_distMatrix_' + fn + '.txt')
        		f_rlbcoor.close()
        		os.remove('.' + outputDir + '/' + fn + '.rlbcoor')
        		f_SS.close()
        		os.remove('.' + outputDir + '/chimeraVOiD_SS_' + fn + '.txt')
        		f_hydro.close()
        		os.remove('.' + outputDir + '/chimeraVOiD_hydro_' + fn + '.txt')
        		f_surfaceMC.close()
        		os.remove('.' + outputDir + '/chimeraVOiD_surfaceMC_' + fn + '.txt')
        		f_cum_hydro.close()
        		os.remove('.' + outputDir + '/chimeraVOiD_cum_hydro_' + fn + '.txt')
        		f_num_rkpt_contacts.close()
        		os.remove('.' + outputDir + '/chimeraVOiD_num_rkpt_contacts_' + fn + '.txt')
        		f_num_contacts.close()
        		os.remove('.' + outputDir + '/chimeraVOiD_num_contacts_' + fn + '.txt')
        		f_AA_factors.close()
        		os.remove('.' + outputDir + '/chimeraVOiD_AA_factors_' + fn + '.txt')
        		f_charge.close()
        		os.remove('.' + outputDir + '/chimeraVOiD_net_charge_' + fn + '.txt')
        		f_OHrxnConst.close()
        		os.remove('.' + outputDir + '/chimeraVOiD_OH_rxn_const_' + fn + '.txt')
        		f_reactivity.close()
        		os.remove('.' + outputDir + '/chimeraVOiD_reactivity_' + fn + '.txt')
        		f_MC_contacts.close()
        		os.remove('.' + outputDir + '/chimeraVOiD_surface_MCs_contacts_' + fn + '.txt')
        #		f_error_file.close()
        #		os.remove('.' + outputDir + '/chimeraVOiD_errors_' + fn + '.txt')

        		rc("close all")
        		continue

                rc("~select") # deselect the surface

        	# Add hydrogen atoms/charges, mutate MSE and other non-standard AAs to standard, and perform other basic DockPrep operations.
        	models = chimera.openModels.list(modelTypes=[chimera.Molecule])
        	prep(models, addHFunc=AddH.hbondAddHydrogens, hisScheme=None, mutateMSE=True, mutate5BU=True, mutateUMS=True, mutateCSL=True, delSolvent=True, delIons=False, delLigands=False, delAltLocs=True, incompleteSideChains="rotamers", nogui=True, rotamerLib=defaults[INCOMPLETE_SC], rotamerPreserve=True, memorize=False, memorizeName=None)


                # Compute global secondary structure and solvent accessible surface area attributes.
                sumHelixGlobal = 0
                sumSheetGlobal = 0
                sumSSGlobal = 0
                sumSASGlobal = 0

                rc("select protein")
                all_residues = selection.currentResidues()
                for res in all_residues:
                        if res.isHelix:
                                sumHelixGlobal += 1
                                sumSSGlobal += 1
                        if res.isSheet:
        		        sumSheetGlobal += 1
                                sumSSGlobal += 1
                        sumSASGlobal += res.areaSAS
                rc("~select") # deselect the surface

                # Select appropriate RKPT(W) atoms that are carbonylatable. R,K,P come from Requena 2003. T infered based on conversion of threonine to 2-amino-3-ketobutyric acid as noted in Levine 2001 and elsewhere. W comes from Taylor 2003.	
                # 'R_CD', 'K_CE', 'P_CD', 'T_CB', 'W_CD1'	
        #        rc("select :arg@cd|:lys@ce|:mly@ce|:kcx@ce|:pro@cd|:thr@cb|:trp@cd1")	
        #        rc("select protein") # old version that loops through all atoms but will be much slower
                rc("select :arg@cd|:lys@ce|:mly@ce|:kcx@ce|:pro@cd|:thr@cb z<" + str(radius) + ".1") # selects just atoms and residues within "radius".1 angstroms of carbonylatable RKPT atoms in order to cut down substantially on computation time when looping through to residues to find contacts for these atoms.
                all_atoms = selection.currentAtoms()
                all_residues = selection.currentResidues()	
                num_residues = len(all_residues)	
        	
                rc("select :arg@cd|:lys@ce|:mly@ce|:kcx@ce|:pro@cd|:thr@cb")	
                atoms = selection.currentAtoms()	
        	
        	
        	
                # COMPUTE SOLVENT ACCESSIBLE SURFACE AREA FOR EACH ATOM AND SS ASSIGNMENT AND HYDROPHOBICITY OF RESIDUE INCLUDING ATOM.	
                atomsList = []	
                	
                areaSASDict = dict({})	
                helixDict = dict({})	
                sheetDict = dict({})	
                hydroDict = dict({})	
        	
                num_rkpt_contacts_dict = dict({})	
                num_r_contacts_dict = dict({})	
                num_k_contacts_dict = dict({})	
                num_p_contacts_dict = dict({})	
                num_t_contacts_dict = dict({})	
        	
                for a in atoms:
        		areaSAS = a.areaSAS
        		residue = a.residue
        		type = residue.type
        		resHelix = residue.isHelix
        		resSheet = residue.isSheet
        		resHydro = residue.kdHydrophobicity
        		resNum = str(residue)
        		resNum = re.sub(r"#.+ .+ ", "", resNum)
        		# f_depth.write(resNum)
        		# f_depth.write("\n")
                        resDict[resNum] = str(aaDict[type])	
                        name = a.name	
                        altLoc = a.altLoc	
                        if altLoc == "":	
                                string = str(residue) + " " + str(aaDict[type]) + " " + str(name) + altLoc #what exactly is the name	
                        else:	
                                string = str(residue) + " " + str(aaDict[type]) + " " + str(name) + '.' + altLoc	
                        string = re.sub(r" ARG ", ":", string)	
                        string = re.sub(r" LYS ", ":", string)	
                        string = re.sub(r" MLY ", ":", string)	
                        string = re.sub(r" KCX ", ":", string)	
                        string = re.sub(r" PRO ", ":", string)	
                        string = re.sub(r" THR ", ":", string)	
                        #string = re.sub(r" TRP ", ":", string)
                        areaSASDict[string] = areaSAS	
                        helixDict[string] = resHelix	
                        sheetDict[string] = resSheet	
                        hydroDict[string] = resHydro	
                        atomsList.append(string)	
                atomsList = sorted(atomsList) #this doesn't seem to be working correctly	
        	
                # COMPUTE DEPTH OF EACH ATOM (i.e. minimum distance to surface)	
                depthDict = dict({}) #add methionine and cysteine in seperate list	
                minDepth = 1000000 # find what constitutes a depth of zero, should approximate radius of surface defining algorithm (1.8)
                rc("select #0")
                rc("~select #0:@") #deselecting the atoms, so we just get the surface	
        #        rc("measure distance :arg@cd|:lys@ce|:mly@ce|:kcx@ce|:pro@cd|:thr@cb|:trp@cd1 selection multiple true show true")	
                rc("measure distance :arg@cd|:lys@ce|:mly@ce|:kcx@ce|:pro@cd|:thr@cb selection multiple true show true") #from surface to all of these atoms	
                r = dialogs.find('reply')	
                text = r.text.get('1.0', 'end')	
                r.Clear()	
                f_temp = open('./chimeraVOiD_temp_' + fn + '.txt', 'w')	
                f_temp.write(text)	
                f_temp.close()	
        	
                with open('./chimeraVOiD_temp_' + fn + '.txt','r') as f_temp:	
                        for line in f_temp:	
                                line = line.rstrip()	
                                if line[:21] == 'minimum distance from':	

                                        line = line.replace('minimum distance from ', '')	
                                        line = re.sub(r" to #0:\? = ", "\t", line)	
                                        resNum = line	
                                        resNum = re.sub(r"#.+:", "", resNum)	
                                        resNum = re.sub(r"@.+", "", resNum)	
                                        line = re.sub(r'@', (" " + resDict[resNum] + " "), line)	
                                        string = re.sub(r"\t.+", "", line)
                                        depth = re.sub(r".+\t", "", line)
                                        depthDict[string] = depth
        	
                                        if float(depth) < minDepth:	
                                                minDepth = float(depth)	
                                                #print "minDepth is:" + str(minDepth)	
                os.remove('./chimeraVOiD_temp_' + fn + '.txt')	

                # COMPUTE INTERATOMIC DISTANCES. (Notes that the order of atoms selected is not deterministic and so varies across runs.)	
                distMatrix = [[0 for x in range(len(atoms))] for x in range(len(atoms))]	
                atomsIndex = dict({})	
                resTypeChainIndex = dict({})	
                xformCoordIndex = dict({})
                counter = 0	
                for a1 in atoms:	
                        residue = a1.residue	
                        type = residue.type	
                        name = a1.name	
                        altLoc = a1.altLoc	
                        if altLoc == "":	
                                string = str(residue) + " " + str(aaDict[type]) + " " + str(name) + altLoc	
                        else:	
                                string = str(residue) + " " + str(aaDict[type]) + " " + str(name) + '.' + altLoc	
                        string = re.sub(r" ARG ", ":", string)	
                        string = re.sub(r" LYS ", ":", string)	
                        string = re.sub(r" MLY ", ":", string)	
                        string = re.sub(r" KCX ", ":", string)	
                        string = re.sub(r" PRO ", ":", string)	
                        string = re.sub(r" THR ", ":", string)	
                        #string = re.sub(r" TRP ", ":", string)	
        	
                        atomsIndex[string] = counter	
                        counter = counter + 1	
        	
                counter1 = 0	
                for a1 in atoms:	
                        residue = a1.residue	
                        type = residue.type	
                        name = a1.name	
                        altLoc = a1.altLoc	
                        if altLoc == "":	
                                string = str(residue) + " " + str(aaDict[type]) + " " + str(name) + altLoc	
                        else:	
                                string = str(residue) + " " + str(aaDict[type]) + " " + str(name) + '.' + altLoc	
                        string = re.sub(r" ARG ", ":", string)	
                        string = re.sub(r" LYS ", ":", string)	
                        string = re.sub(r" MLY ", ":", string)	
                        string = re.sub(r" KCX ", ":", string)	
                        string = re.sub(r" PRO ", ":", string)	
                        string = re.sub(r" THR ", ":", string)	
                        #string = re.sub(r" TRP ", ":", string)	
        	
        	
                        # For rlbcoor data.	
                        resTypeChain = str(residue) + " A"	
                        resTypeChain = re.sub(r"#0 ", "",resTypeChain)	
                        resTypeChainIndex[string] =  str(resTypeChain)	
                        xformCoordIndex[string] = str(a1.xformCoord())	
                        counter2 = 0	
                        for a2 in atoms:	
                                distance = chimera.distance(a1.xformCoord(), a2.xformCoord())	
                                distMatrix[counter2][counter1] = str(distance)	
        	
                                counter2 = counter2 + 1	
        	
                        counter1 = counter1 + 1	
        	
                # PRINT OUT ALL RESULTS.	
                f_distMatrix.write("\t")	
                for a1 in atomsList:	
                        atom1 = str(a1)	
                        f_distMatrix.write(atom1)	
                        f_distMatrix.write("\t")	
                f_distMatrix.write("\n")	
        	
                for a1 in atomsList:	
                        atom1 = str(a1)	
                        depth = depthDict[atom1]	
                        hydro = hydroDict[atom1]	
        	
                        num_rkpt_contacts_dict[atom1] = 0 #initializes the values	
                        num_r_contacts_dict[atom1] = 0	
                        num_k_contacts_dict[atom1] = 0	
                        num_p_contacts_dict[atom1] = 0	
                        num_t_contacts_dict[atom1] = 0	
                        		
                        f_hydro.write(atom1)	
                        f_hydro.write("\t")	
                        f_hydro.write(str(hydro))	
                        f_hydro.write("\n")	
        	
                        f_depth.write(atom1)	
                        f_depth.write("\t")	
                        f_depth.write(str(depth))	
                        f_depth.write("\n")	
        	
                        f_distMatrix.write(atom1)	
                        f_distMatrix.write("\t")	
        	
                        for a2 in atomsList:	
                                atom2 = str(a2)	
                                a1Index = atomsIndex[atom1]	
                                a2Index = atomsIndex[atom2]	
                                distance = distMatrix[a2Index][a1Index]	
                                f_distMatrix.write(str(distance))	
                                f_distMatrix.write("\t")	
        	
                                dist = float(distance)	
        	
                                if dist > 0 and dist <= radius: 	
                                        num_rkpt_contacts_dict[atom1] += 1	
                                        #create num_k/r/p/t_contacts files too, just add 4 more conditionals	
                                        aa_letter = re.split(' ',atom2)[1]	
                                        if aa_letter == 'R':	
                                                num_r_contacts_dict[atom1] += 1	
                                        if aa_letter == 'K':	
                                                num_k_contacts_dict[atom1] += 1	
                                        if aa_letter == 'P':	
                                                num_p_contacts_dict[atom1] += 1	
                                        if aa_letter == 'T':	
                                                num_t_contacts_dict[atom1] += 1	
        	
                        f_distMatrix.write("\n")	
        	
                #write to num_rkpt_contacts	
                for a1 in atomsList:	
                        atom1 = str(a1)	
                	
                        num_rkpt_contacts = num_rkpt_contacts_dict[atom1]	
                        num_r_contacts = num_r_contacts_dict[atom1]	
                        num_k_contacts = num_k_contacts_dict[atom1]	
                        num_p_contacts = num_p_contacts_dict[atom1]	
                        num_t_contacts = num_t_contacts_dict[atom1]	
        	
                        f_num_rkpt_contacts.write(atom1)	
                        f_num_rkpt_contacts.write("\t")	
                        f_num_rkpt_contacts.write(str(num_rkpt_contacts))	
                        f_num_rkpt_contacts.write("\t")	
                        f_num_rkpt_contacts.write(str(num_r_contacts))	
                        f_num_rkpt_contacts.write("\t")	
                        f_num_rkpt_contacts.write(str(num_k_contacts))	
                        f_num_rkpt_contacts.write("\t")	
                        f_num_rkpt_contacts.write(str(num_p_contacts))	
                        f_num_rkpt_contacts.write("\t")	
                        f_num_rkpt_contacts.write(str(num_t_contacts))	
                        f_num_rkpt_contacts.write("\n")	
        		
                f_rlbcoor.write("\n")	
                f_rlbcoor.write("\n")	
                f_rlbcoor.write(" loop_\n")	
                f_rlbcoor.write("   _relibase_surface_pseudocentre_id\n")	
                f_rlbcoor.write("   _relibase_surface_pseudocentre_name\n")	
                f_rlbcoor.write("   _relibase_surface_pseudocentre_number\n")	
                f_rlbcoor.write("   _relibase_surface_pseudocentre_residue_type\n")	
                f_rlbcoor.write("   _relibase_surface_pseudocentre_residue_number\n")	
                f_rlbcoor.write("   _relibase_surface_pseudocentre_chainid\n")	
                f_rlbcoor.write("   _relibase_surface_pseudocentre_backbone\n")	
                f_rlbcoor.write("   _relibase_surface_pseudocentre_x\n")	
                f_rlbcoor.write("   _relibase_surface_pseudocentre_y\n")	
                f_rlbcoor.write("   _relibase_surface_pseudocentre_z\n")	
                counter2 = 0	
                for a1 in atomsList:	
                        atom1 = str(a1)	
                        depth = depthDict[atom1]	
                        label = ''	
                        num = 0	
                        backbone = 0 ##	
                        coords = xformCoordIndex[atom1]	
        	
                        xyzCoords = coords.split(' ',2)	
        	
                        if xyzCoords[0].find('.') > -1:	
                                xDecimal = xyzCoords[0].split('.',1)[1]	
                                xAddZeros = 4 - len(xDecimal)	
                                for i in range(xAddZeros):	
                                        xyzCoords[0] = xyzCoords[0] + "0"	
                        else:	
                                xyzCoords[0] = xyzCoords[0] + ".0000"	
                        if xyzCoords[1].find('.') > -1:	
                                yDecimal = xyzCoords[1].split('.',1)[1]	
                                yAddZeros = 4 - len(yDecimal)	
                                for i in range(yAddZeros):	
                                        xyzCoords[1] = xyzCoords[1] + "0"	
                        else:	
                                xyzCoords[1] = xyzCoords[1] + ".0000"	
                        if xyzCoords[2].find('.') > -1:	
                                zDecimal = xyzCoords[2].split('.',1)[1]	
                                zAddZeros = 4 - len(zDecimal)	
                                for i in range(zAddZeros):	
                                        xyzCoords[2] = xyzCoords[2] + "0"	
                        else:	
                                xyzCoords[2] = xyzCoords[2] + ".0000"	
        	
                        if float(depth) == float(minDepth):	
                                label = 'Donor'	
                                num = 0	
                        else:	
                                label = 'Acceptor'	
                                num = 1	
                        string = str(num) + " " + str(label) + " " + str(counter2) + " " + str(resTypeChainIndex[str(atom1)]) + " " + str(backbone) + " " + str(xyzCoords[0]) + " " + str(xyzCoords[1]) + " " + str(xyzCoords[2]) + "\n"	
                        f_rlbcoor.write(string)	
        	
                        counter2 = counter2 + 1	
        	
                f_distMatrix.close()	
                f_temp.close()	
                f_depth.close()	
                f_hydro.close()	
                f_num_rkpt_contacts.close()	
        	
        	
                # Compute contacts between RKPT atoms and any residue:	
                contacts_dict = dict({})	
                num_contacts_dict = dict({})	
        	
                error_count = 0	
                rc("select :arg@cd|:lys@ce|:mly@ce|:kcx@ce|:pro@cd|:thr@cb")
                atoms = selection.currentAtoms()
                rc("select :arg@cd|:lys@ce|:mly@ce|:kcx@ce|:pro@cd|:thr@cb z<12.1")
                all_residues = selection.currentResidues()
                
                for a1 in atoms:	
                        residue = a1.residue	
                        type = residue.type	
                        resNum = str(residue)	
                        resNum = re.sub(r"#.+ .+ ", "", resNum)	
        	
                        name = a1.name	
                        altLoc = a1.altLoc	
                        if altLoc == "":	
                                string = str(residue) + " " + str(aaDict[type]) + " " + str(name) + altLoc	
                        else:	
                                string = str(residue) + " " + str(aaDict[type]) + " " + str(name) + '.' + altLoc	
                        string = re.sub(r" ARG ", ":", string)	
                        string = re.sub(r" LYS ", ":", string)	
                        string = re.sub(r" MLY ", ":", string)	
                        string = re.sub(r" KCX ", ":", string)	
                        string = re.sub(r" PRO ", ":", string)	
                        string = re.sub(r" THR ", ":", string)	
                        #string = re.sub(r" TRP ", ":", string)	
        	
                        contacts_count = -1	
                        contacts_list = []	
        	
                        for res in all_residues: #Note that this inlcudes the a1.residue itself, but that is okay as long as we set contacts = -1 initially
                                res_atoms = res.atoms 	
                                for res_atom in res_atoms:	
                                        distance = chimera.distance(a1.xformCoord(), res_atom.xformCoord())	
                                        if distance <= radius and 'MSE' not in str(res) and 'CSX' not in str(res) and 'UNK' not in str(res) and 'CSD' not in str(res) and 'CSO' not in str(res) and 'KCX' not in str(res):	
                                                contacts_count += 1	
                                                contacts_list.append(res_atom.residue)	
                                               # if "116" in string:
                                                        # if 'PRO 18' in str(res):
                                                        #         print("""
                                                        #         The distance between {residue_atom} and {a1} is
                                                        #         {distance}. 
                                                        #         """.format(residue_atom = str(res_atom.residue),
                                                        #                 a1 = string,
                                                        #                 distance = distance))
                                                        #         print(a1.xformCoord())
                                                        #         print(res_atom.xformCoord())
                                                break #need to break the first enclosing for loop	
        	
                        num_contacts_dict[string] = contacts_count	
                        contacts_dict[string] = contacts_list	
                        print("On august 16 at 1146")
                        print(string)
                        print(contacts_count)
        	
        #                f_error_file.write(string)	
        #                f_error_file.write("\t")	
                        for res in contacts_dict[string]:	
        #                	f_error_file.write(str(res))
        #                	f_error_file.write("\t")
                        	error_count += 1
        #                f_error_file.write(str(error_count))	
        #                f_error_file.write("\n")	
        	
        #        f_error_file.close()	
        	
        	
                # Compute Bubble features
                cum_hydro_dict = dict({})	
                factors_dict = dict({})	
                net_charge_dict = dict({})
                pos_charge_dict = dict({})
                neg_charge_dict = dict({})
                OHrxnConst_dict = dict({})
                sumHelixLocal_dict = dict({})
                sumSheetLocal_dict = dict({})
                sumSSLocal_dict = dict({})
                sumSASLocal_dict = dict({})
                reactivity_dict = dict({})
        	
                for a1 in atomsList:	
                        atom1 = str(a1)	
                        aa_letter = re.split(' ',atom1)[1]	
                        print aa_letter	
        	
                        cum_hydro = 0 # include the atom's own hydrophobicity. Note that this works because contacts_list includes the residue itself.
                        factor1 = 0	
                        factor2 = 0	
                        factor3 = 0	
                        factor4 = 0	
                        factor5 = 0	
        	
                        net_charge = 0
                        pos_charge = 0
                        neg_charge = 0
                        OHrxnConst_sum = 0
                        sumHelixLocal = 0
                        sumSheetLocal = 0
                        sumSSLocal = 0
                        sumSASLocal = 0
                        reactivity_sum = 0
        	
                        all_residues = selection.currentResidues()
                        for res in all_residues:
                        #for res in contacts_dict[atom1]: #include a1's residue

        			kdHydro = res.kdHydrophobicity
        			# Handle exception non-standard amino acids that don't have kdHydrophobicity.
        			if str(res.type) == "MLY":
        				kdHydro = -3.9
        			if str(res.type) == "PTR":
        				kdHydro = -1.3
                                cum_hydro += kdHydro
                                print cum_hydro	
        	
                                res_letter = str(aaDict[res.type]) #in order to get right format for accessing dict, essentially same as aa_letter
                                print "res letter is:" + res_letter	
        	
                                net_charge += charge_index_dict[res_letter]

                                if charge_index_dict[res_letter] > 0:
        				pos_charge += charge_index_dict[res_letter]
                                if charge_index_dict[res_letter] < 0:
        				neg_charge += charge_index_dict[res_letter]

                                factor1 += AAindex_dict[res_letter][0]	
                                factor2 += AAindex_dict[res_letter][1]	
                                factor3 += AAindex_dict[res_letter][2]	
                                factor4 += AAindex_dict[res_letter][3]	
                                factor5 += AAindex_dict[res_letter][4]	

                                OHrxnConst_sum += OHrxnConst_index_dict[res_letter]

                                if res.isHelix:
        				sumHelixLocal += 1
        				sumSSLocal += 1
                                if res.isSheet:
        				sumSheetLocal += 1
        				sumSSLocal += 1

                                sumSASLocal += res.areaSAS

                                reactivity_sum += reactivity_index_dict[res_letter]

                        factors_dict[atom1] = [factor1, factor2, factor3, factor4, factor5]	
                        cum_hydro_dict[atom1] = cum_hydro	
                        net_charge_dict[atom1] = net_charge
                        pos_charge_dict[atom1] = pos_charge
                        neg_charge_dict[atom1] = neg_charge
                        OHrxnConst_dict[atom1] = OHrxnConst_sum
                        sumHelixLocal_dict[atom1] = sumHelixLocal
                        sumSheetLocal_dict[atom1] = sumSheetLocal
                        sumSSLocal_dict[atom1] = sumSSLocal
                        sumSASLocal_dict[atom1] = sumSASLocal
                        reactivity_dict[atom1] = reactivity_sum
                                                	
                for a1 in atomsList:	
                        atom1 = str(a1)	
                        num_contacts = num_contacts_dict[atom1]	
        	
                        cum_hydro = cum_hydro_dict[atom1]	
                        avg_cum_hydro = cum_hydro_dict[atom1]/(num_contacts+1)	
        	
                        factors = factors_dict[atom1]	
                        avg_factors = [factor / (num_contacts+1) for factor in factors_dict[atom1]]	
        	
                        net_charge = net_charge_dict[atom1]
                        avg_net_charge = float(net_charge_dict[atom1])/float((num_contacts+1))
                        pos_charge = pos_charge_dict[atom1]
                        avg_pos_charge = float(pos_charge_dict[atom1])/float((num_contacts+1))
                        neg_charge = neg_charge_dict[atom1]
                        avg_neg_charge = float(neg_charge_dict[atom1])/float((num_contacts+1))

                        OHrxnConst_sum = OHrxnConst_dict[atom1]
                        avg_OHrxnConst = float(OHrxnConst_dict[atom1])/float((num_contacts+1))

                        helix = helixDict[atom1]	
                        sheet = sheetDict[atom1]
                        sumHelixLocal = sumHelixLocal_dict[atom1]
                        sumSheetLocal = sumSheetLocal_dict[atom1]
                        sumSSLocal = sumSSLocal_dict[atom1]
                        avgHelixLocal = float(sumHelixLocal_dict[atom1])/float((num_contacts+1))
                        avgSheetLocal = float(sumSheetLocal_dict[atom1])/float((num_contacts+1))
                        avgSSLocal = float(sumSSLocal_dict[atom1])/float((num_contacts+1))

                        areaSAS = areaSASDict[atom1]
                        sumSASLocal = sumSASLocal_dict[atom1]
                        avgSASLocal = float(sumSASLocal_dict[atom1])/float((num_contacts+1))

                        reactivity_sum = reactivity_dict[atom1]
                        avg_reactivity = float(reactivity_dict[atom1])/float((num_contacts+1))
        	
                        f_cum_hydro.write(atom1)	
                        f_cum_hydro.write("\t")	
                        f_cum_hydro.write(str(cum_hydro))	
                        f_cum_hydro.write("\t")	
                        f_cum_hydro.write(str(avg_cum_hydro))	
                        f_cum_hydro.write("\n")	
        	
                        f_num_contacts.write(atom1)	
                        f_num_contacts.write("\t")	
                        f_num_contacts.write(str(num_contacts))	
                        f_num_contacts.write("\n")	
        	
                        f_AA_factors.write(atom1)	
                        f_AA_factors.write("\t")	
                        f_AA_factors.write(str(factors[0]))	
                        f_AA_factors.write("\t")	
                        f_AA_factors.write(str(factors[1]))	
                        f_AA_factors.write("\t")	
                        f_AA_factors.write(str(factors[2]))	
                        f_AA_factors.write("\t")	
                        f_AA_factors.write(str(factors[3]))	
                        f_AA_factors.write("\t")	
                        f_AA_factors.write(str(factors[4]))	
                        f_AA_factors.write("\t")	
                        f_AA_factors.write(str(avg_factors[0]))	
                        f_AA_factors.write("\t")	
                        f_AA_factors.write(str(avg_factors[1]))	
                        f_AA_factors.write("\t")	
                        f_AA_factors.write(str(avg_factors[2]))	
                        f_AA_factors.write("\t")	
                        f_AA_factors.write(str(avg_factors[3]))	
                        f_AA_factors.write("\t")	
                        f_AA_factors.write(str(avg_factors[4]))	
                        f_AA_factors.write("\n")	
        	
                        f_charge.write(atom1)	
                        f_charge.write("\t")	
                        f_charge.write(str(net_charge))	
                        f_charge.write("\t")	
                        f_charge.write(str(avg_net_charge))	
                        f_charge.write("\t")	
                        f_charge.write(str(pos_charge))	
                        f_charge.write("\t")	
                        f_charge.write(str(avg_pos_charge))	
                        f_charge.write("\t")	
                        f_charge.write(str(neg_charge))	
                        f_charge.write("\t")	
                        f_charge.write(str(avg_neg_charge))	
                        f_charge.write("\n")	

                        f_OHrxnConst.write(atom1)
                        f_OHrxnConst.write("\t")
                        f_OHrxnConst.write(str(OHrxnConst_sum))
                        f_OHrxnConst.write("\t")
                        f_OHrxnConst.write(str(avg_OHrxnConst))
                        f_OHrxnConst.write("\n")

                        f_SS.write(atom1)
                        f_SS.write("\t")
                        f_SS.write(str(helix))
                        f_SS.write("\t")
                        f_SS.write(str(sheet))
                        f_SS.write("\t")
                        f_SS.write(str(sumHelixLocal))
                        f_SS.write("\t")
                        f_SS.write(str(avgHelixLocal))
                        f_SS.write("\t")
                        f_SS.write(str(sumSheetLocal))
                        f_SS.write("\t")
                        f_SS.write(str(avgSheetLocal))
                        f_SS.write("\t")
                        f_SS.write(str(sumSSLocal))
                        f_SS.write("\t")
                        f_SS.write(str(avgSSLocal))
                        f_SS.write("\t")
                        f_SS.write(str(sumHelixGlobal))
                        f_SS.write("\t")
                        f_SS.write(str(sumSheetGlobal))
                        f_SS.write("\t")
                        f_SS.write(str(sumSSGlobal))
                        f_SS.write("\n")

                        f_areaSAS.write(atom1)
                        f_areaSAS.write("\t")
                        f_areaSAS.write(str(areaSAS))
                        f_areaSAS.write("\t")
                        f_areaSAS.write(str(sumSASLocal))
                        f_areaSAS.write("\t")
                        f_areaSAS.write(str(avgSASLocal))
                        f_areaSAS.write("\t")
                        f_areaSAS.write(str(sumSASGlobal))
                        f_areaSAS.write("\n")

                        f_reactivity.write(atom1)
                        f_reactivity.write("\t")
                        f_reactivity.write(str(reactivity_sum))
                        f_reactivity.write("\t")
                        f_reactivity.write(str(avg_reactivity))
                        f_reactivity.write("\n")

                f_num_contacts.close()
                f_cum_hydro.close()
                f_AA_factors.close()
                f_charge.close()
                f_OHrxnConst.close()
                f_SS.close()
                f_areaSAS.close()
                f_reactivity.close()
        		
                #Surface Methionines and Cysteines, S atom of methionine and cysteine gets oxidized
        	
                #Must select residues and give them names according to convention	
                rc("select :met@SD|:cys@SG")	
                MC_atoms = selection.currentAtoms()	
        	
                MC_atomsList = []
                MorC_dict = dict({})
        	
                for a in MC_atoms:	
                        residue = a.residue	
                        type = residue.type	
                        resNum = str(residue)	
                        resNum = re.sub(r"#.+ .+ ", "", resNum)	
        	
                        MC_resDict[resNum] = str(aaDict[type])	
                        name = a.name	
                        altLoc = a.altLoc	
                        if altLoc == "":	
                                string = str(residue) + " " + str(aaDict[type]) + " " + str(name) + altLoc	
                        else:	
                                string = str(residue) + " " + str(aaDict[type]) + " " + str(name) + '.' + altLoc	
                        string = re.sub(r" MET ", ":", string)	
                        string = re.sub(r" CYS ", ":", string)	
        	
                        MC_atomsList.append(string)	
                        MorC_dict[string] = str(aaDict[type])
        	
        	
                MC_atomsList = sorted(MC_atomsList) #this doesn't seem to be working correctly      	
        	
                ##COMPUTE DEPTH OF EACH ATOM (i.e. minimum distance to surface)
                if len(MC_atomsList) == 0:	
                	surface_MCs = []
        	        surface_Ms = 0
        	        surface_Cs = 0
                if len(MC_atomsList) > 0:	
        	        MC_depthDict = dict({}) #add methionine and cysteine in seperate list
        	        MC_minDepth = 1000000 # find min depth for MCs
        	        rc("select #0")
        	        rc("~select #0:@") #deselecting the atoms, so we just get the surface
        	#        rc("measure distance :arg@cd|:lys@ce|:mly@ce|:kcx@ce|:pro@cd|:thr@cb|:trp@cd1 selection multiple true show true")
        	        rc("measure distance :met@S=|:cys@S= selection multiple true show true") #from surface to all of these atoms
        	        r_mc = dialogs.find('reply')
        	        MC_text = r_mc.text.get('1.0', 'end')
        	        #r_mc.Clear()
        	        f_temp2 = open('./chimeraVOiD_temp2_' + fn + '.txt', 'w')
        	        f_temp2.write(MC_text)
        	        f_temp2.close()
        	
        	        with open('./chimeraVOiD_temp2_' + fn + '.txt','r') as f_temp2:
        	                for line in f_temp2:
        	                        line = line.rstrip()
        	                        if line[:21] == 'minimum distance from':
        	                                line = line.replace('minimum distance from ', '')
        	                                line = re.sub(r" to #0:\? = ", "\t", line)
        	                                resNum = line
        	                                resNum = re.sub(r"#.+:", "", resNum)
        	                                resNum = re.sub(r"@.+", "", resNum)
        	                                line = re.sub(r'@', (" " + MC_resDict[resNum] + " "), line)
        	                                string = re.sub(r"\t.+", "", line)
        	                                depth = re.sub(r".+\t", "", line)
        	                                MC_depthDict[string] = depth
        	                                if float(depth) < MC_minDepth:
        	                                        MC_minDepth = float(depth)
        	        os.remove('./chimeraVOiD_temp2_' + fn + '.txt')
        	
        	        surface_MCs =[]
        	        surface_Ms = 0
        	        surface_Cs = 0
        	        for a1 in MC_atomsList:
        	                atom1 = str(a1)
        	                MC_depth = MC_depthDict[atom1]
        	                if float(MC_depth) == MC_minDepth and float(MC_depth) <= 1.88:
        	                        f_surfaceMC.write(atom1)
        	                        f_surfaceMC.write("\t")
        	                        f_surfaceMC.write(str(MC_depth))
        	                        f_surfaceMC.write("\n")
        	
        	                        surface_MCs.append(atom1)

        	                        if MorC_dict[atom1] == 'M':
        	                                surface_Ms += 1
        	                        if MorC_dict[atom1] == 'C':
        	                                surface_Cs += 1
        	
        	        f_surfaceMC.close()
        	
                MC_contacts_dict = dict({})	
                MC_num_contacts_dict = dict({})	
                M_num_contacts_dict = dict({})	
                C_num_contacts_dict = dict({})	
                for a1 in atoms:	
                        residue = a1.residue	
                        type = residue.type	
                        resNum = str(residue)	
                        resNum = re.sub(r"#.+ .+ ", "", resNum)	
        	
                        name = a1.name	
                        altLoc = a1.altLoc	
                        if altLoc == "":	
                                string = str(residue) + " " + str(aaDict[type]) + " " + str(name) + altLoc	
                        else:	
                                string = str(residue) + " " + str(aaDict[type]) + " " + str(name) + '.' + altLoc	
                        string = re.sub(r" ARG ", ":", string)	
                        string = re.sub(r" LYS ", ":", string)	
                        string = re.sub(r" MLY ", ":", string)	
                        string = re.sub(r" KCX ", ":", string)	
                        string = re.sub(r" PRO ", ":", string)	
                        string = re.sub(r" THR ", ":", string)	
                        #string = re.sub(r" TRP ", ":", string)
                        atom_string = string
        	
        	
                        MC_contacts_count = 0	
                        M_contacts_count = 0	
                        C_contacts_count = 0	
                        MC_contacts_list = []	
        	
                        for b1 in MC_atoms:	
                                distance = chimera.distance(a1.xformCoord(), b1.xformCoord())	
                                residue = b1.residue	
                                type = residue.type	
                                resNum = str(residue)	
                                resNum = re.sub(r"#.+ .+ ", "", resNum)	
        	
                                MC_resDict[resNum] = str(aaDict[type])	
                                name = b1.name	
                                altLoc = b1.altLoc	
                                if altLoc == "":	
                                        string = str(residue) + " " + str(aaDict[type]) + " " + str(name) + altLoc #what exactly is the name	
                                else:	
                                        string = str(residue) + " " + str(aaDict[type]) + " " + str(name) + '.' + altLoc	
        	
                                string = re.sub(r" MET ", ":", string)	
                                MC_string = re.sub(r" CYS ", ":", string)	
                                        	
                                if distance <= radius and MC_string in surface_MCs:
                                        #print distance	
                                        MC_contacts_count += 1	
                                        MC_contacts_list.append(b1)	
        	
                                        print "MC_string:" + MC_string	
        	
                                        #specify if cysteine or methionine	
                                        if re.split(' ',MC_string)[1] == 'M':	
                                                M_contacts_count += 1	
                                        if re.split(' ',MC_string)[1] == 'C':	
                                                C_contacts_count += 1	
                                        #specify if cysteine or methionine	
        	
                        MC_num_contacts_dict[atom_string] = MC_contacts_count	
                        M_num_contacts_dict[atom_string] = M_contacts_count	
                        C_num_contacts_dict[atom_string] = C_contacts_count	
                        MC_contacts_dict[atom_string] = MC_contacts_list	
        	
                for a1 in atomsList:	
                        atom1 = str(a1)	
                        MC_num_contacts = MC_num_contacts_dict[atom1]	
                        M_num_contacts = M_num_contacts_dict[atom1]	
                        C_num_contacts = C_num_contacts_dict[atom1]	
        	
                        f_MC_contacts.write(atom1)	
                        f_MC_contacts.write("\t")	
                        f_MC_contacts.write(str(MC_num_contacts))	
                        f_MC_contacts.write("\t")	
                        f_MC_contacts.write(str(M_num_contacts))	
                        f_MC_contacts.write("\t")	
                        f_MC_contacts.write(str(C_num_contacts))	
                        f_MC_contacts.write("\t")	
                        f_MC_contacts.write(str(len(surface_MCs)))	
                        f_MC_contacts.write("\t")	
                        f_MC_contacts.write(str(surface_Ms))	
                        f_MC_contacts.write("\t")	
                        f_MC_contacts.write(str(surface_Cs))	
                        f_MC_contacts.write("\n")	
        	
                f_MC_contacts.close()	
                if len(MC_atomsList) > 0:	
                	f_temp2.close()
        	
        	
                #throwerror()	
        	

        #	# FOR DEBUGGING, OTHERWISE COMMENT OUT	
        #	f_debug.close()

        	
                rc("close all")	

                # Add success to log
                log = open(logfile, 'a')
                log.write("Success: " + fn + "\n")
                log.close()
#rc("stop now")	
	

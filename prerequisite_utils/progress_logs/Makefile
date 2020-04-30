all: total_fragments_count total_proteins_count

#### Step 0: Update all completed pdb files
update_confirmation:
	bash /home/js741/backup_scripts/copy_model3.sh >update_confirmation

##### Step 1: find all finished pdb files, put paths in a text file

find_pdb_all: update_confirmation
	find ../human/ -type f -name *.pdb > find_pdb_all

##### Step 2: find all runscripts in "Final_Push" for I-TASSER and QUARK

find_runscripts_finalpush_itasser:
	find /n/groups/drad/I-TASSER4.3/BatchJobs/Final_Push/ -maxdepth 4 -mindepth 2 -type f -name "*.sh" >find_runscripts_finalpush_itasser


find_runscripts_finalpush_quark:
	find /n/groups/drad/QUARKmod/BatchJobs/Final_Push/ -maxdepth 4 -mindepth 2 -type f -name "*.sh" >find_runscripts_finalpush_quark

##### Step 3: Find only model1.pdb lines from the finished pdb files
find_all_pdb_itasser_model1: find_pdb_all
	grep I-TASSER find_pdb_all | grep model1.pdb > find_all_pdb_itasser_model1

find_all_pdb_quark_model1: find_pdb_all
	grep QUARK find_pdb_all | grep model1.pdb > find_all_pdb_quark_model1


##### Step 4: Get all unique proteins and fragments that have been completed. 

# Cut out unique proteins and fragments from model1 files
pdb_all_itasser_model1_fragments: find_all_pdb_itasser_model1
	cut -d"/" -f4 find_all_pdb_itasser_model1 | sort | uniq >pdb_all_itasser_model1_fragments

pdb_all_itasser_model1_proteins: find_all_pdb_itasser_model1
	cut -d"/" -f3 find_all_pdb_itasser_model1 | sort | uniq >pdb_all_itasser_model1_proteins

pdb_all_quark_model1_fragments: find_all_pdb_quark_model1
	cut -d"/" -f4 find_all_pdb_quark_model1 | sort | uniq >pdb_all_quark_model1_fragments

pdb_all_quark_model1_proteins: find_all_pdb_quark_model1
	cut -d"/" -f3 find_all_pdb_quark_model1 | sort | uniq >pdb_all_quark_model1_proteins

# Combine counts from QUARK and I-TASSER to get overall counts
pdb_all_model1_fragments: pdb_all_quark_model1_fragments pdb_all_itasser_model1_fragments
	cat pdb_all_itasser_model1_fragments pdb_all_quark_model1_fragments | sort | uniq >pdb_all_model1_fragments

pdb_all_model1_proteins: pdb_all_itasser_model1_proteins pdb_all_quark_model1_proteins
	cat pdb_all_itasser_model1_proteins pdb_all_quark_model1_proteins | sort | uniq >pdb_all_model1_proteins

##### Step 5: Get all unique proteins and fragments that have been scripts in the works. 

# Repeat the process, this time on scripts in the finalpush files
finalpush_all_itasser_fragments: find_runscripts_finalpush_itasser
	cut -d"/" -f9 find_runscripts_finalpush_itasser | sed 's/.sh//' | sort | uniq >finalpush_all_itasser_fragments

finalpush_all_itasser_proteins: find_runscripts_finalpush_itasser
	cut -d"/" -f9 find_runscripts_finalpush_itasser | sed 's/.sh//' | cut -d'_' -f1 | sort | uniq >finalpush_all_itasser_proteins

finalpush_all_quark_fragments: find_runscripts_finalpush_quark
	cut -d"/" -f9 find_runscripts_finalpush_quark | sed 's/.sh//' | sort | uniq >finalpush_all_quark_fragments

finalpush_all_quark_proteins: find_runscripts_finalpush_quark
	cut -d"/" -f9 find_runscripts_finalpush_quark | sed 's/.sh//' | cut -d'_' -f1 | sort | uniq >finalpush_all_quark_proteins

# Combine QUARK and I-TASSER to get overall counts from scripts
finalpush_all_fragments: finalpush_all_quark_fragments finalpush_all_itasser_fragments
	cat finalpush_all_quark_fragments finalpush_all_itasser_fragments | sort | uniq >finalpush_all_fragments

finalpush_all_proteins: finalpush_all_quark_proteins finalpush_all_itasser_proteins
	cat finalpush_all_quark_proteins finalpush_all_itasser_proteins | sort | uniq >finalpush_all_proteins

##### Step 6: Get a combined list of all unique proteins and fragments that are either completed or queued up

total_fragments: finalpush_all_fragments pdb_all_model1_fragments
	cat finalpush_all_fragments pdb_all_model1_fragments | sort | uniq >total_fragments

total_proteins: finalpush_all_proteins pdb_all_model1_proteins
	cat finalpush_all_proteins pdb_all_model1_proteins | sort | uniq >total_proteins

##### Step 7: Make a file with fragment and protein counts
total_proteins_count: total_proteins
	bash get_total_protein_counts.sh

total_fragments_count: total_fragments
	bash get_total_fragment_counts.sh

clean:
	bash clean.sh
	

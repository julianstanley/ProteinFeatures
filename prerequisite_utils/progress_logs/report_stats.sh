itasser_fragments=$(wc -l < pdb_all_itasser_model1_fragments)
itasser_proteins=$(wc -l < pdb_all_itasser_model1_proteins)
quark_fragments=$(wc -l < pdb_all_quark_model1_fragments)
quark_proteins=$(wc -l < pdb_all_quark_model1_proteins)

echo "There are ${itasser_fragments} unique fragments from I-TASSER completed, from ${itasser_proteins} proteins."
echo "There are ${quark_fragments} unique fragments from QUARK completed, from ${quark_proteins} proteins."

total_fragments=$(wc -l < pdb_all_model1_fragments)
total_proteins=$(wc -l < pdb_all_model1_proteins)

echo "There are ${total_fragments} unique fragments in total completed, from ${total_proteins} proteins. Those proteins and fragments are in: 'pdb_all_model1_proteins' and 'pdb_all_model1_fragments', respectively."

# Find protein and fragment counts for QUARK and I-TASSER, this time with the runscript data
itasser_fragments_scripts=$(wc -l < finalpush_all_itasser_fragments)
itasser_proteins_scripts=$(wc -l < finalpush_all_itasser_proteins)
quark_fragments_scripts=$(wc -l < finalpush_all_quark_fragments)
quark_proteins_scripts=$(wc -l < finalpush_all_quark_proteins)

echo "There are ${itasser_fragments_scripts} unique fragments from I-TASSER with runscripts, from ${itasser_proteins_scripts} proteins."
echo "There are ${quark_fragments_scripts} unique fragments from QUARK with runscripts, from ${quark_proteins_scripts} proteins."

total_fragments_scripts=$(wc -l < finalpush_all_fragments)
total_proteins_scripts=$(wc -l < finalpush_all_proteins)

echo "There are ${total_fragments_scripts} unique fragments in total with runscripts, from ${total_proteins_scripts} proteins. Those proteins and fragments are in: 'finalpush_all_proteins' and 'finalpush_all_fragments', respectively."

all_fragments_count=$(wc -l <total_fragments)
all_proteins_count=$(wc -l <total_proteins)

echo "There are ${all_fragments_count} unique fragments that have been completed or queued up, from ${all_proteins_count} total proteins. Those proteins and fragments are in: "total_proteins" an "total_fragments", respectively."
cat "$1.csv" | awk '{split($0,a,","); print a[1]"\n"a[2]}' > "$1.fasta"

import argparse
import pandas as pd
import csv

doorfile = "door_operon_J2315.tabular"
csv_file = "results.csv"
fasta_file = "results.fasta"
annotation_file = "k52.gb"
start_operons = []
start_operon_ids = []
bad_operons = []
ranges53 = []
ranges35 = []
discard_tags = False

parser = argparse.ArgumentParser(description='Annotate and filter the results of sequence-finder.py')
parser.add_argument('doorfile', metavar='doorfile', type=str, help='door file with operon information')
parser.add_argument('csv_file', metavar='csv_file', type=str, help='resulting csv from sequence-finder.py')
parser.add_argument('fasta_file', metavar='fasta_file', type=str, help='resulting fasta file from sequence-finder.py')
parser.add_argument('annotation_file', metavar='annotation_file', type=str, help='genebank file with locus tags for data')
parser.add_argument('discard_tags', metavar='discard_tags', type=bool, help='discard locus tags from filter file?')

try:
    args = parser.parse_args()
    csv_file = args.csv_file
    doorfile = args.doorfile
    fasta_file = args.fasta_file
    annotation_file = args.annotation_file
    discard_tags = args.discard_tags
except:
    print()
    print("One or more file names not provided, attempting to use default filenames.")

    '''
	Goals:
		1. Add locus tag information to each gene based on position.
		2. Keep only rows which:
			- Have a locus tag present in the provided .csv
				or
			- Have a locus tag matching the first gene in an operon
			from the door file
    '''

# Step 1. genebank file
locus_tags_pos = {}
locus_tags_neg = {}
locus_tag_conversion = {}

def parse_gb():
    with open(annotation_file) as f:
        overlaps = 0
        for line in f:
            if 'gene' in line and '..' in line:
                strand = '+'
                start = line.strip().split()[1].split('..')[0]
                end = line.strip().split()[1].split('..')[1]

                if 'complement' in start:
                    start = start[11:]
                    end = end[:-1]
                    strand = '-'

                if '<' in start:
                    start = start[1:]
                if '>' in end:
                    end = end[1:]

                start = int(start)
                end = int(end)

                while '/old_locus_tag' not in line and '/locus_tag' not in line:
                    line = f.readline()

                # Sometimes there is no old_locus_tag, in which case we will hit the locus tag first
                if '/old_locus_tag' in line:
                    old_locus_tag = line.strip().split('=')[1][1:-1]
                else:
                    old_locus_tag = 'UNK'

                # If there was no old_locus_tag then this will hit immediately
                while '/locus_tag' not in line:
                    line = f.readline()
                locus_tag = line.strip().split('=')[1][1:-1]

                locus_tag_conversion[locus_tag] = old_locus_tag
                locus_tag_conversion[old_locus_tag] = locus_tag

                if strand == '+':
                    for pos in range(start-50, end):
                        if pos in locus_tags_pos:
                            locus_tags_pos[pos].append([locus_tag, start])
                        else:
                            locus_tags_pos[pos] = [[locus_tag, start]]
                else:
                    for pos in range(start, end+50):
                        if pos in locus_tags_neg:
                            locus_tags_neg[pos].append([locus_tag, end])
                        else:
                            locus_tags_neg[pos] = [[locus_tag, end]]

# Step 2. use door file to create a dictionary converting genes of operons to the start gene of the operon
# For any gene not in an operon, this dictionary will by default map a gene to itself.
operon_start_gene = {}

def parse_door():
    # initially map every gene to itself
    for gene_list in locus_tags_pos:
        gene_list = locus_tags_pos[gene_list]
        for gene in gene_list:
            operon_start_gene[gene[0]] = gene[0]
    for gene_list in locus_tags_neg:
        gene_list = locus_tags_neg[gene_list]
        for gene in gene_list:
            operon_start_gene[gene[0]] = gene[0]

    with open(doorfile) as f:
        reader = csv.reader(f, delimiter='\t')
        all_operons = {}
        prev_operon = -1
        operon_genes = []
        for row in reader:
            if len(row) > 0 and row[3] != "Start":
                operon_id = row[0]
                old_locus_tag = row[2]
                strand = row[5]


                # Add proper locus tag to operon
                if old_locus_tag in locus_tag_conversion:
                    locus_tag = locus_tag_conversion[old_locus_tag]
                else:
                    locus_tag = 'UNK'
                
                # on the positive strand, the first gene in an operon is our start,
                # for the negative strand, the last gene is the start.
                if operon_id != prev_operon:
                    if len(operon_genes) > 0:
                        if operon_genes[0]['strand'] == '+':
                            for gene in operon_genes:
                                operon_start_gene[gene['locus_tag']] = operon_genes[0]['locus_tag']
                        else:
                            for gene in operon_genes:
                                operon_start_gene[gene['locus_tag']] = operon_genes[-1]['locus_tag']
                    operon_genes = []

                prev_operon = operon_id
                operon_genes.append({'strand': strand, 'locus_tag': locus_tag})

        # repeat this once for the last operon in the file....
        if len(operon_genes) > 0:
            if operon_genes[0]['strand'] == '+':
                for gene in operon_genes:
                    operon_start_gene[gene['locus_tag']] = operon_genes[0]['locus_tag']
            else:
                for gene in operon_genes:
                    operon_start_gene[gene['locus_tag']] = operon_genes[-1]['locus_tag']


# Step 3. get tags to filter from locus_tags_to_filter.csv
def get_filter_tags():
    found_tags = []

    try:
        with open("locus_tags_to_filter.csv") as f:
            filter_df = pd.read_csv(f)
            tags = filter_df['Locus tag']
            new_tags = []

            '''
            Here we are making sure that if a locus tag from a later gene in an operon
            is present in the .csv file, then we add the first gene from that operon to
            keep.
            '''


            # Since any gene not in an operon marks itself as its own start gene, this will
            # include both the first gene of any operon and any gene which is not part of an
            # operon.
            for tag in tags:
                if tag in operon_start_gene:
                    found_tags.append(operon_start_gene[tag])
                elif str(tag) != "nan":
                    print("[!] Not found in genebank file: ", tag)

            found_tags = set(found_tags)

    except:
        print("Either <locus_tags_to_filter.csv> was not provided, or it was missing column: <Locus tag>")
        found_tags = []

    return found_tags

# ==================================================

print("Parsing genebank file...")
parse_gb()
print("Parsing door file...")
parse_door()
print("Parsing inclusion csv...")
filter_tags = get_filter_tags()

print("Filtering and saving .csv")
starts = {} # Used to filter the .fasta more easily
with open("results.csv") as f:
    outrows = []
    reader = csv.reader(f)
    missing = 0
    mislabeled = 0
    rows = []

    # add locus tags to all rows
    header = None
    for row in reader:
        if header is None:
            header = row
        else:
            if row[2] == '+' and int(row[0]) in locus_tags_pos:
                row.append(locus_tags_pos[int(row[0])])
                rows.append(row)
            elif row[2] == '-' and int(row[0]) in locus_tags_neg:
                row.append(locus_tags_neg[int(row[0])])
                rows.append(row)
            else:
                missing += 1
                print(row)

    final_rows = []

    for row in rows:
        for tag in row[3]:
           if (tag[0] in filter_tags and not discard_tags) or (tag[0] not in filter_tags and discard_tags):
               row[3] = tag
               final_rows.append(row)
               starts[int(row[0])] = tag[0]
               break

    with open("results-annotated.csv", 'w+') as f:
        wr = csv.writer(f)
        wr.writerow(header + ['locus_tag_and_gene_start'])
        for row in final_rows:
            wr.writerow(row)

    print(len(rows), " initial rows.")
    print(len(final_rows), " final rows.")

# -----------------------------------------

print("Filtering and saving .fasta")
with open("results.fasta") as f:
    lines = []
    for line in f:
        lines.append(line)

    index = 0
    final_lines = []

    for i in range(len(lines))[::2]:
        line = lines[i]
        start = int(line[1:-2])
        if start in starts:
            if '+' in line:
                line = line[:-1] + '|' + starts[start] + '\n'
            elif '-' in line:
                line = line[:-1] + '|' + starts[start] + '\n'
            final_lines.append(line)
            final_lines.append(lines[i+1])

    with open("results-annotated.fasta", 'w+') as f:
        f.writelines(final_lines)

print("All results saved to <results-annotated.csv> and <results-annotated.fasta>")



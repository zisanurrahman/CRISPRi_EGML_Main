'''
    closest_results.py

    Filters the input .csv file to contain only one result per locus tag,
    such that the chosen result is the one closest to the ATG.
'''

def process_csv():
    infile = "results-annotated.csv"
    outfile = "results-nearest.csv"
    cols = []
    final_dict = {}
    with open(infile) as f:
        for line in f:
            if len(cols) < 1:
                cols = line.split(',')
                cols = cols + ['gene_start']
            else:
                vals = line.split(',')
                vals[-1] = vals[-1][1:-3]
                vals[-2] = vals[-2][3:-1]
                locus_tag = vals[-2]
                locus_tag_start = vals[-1]

                attrs = {} 
                distance = int(locus_tag_start) - int(vals[0])

                if distance < 0:
                    distance = -distance

                for i in range(len(cols)):
                    attrs[cols[i]] = vals[i]

                if locus_tag in final_dict:
                    cur_distance = float(final_dict[locus_tag][cols[-1]]) - float(final_dict[locus_tag][cols[0]])
                    if cur_distance < 0:
                        cur_distance = -cur_distance

                    if distance < cur_distance:
                        final_dict[locus_tag] = attrs
                else:
                    final_dict[locus_tag] = attrs

    with open(outfile, 'w+') as f:
        header = ""
        lines = ""
        header = header[:-1]
        keys = [key for key in final_dict]
        entry = final_dict[keys[0]]
        keys = [key for key in entry]
        keys[-2] = "locus_tag"
        keys[-1] = "gene_start" 
        for key in keys:
            key = key.strip()
            header = header + key + ","

        header = header[:-1]
        print(header)
        for key in final_dict:
            entry = final_dict[key]
            line = ""
            for col in cols:
                line = line + entry[col] + ","
            line = line[:-1]
            line = line + '\n'
            lines = lines + line
        
        f.write(header)
        f.write("\n")
        f.writelines(lines)
    return final_dict

def process_fasta(key_dict):
    infile = "results-annotated.fasta"
    outfile = "results-nearest.fasta"
    out = ""
    with open( infile ) as f:
        alt_line = False
        add_line = False
        for line in f:
            if not alt_line:
                locus_tag = line.split('|')[1][:-1]
                position = line.split('|')[0][1:-1]
                alt_line = True

                if key_dict[locus_tag]['start'] == position:
                    out = out + line
                    add_line = True
            else:
                if add_line:
                    out = out + line
                alt_line = False
                add_line = False

    with open(outfile, 'w+') as f:
        f.write(out)


final_csv_dict = process_csv()
process_fasta(final_csv_dict)



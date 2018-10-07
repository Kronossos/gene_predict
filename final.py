import collections
import re
import sys

# Function to compare whether a gene with a given begin and end has appeared before.
def compare(current_data, last_gene, data_dict, limit=800):
    for compared_data in data_dict:
        if compared_data == current_data: continue
        for first_number in data_dict[compared_data]:
            if abs(data_dict[compared_data][first_number]["bg"] - data_dict[current_data][last_gene]["bg"]) > limit: continue

            if abs(data_dict[compared_data][first_number]["end"] - data_dict[current_data][last_gene]["end"]) > limit: continue

            data_dict[compared_data][first_number]["frequency"] += 1
            data_dict[current_data][last_gene]["frequency"] += 1

            data_dict[compared_data][first_number]["Other match"] += current_data + "_" + last_gene + " "
            data_dict[current_data][last_gene]["Other match"] += compared_data + "_" + first_number + " "

# Function to read data from given files.
def create_dict(file, data_dict):
    if file not in data_dict:
        data_dict[file] = collections.OrderedDict()

    last_gene = ""
    with open("temp/clear/" + str(file)) as data:
        for line in data:

            gene, sign, bg, end = line.replace("\n", "").split("\t")

            if gene not in data_dict[file]:
                data_dict[file][gene] = {"sign": sign, "bg": int(bg), "end": 0, "Other match": "", "frequency": 0}

            if gene != last_gene and last_gene:
                data_dict[file][last_gene]["end"] = int(last_end)
                compare(file, last_gene, data_dict)

            last_gene = gene
            last_end = end
        else:
            data_dict[file][last_gene]["end"] = int(last_end)
            compare(file, last_gene, data_dict)

# Fuction to load seqence from file into memory.
def read_seq(seq_file):
    seq = ""
    with open(str(seq_file)) as data_seqence:
        for line in data_seqence:
            if re.search(">", line): continue
            seq += line.replace("\n", "")
    return seq

# Write output in csv. For more see "usage" in "run" file.
def write(data_dict):
    with open("result/result.csv", "w") as result_file:
        for sheet_name in data_dict:
            result_file.write(sheet_name + "\t" + "+\-" + "\t" + "START" + "\t" + "STOP" + "\t" + "Other match" + "\n")
            for gene_number in data_dict[sheet_name]:
                result_file.write(
                    gene_number + "\t" + data_dict[sheet_name][gene_number]["sign"] + "\t" + str(data_dict[sheet_name][gene_number]["bg"]) + "\t" + str(
                        data_dict[sheet_name][gene_number]["end"]) + "\t" + str(data_dict[sheet_name][gene_number]["Other match"]) + "\n")

# Save the most common prediction in separate files. This step is necessary to use NCBI tool.
def select_predictions(data_dict):
    result = 1
    for sheet_name in data_dict:
        for gene_number in data_dict[sheet_name]:
            if data_dict[sheet_name][gene_number]["frequency"] == 2:
                similar = data_dict[sheet_name][gene_number]["Other match"][:-1].split(" ")
                for element in similar:
                    tool, tool_number = element.split("_")
                    data_dict[tool][tool_number]["frequency"] = 0

                with open("temp/result/" + str(result), "w") as result_file:
                    result_file.write(seq[data_dict[sheet_name][gene_number]["bg"]:data_dict[sheet_name][gene_number]["end"]])
                result += 1


sequence_file = sys.argv[1]
seq = read_seq(sequence_file)

data = collections.OrderedDict()

# Tools names (same name as file).
for name in ["Geneid", "Genescan", "Genmark"]:
    create_dict(name, data)

write(data)
select_predictions(data)





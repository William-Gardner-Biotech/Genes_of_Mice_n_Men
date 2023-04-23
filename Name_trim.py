import re

# Blast results path goes here
results_location = "~/$YOURPATH$/mouse_vs_human_blast_results.tab"

results = open(results_location, "r")
out_list = ["Mouse Gene, Human Gene"]
for line in results:
    splits = line.split("\t")
    mouse = (splits[0])
    subgroups_mouse = re.search("((\w*)$)", mouse)
    mouse = (subgroups_mouse.group(0))
    human = (splits[1])
    subgroups_human = re.search("((\w*)$)", human)
    human = (subgroups_human.group(0))
    addition = "\n"+mouse+","+human
    out_list.append(addition)

outfile = open("~$YOURPATH$/Mouse_Human.csv", "w")

for i in out_list:
    outfile.write(i)

outfile.close()

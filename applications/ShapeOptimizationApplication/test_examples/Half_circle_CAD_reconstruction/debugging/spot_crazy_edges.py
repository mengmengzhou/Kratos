# compute edges 3418
#       process mdpa
import csv
input_file = "optimized_shell_3418_elements.mdpa"
edges = []
with open(input_file) as file:
    reader = csv.reader(file, delimiter=" ")
    while not 'Conditions' in reader.__next__():
        pass
    for row in reader:
        if 'End' in row:
            break
        values = [x for x in row if x!= '']
        nodes = [values[i] for i in [2,3,4]]
        nodes.sort()
        [n1,n2,n3] = nodes
        e1 = [n1,n2]
        e2 = [n1,n3]
        e3 = [n2,n3]
        if e1 not in edges:
            edges.append(e1)
        if e2 not in edges:
            edges.append(e2)
        if e3 not in edges:
            edges.append(e3)
for x in edges:
    print(x)

# compute corresponding lengths
#       init a mapper with updated json
# mapper = ...

# mapper.compute_length

# get top 20 lengths

# write to file xyz of corresponding FE-nodes
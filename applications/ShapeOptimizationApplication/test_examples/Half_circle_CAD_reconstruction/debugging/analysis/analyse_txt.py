from math import sqrt

class Point:
    def __init__(self, string_list):
        string_list = [float(x) for x in string_list]
        if len(string_list) == 3:
            self.type = 0
            [self.x, \
             self.y, \
             self.z] = string_list

        elif len(string_list) == 10:
            self.type = 1
            [self.patch, \
             self.u, \
             self.v, \
             self.u_min, \
             self.u_max, \
             self.v_min, \
             self.v_max, \
             self.x, \
             self.y, \
             self.z] = string_list
            self.patch = int(self.patch)
        elif len(string_list) == 6:
            self.type = 2
            [self.patch, \
             self.u, \
             self.v, \
             self.x, \
             self.y, \
             self.z] = string_list            
            self.patch = int(self.patch)
        else:
            print("len of string is", len(string_list))
    def distance_from(self, other):
        return sqrt( (self.x - other.x)**2 +
                     (self.y - other.y)**2 +
                     (self.z - other.z)**2
                     )

def check_patch(FE_node, neighbour, newton_raphson):
    if bla < FE_node.x < bla and \
       bla < FE_node.y < bla and \
       bla < FE_node.z < bla:
       correct_patch = bla
       return [correct_patch == neighbour.patch, correct_patch == newton_raphson.patch]

def print_mapping_limits(fe_nodes_l, neighbours, newton_rap):
    assert len(fe_nodes_l) == len(neighbours) == len(newton_rap)
    for i in range(4):
            print("patch", i)
            print_patch_mapping_limits(i, fe_nodes_l, neighbours, newton_rap)

def print_patch_mapping_limits(num, fe_nodes_l, neighbours, newton_rap):
    patch = [(fe, ne, nr) \
               for fe, ne, nr in zip(fe_nodes_l, neighbours, newton_rap) \
               if ne.patch == num \
               and same_patch(ne,nr)]
    x_min_tuple = min(patch, key = lambda tup: tup[0].x)
    x_max_tuple = max(patch, key = lambda tup: tup[0].x)
    y_min_tuple = min(patch, key = lambda tup: tup[0].y)
    y_max_tuple = max(patch, key = lambda tup: tup[0].y)
    z_min_tuple = min(patch, key = lambda tup: tup[0].z)
    z_max_tuple = max(patch, key = lambda tup: tup[0].z)

    x_min = x_min_tuple[0].x
    x_max = x_max_tuple[0].x
    y_min = y_min_tuple[0].y
    y_max = y_max_tuple[0].y
    z_min = z_min_tuple[0].z
    z_max = z_max_tuple[0].z
    print(x_min, x_max, y_min, y_max, z_min,z_max)
    print()


def same_patch(neighbo, newton_rap):
    assert neighbo.patch == newton_rap.patch
    return True

def in_out_of_patch(fe_nodes_l, neighbours, newton_rap):
    inside = list()
    outside = list()
    for i, [fe, ne, nr] in enumerate(zip(fe_nodes_l, neighbours, newton_rap)):
        if ne.u_min < nr.u < ne.u_max and \
           ne.v_min < nr.v < ne.v_max:
           inside.append((fe,ne,nr))
        else:
           outside.append((fe,ne,nr))
    print("{0:.2f}% of the points are out of patch".format( len(outside)/(len(inside)+len(outside)) ) )
    return [inside, outside]

def bad_reconstructed(fe_nodes_l, neighbours, newton_rap, threshold):
    return [(fe, ne, nr) \
            for fe, ne, nr in zip(fe_nodes_l, neighbours, newton_rap) \
            if fe.distance_from(nr) > threshold]

def print_to_file(filename, tuple_list, id):
    print("Printing to file...")
    with open(filename, "w") as output_file:
        for [fe, ne, nr] in tuple_list:
            if id == 0:
                # print fe
                print(fe.x, fe.y, fe.z, file = output_file)
            elif id == 1:
                # print ne
                print(ne.x, ne.y, ne.z, file = output_file)
            elif id == 2:
                # print nr
                print(nr.x, nr.y, nr.z, file = output_file)
            else:
                raise("Improper use of function")

####################################################################################################################### 
import csv

fe_nodes_l = list()
neighbours = list()
newton_rap = list()
fe_new_lis = list()
nr_new_lis = list()

# read original data
input_file = "500x500.txt"
with open(input_file) as file:
    reader = csv.reader(file, delimiter=" ")
    # ignore two header lines
    reader.__next__()
    reader.__next__()
    for row in reader:
        fe_node = row[:3]
        neighbo = row[3:13]
        new_rap = row[13:]
        fe_nodes_l.append(Point(fe_node))
        neighbours.append(Point(neighbo))
        newton_rap.append(Point(new_rap))

# read updated data
input_file = "500x500_updated.txt"
with open(input_file) as file:
    reader = csv.reader(file, delimiter=" ")
    # ignore two header lines
    reader.__next__()
    reader.__next__()
    for row in reader:
        fe_node = row[:3]
        new_rap = row[13:]
        fe_new_lis.append(Point(fe_node))
        nr_new_lis.append(Point(new_rap))

####################################################################################################################### 

print_mapping_limits(fe_nodes_l, neighbours, newton_rap)
[inside, outside] = in_out_of_patch(fe_nodes_l,neighbours, newton_rap)
filename = "out_of_patch_points_nr.txt"
print_to_file(filename, outside, 2)

for t in  [0.001*x for x in range(11)]:
    br_list = bad_reconstructed(fe_new_lis,  neighbours, nr_new_lis, t)
    print(len(br_list), "points with threshold =", t)
    filename = "fe_points_bad_reconstructed_{0}.txt".format(t)
    print_to_file(filename, br_list, 0)
    filename = "nr_points_bad_reconstructed_{0}.txt".format(t)
    print_to_file(filename, br_list, 2)

print_to_file("fe_points.txt", list(zip(fe_new_lis, neighbours, nr_new_lis)), 0)
print_to_file("nr_points.txt", list(zip(fe_new_lis, neighbours, nr_new_lis)), 2)

distances = [fe.distance_from(nr) for fe, nr in zip(fe_new_lis, nr_new_lis)]
print("min:", min(distances))
print("max:", max(distances))
print("objective:", sum((x**2 for x in distances)))


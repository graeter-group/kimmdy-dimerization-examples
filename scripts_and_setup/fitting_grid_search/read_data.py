def read_angle_file(file_name_in):
    file = open(file_name_in)
    data = False
    data_angles = {}

    for line in file:
        if "@TYPE xy" in line and not data:
            data = True
            continue
        if data:
            line = [a for a in line.strip("\n").split(" ") if a]
            time = int(float(line[0]))  # in ps
            values = [float(a) for a in line[1:]]
            data_angles[time] = values
    return data_angles


def read_distance_file(file_name_in):
    file = open(file_name_in)
    data = False
    data_distances = {}

    for line in file:
        if "@TYPE xy" in line and not data:
            data = True
            continue
        if data:
            line = [a for a in line.strip("\n").split(" ") if a]
            time = int(float(line[0]))  # in ps
            residues = int((len(line)-1)/2 + 1)
            values = [(float(line[2*i-1]) + float(line[2*i])) / 2 for i in range(1, residues)]  # in nm
            data_distances[time] = values
    return data_distances

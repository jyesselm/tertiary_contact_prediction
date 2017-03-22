from rnamake import structure, util

def contact_map_from_pdbs(pdbs, res_seperation=8, max_dist=14, max=100000):
    if len(pdbs) == 0:
        raise ValueError("no pdbs supplied")

    s = structure.structure_from_pdb(pdbs[0])
    num_res = len(s.residues())
    matrix = []

    for i in range(num_res):
        matrix.append([0 for j in range(num_res)])

    for k, pdb in enumerate(pdbs):
        if k > max:
            break
        s = structure.structure_from_pdb(pdb)
        residues = s.residues()

        for i, r1 in enumerate(residues):
            try:
                beads = r1.get_beads()
                a1 = beads[-1]
            except:
                continue
            for j, r2 in enumerate(residues):
                if abs(i - j) < res_seperation:
                    continue
                if i >= j:
                    continue
                try:
                    beads = r2.get_beads()
                    a2 = beads[-1]
                except:
                    continue

                dist = util.distance(a1.center, a2.center)
                if dist < max_dist:
                    matrix[i][j] += 1
    return matrix

def write_contact_map(filename, matrix):
    num_res = len(matrix)
    f = open(filename, "w")
    for i in range(num_res):
        for j in range(num_res):
            f.write(str(matrix[i][j]) + " ")
        f.write("\n")
    f.close()

def read_contact_map(filename):
    f = open(filename)
    lines = f.readlines()
    f.close()

    matrix = []
    for l in lines:
        row = [float(x) for x in l.split()]
        matrix.append(row)

    return matrix
def interpolate(x1, x2, y1, y2, x):
    return y2 + (x - x2) * (y2 - y1) / (x2 - x1)

def map_positions(map, positions):
    """
    Iterate through genetic positions and interpolate coordinates in cM.
    """
    genetic_dist = []
    idx = 0
    n_row = len(map.index)
    for p in positions:
        while True:
            row = map.loc[idx]
            if p == row.Pos:
                genetic_dist.append(row.GeneticDist)
                break
            elif (p < row.Pos and idx != 0) or (idx + 1 == n_row):
                interp = interpolate(map.loc[idx-1, "Pos"], row.Pos, map.loc[idx-1, "GeneticDist"], row.GeneticDist, p)
                genetic_dist.append(interp)
                break
            elif p < row.Pos and idx == 0:
                interp = interpolate(row.Pos, map.loc[idx+1, "Pos"], row.GeneticDist, map.loc[idx+1, "GeneticDist"], p)
                genetic_dist.append(interp)
                break
            else:
                idx += 1
    return(genetic_dist)

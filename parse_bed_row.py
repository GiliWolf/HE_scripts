passed_es_map = {19: 40, 22: 40, 24: 40, 25: 41}

def parse_bed_row(read, passed_es_map):
    genomic_blocks = [(114507396, 114507422), (114667443, 114667457)]
    read_blocks = [(0, 25), (26, 39)]
    rows = []
    for i, block in enumerate(read_blocks):
        block_start, block_end = block[0], block[1]
        genomic_start = genomic_blocks[i][0]

        # find cluster: first es in the block, last es in the block, and number of es in the block
        first_fit_pos = None  # Initialize to None to check later
        last_fit_pos = None
        cluster_score = 0
        for pos in passed_es_map:
            if block_start <= pos <= block_end:
                cluster_score += 1  # Increment the sum for each matching position
                if first_fit_pos is None:  # If first_fit_pos hasn't been set yet
                    first_fit_pos = pos
                last_fit_pos = pos  # update last_fit_pos to the current match
        print(f"block ${i}: first_fit_pos = ${first_fit_pos}")
        print(f"block ${i}: last_fit_pos = ${last_fit_pos}")

        if cluster_score==0:
            continue
        rows.append({
            'chr': 1,
            'start': genomic_start + first_fit_pos - block_start,
            'end': genomic_start + last_fit_pos - block_start,
            'name': 1,
            'score': cluster_score,
            'strand':1
        })  
    return rows

rows = parse_bed_row(0, passed_es_map)
print(rows)
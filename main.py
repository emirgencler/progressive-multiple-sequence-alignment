# 230601001 
# Emir Gençler 
# Computer Engineering


# Read protein sequences from FASTA file
def load_sequences(fname):
    names ,sequences = [], []
    with open(fname, encoding="utf-8") as file:
        for raw_line in file:
            line = raw_line.strip()
            if line.startswith(">"):
                names.append(line[1:].strip())
                sequences.append("")
            elif line:
                sequences[-1] += line.upper()
    return names, sequences


def load_scores(fname):
    data = []
    with open(fname, encoding="utf-8") as file:
        for raw_line in file:
            line = raw_line.strip()
            if line and not line.startswith("#"):
                data.append(line.split())

    letters = data[0]
    original = {}
        
    for row in data[1:]:
        symbol = row[0]
        original[symbol] = []

        for value in row[1:]:
            original[symbol].append(int(value))

    # Transpose BLOSUM62 and add 2 to diagonal
    scores = {}

    for i in range(len(letters)):
        for j in range(len(letters)):
            value = original[letters[j]][i]

            if i == j:
                value += 2

            scores[letters[i], letters[j]] = value

    return scores


def make_table(rows, columns, initial):
    table = []
    for i in range(rows):
        table.append([initial] * columns)
    return table


def choose_move(diagonal, up, left):
    # Tie priority: diagonal, up, left.
    # D: diagonal, U: gap in second sequence, L: gap in first sequence
    best = max(diagonal, up, left)
    if diagonal == best:
        return best, "D"
    if up == best:
        return best, "U"
    return best, "L"

# Global alignment using dynamic programming
def align_two(first, second, scores, gap):
    n = len(first)
    m = len(second)
    dp = make_table(n + 1, m + 1, 0)
    back = make_table(n + 1, m + 1, "")

    for i in range(1, n + 1):
        dp[i][0], back[i][0] = i * gap, "U"
    for j in range(1, m + 1):
        dp[0][j], back[0][j] = j * gap, "L"

    # Find the best score for each position
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            d = dp[i-1][j-1] + scores[first[i-1], second[j-1]]
            u = dp[i-1][j] + gap
            l = dp[i][j-1] + gap
            dp[i][j], back[i][j] = choose_move(d, u, l)
    # Trace back to construct the aligned sequences
    a, b, i, j = [], [], n, m
    while i > 0 or j>0 :
        move = back[i][j]
        if move != "L":
            a.append(first[i-1])
            i -= 1
        else: a.append("-")
        if move != "U":
            b.append(second[j-1])
            j -= 1
        else: b.append("-")

    return dp[n][m], "".join(reversed(a)), "".join(reversed(b))

# Calculate pairwise similarity scores
def compare_pairs(names, sequences, scores, gap):
    k = len(names)
    #initialize 
    similarities = make_table(k, k, 0.0)
    for i in range(k):
        similarities[i][i] = 1.0

    results = []
    # Align each unique pair of sequences
    for i in range(k):
        for j in range(i + 1, k):
            score, first, second = align_two(sequences[i], sequences[j], scores, gap)
            matches = 0
            for a, b in zip(first, second):
                if a == b and a != "-":
                    matches += 1
            length = len(first)
            # Similarity is the fraction of identical positions
            similarity = matches / length
            similarities[i][j] = similarity
            similarities[j][i] = similarity
            results.append((i, j, score, length, matches, similarity, first, second))
    return similarities, results


def display_matrix(labels, values):
    # Display similarities with five decimal places.
    width = 9
    for name in labels:
        if len(name) + 2 > width:
            width = len(name) + 2
    print(" " * width, end="")
    for name in labels:
        print(f"{name:>{width}}", end="")
    print()

    for i in range(len(labels)):
        print(f"{labels[i]:<{width}}", end="")
        for value in values[i]: print(f"{value:>{width}.5f}", end="")
        print()


def format_tree(node):
    if isinstance(node, str):
        return node
    return "(" + format_tree(node[0]) + "," + format_tree(node[1]) + ")"


def cluster_similarity(first, second, original_matrix):
    total = 0
    for i in first:
        for j in second: total += original_matrix[i][j]
    # Average similarity between two clusters
    return total / (len(first) * len(second))

# Build the guide tree using UPGMA
def build_guide_tree(names, matrix):
    groups, history = [], []
    for i in range(len(names)): groups.append(((i,), names[i]))

    while len(groups) > 1:
        # Find the two most similar clusters
        best_key = None
        best_i, best_j = 0, 1

        for i in range(len(groups)):
            for j in range(i + 1, len(groups)):
                sim = cluster_similarity(groups[i][0], groups[j][0], matrix)
                key = (-sim, min(groups[i][0]), min(groups[j][0]))
                if best_key is None or key < best_key:
                    best_key, best_i, best_j = key, i, j

        a, b, merged, remaining = groups[best_i], groups[best_j], None, []
        merged = (tuple(sorted(a[0] + b[0])), (a[1], b[1]))

        remaining = []
        for i in range(len(groups)):
            if i != best_i and i != best_j: remaining.append(groups[i])
        remaining.append(merged)
        groups = sorted(remaining, key=lambda g: min(g[0]))

        labels = []
        for group in groups:
            members = []
            for i in group[0]: members.append(names[i])
            labels.append("{" + ",".join(members) + "}")
        # Update the similarity matrix after merging
        updated = make_table(len(groups), len(groups), 0.0)
        for i in range(len(groups)):
            for j in range(len(groups)):
                if i == j: updated[i][j] = 1.0
                else: updated[i][j] = cluster_similarity(
                    groups[i][0], groups[j][0], matrix)

        history.append((a[1], b[1], -best_key[0], labels, updated))

    return groups[0][1], history

# Calculate gap cost for a profile column
def score_column_gap(rows, column, gap):
    count = 0
    for row in rows:
        if row[column] != "-": count += 1
    return gap * count / len(rows)


def compare_columns(rows_a, i, rows_b, j, scores, gap):
    # Calculate average score between two profile columns
    total = 0
    for first in rows_a:
        for second in rows_b:
            a, b = first[i], second[j]
            if a == b == "-": value = 0
            elif a == "-" or b == "-": value = gap
            else: value = scores[a, b]
            total += value
    return total / (len(rows_a) * len(rows_b))

# Align two groups of sequences
def merge_alignments(group_a, group_b, scores, gap):
    a, b = list(group_a.values()), list(group_b.values())
    n, m = len(a[0]), len(b[0])
    # Initialize dp and traceback tables
    dp = make_table(n+1, m+1, 0.0)
    back = make_table(n+1, m+1, "")
    
    for i in range(1, n+1):
        dp[i][0], back[i][0] = dp[i-1][0]+score_column_gap(a,i-1,gap), "U"
    for j in range(1, m+1):
        dp[0][j], back[0][j] = dp[0][j-1]+score_column_gap(b,j-1,gap), "L"
    # Find the best alignment between two profiles
    for i in range(1, n+1):
        for j in range(1, m+1):
            d = dp[i-1][j-1]+compare_columns(a,i-1,b,j-1,scores,gap)
            u = dp[i-1][j]+score_column_gap(a,i-1,gap)
            l = dp[i][j-1]+score_column_gap(b,j-1,gap)
            dp[i][j], back[i][j] = choose_move(d,u,l)
    # Trace back to build the merged alignment
    result = {name: [] for name in list(group_a)+list(group_b)}

    # Gaps inserted into a group are added to all its sequences
    i, j = n, m
    while i > 0 or j > 0:
        for name, row in group_a.items():
            result[name].append("-" if back[i][j]=="L" else row[i-1])
        for name, row in group_b.items():
            result[name].append("-" if back[i][j]=="U" else row[j-1])
        move = back[i][j]
        if move != "L": i -= 1
        if move != "U": j -= 1
    for name in result: result[name] = "".join(reversed(result[name]))
    return result
# Merge sequences following the guide tree
def progressive_alignment(tree, sequences, scores, gap, history):
    # Return the sequence if the node is a leaf
    if isinstance(tree, str):
        return {tree: sequences[tree]}
    first = progressive_alignment(tree[0], sequences, scores, gap, history)
    second = progressive_alignment(tree[1], sequences, scores, gap, history)
    # Merge the two groups following the guide tree
    aligned = merge_alignments(first, second, scores, gap)
    history.append((format_tree(tree), aligned))
    return aligned

# Print alignments in blocks of 50 columns
def display_alignment(alignment, names, mark_conserved=False):
    present = []
    for name in names:
        if name in alignment:
            present.append(name)
    length = len(alignment[present[0]])
    total = 0
    for start in range(0, length, 50):
        stop = min(start + 50, length)
        print("  Columns", start + 1, "-", stop)
        for name in present:
            print(f"  {name:<7} {alignment[name][start:stop]}")
        # Mark positions shared by all sequences
        if mark_conserved:
            marks = ""
            for i in range(start, stop):
                first = alignment[present[0]][i]
                same = first != "-"
                for name in present:
                    if alignment[name][i] != first:
                        same = False
                if same:
                    marks += "*"
                    total += 1
                else:
                    marks += " "
            print(" " * 10 + marks)
        print()
    return total


def run_project():
    names, seqs = load_sequences("Input.txt")
    scores = load_scores("Blosum62.txt")
    gap = int(input("Enter gap penalty: "))
    if gap >= 0: raise ValueError("Gap must be negative")
    print("\nSEQUENCE LENGTHS")
    for i in range(len(names)): print(names[i], ":", len(seqs[i]))
    print("\nSTEP 1: PAIRWISE ALIGNMENTS")
    matrix, pairs = compare_pairs(names, seqs, scores, gap)
    for i, j, score, length, matches, sim, a, b in pairs:
        print(names[i], "-", names[j], "Score:", score, "Length:", length)
        print("Matches:", matches, "Similarity:", round(sim, 5))
    best = max(pairs, key=lambda p: p[5])
    a, b = names[best[0]], names[best[1]]
    print("\nMOST SIMILAR PAIR:", a, b)
    display_alignment({a: best[6], b: best[7]}, names)
    print("SIMILARITY MATRIX"); display_matrix(names, matrix)
    print("\nSTEP 2: UPGMA GUIDE TREE")
    tree, merges = build_guide_tree(names, matrix)
    for num, (left, right, sim, labels, updated) in enumerate(merges, 1):
        print("Merge", num, format_tree(left), "+", format_tree(right), round(sim, 5))
        display_matrix(labels, updated)
    print("GUIDE TREE:", format_tree(tree))
    print("\nSTEP 3: PROGRESSIVE MULTIPLE ALIGNMENT")
    seq_dict, steps = dict(zip(names, seqs)), []
    final = progressive_alignment(tree, seq_dict, scores, gap, steps)
    for num, (subtree, aligned) in enumerate(steps, 1):
        print("Alignment after merge", num, subtree)
        display_alignment(aligned, names)
    length = len(final[names[0]])
    for name in names:
        assert len(final[name]) == length
        assert final[name].replace("-", "") == seq_dict[name]
    print("\nFINAL MULTIPLE ALIGNMENT:")
    total = display_alignment(final, names, mark_conserved=True)
    print("Final length:", length, "Conserved columns:", total)
if __name__ == "__main__": 
    run_project()

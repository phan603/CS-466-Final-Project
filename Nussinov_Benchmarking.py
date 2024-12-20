import numpy as np
import matplotlib.pyplot as plt
import time

def nussinov(sequence, match_score=1, mismatch_score=0, min_loop_length=0):
    """
    Nussinov algorithm to predict RNA secondary structure using dynamic programming.
    """
    n = len(sequence)
    dp = np.zeros((n, n), dtype=int)

    # Filling DP table
    for l in range(1, n):  # l is the length of the subproblem
        for i in range(n - l):
            j = i + l
            if j - i > min_loop_length:  # Only consider pairing if j-i > min_loop_length
                match = (sequence[i], sequence[j])
                pair_score = match_score if is_complementary(match) else mismatch_score
                dp[i, j] = max(
                    dp[i + 1, j],  # i is unpaired
                    dp[i, j - 1],  # j is unpaired
                    dp[i + 1, j - 1] + pair_score if j - i > min_loop_length else 0,  # i, j paired
                    max(dp[i, k] + dp[k + 1, j] for k in range(i, j))  # bifurcation
                )
    return dp

def is_complementary(pair):
    """
    Check if a nucleotide pair is complementary (A-U, G-C)
    including wobble pairs (G-U).
    """
    complementary_pairs = {"AU", "UA", "GC", "CG", "GU", "UG"}
    return "".join(pair) in complementary_pairs

def get_memory_usage(dp_table):
    """
    Estimate the memory usage of the DP table in bytes.
    """
    return dp_table.nbytes

# Generate large RNA sequences
def generate_large_rna(length):
    """
    Generate a random RNA sequence of a given length.
    """
    return ''.join(np.random.choice(['A', 'U', 'G', 'C'], size=length))

# Test large RNA strings
sequence_lengths = [100, 500, 1000, 5000, 10000]  # Sizes for testing
times = []
memory_usages = []

for length in sequence_lengths:
    sequence = generate_large_rna(length)
    start_time = time.time()
    dp_table = nussinov(sequence)
    end_time = time.time()

    elapsed_time = end_time - start_time
    memory_usage = get_memory_usage(dp_table)

    times.append(elapsed_time * 1000)  # Convert to milliseconds
    memory_usages.append(memory_usage)

    print(f"Sequence Length: {length} | Time Taken: {elapsed_time*1000:.2f} ms | Memory: {memory_usage / 1024 / 1024:.2f} MB")

# Plotting results
# Length vs Time
plt.figure(figsize=(10, 5))
plt.plot(sequence_lengths, times, marker='o')
plt.title("Sequence Length vs Time Taken")
plt.xlabel("Sequence Length")
plt.ylabel("Time Taken (ms)")
plt.grid(True)
plt.show()
import random

NUM_VALUES = 10 ** 5
RANGE_MIN = 0
RANGE_MAX = 10 ** 10

OUTPUT_FILE = 'input1.txt'
nums = random.sample(range(RANGE_MIN, RANGE_MAX + 1), NUM_VALUES)

with open(OUTPUT_FILE, "w") as f:
    f.write(f"{NUM_VALUES}\n")
    for num in nums:
        f.write(f"{num}\n")

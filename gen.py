import random


TYPE = 3

if TYPE == 1:
    x = int(input())
    RANGE_MIN, RANGE_MAX, NUM_VALUES = 0, 10 ** 10, 10 ** 5
    nums = random.sample(range(RANGE_MIN, RANGE_MAX + 1), NUM_VALUES)
    OUTPUT_FILE = f'test_data/input_rand{x}.in'
elif TYPE == 2:
    x = int(input())
    length = 18 - x
    RANGE_MIN, RANGE_MAX, NUM_VALUES0 = 0, 2 ** length, min(10 ** 5, 2 ** length)

    samples = random.sample(range(RANGE_MIN, RANGE_MAX), NUM_VALUES0)
    bstrs = ["0" * (length - len(bin(sample)[2:])) + bin(sample)[2:] for sample in samples]
    for i, bstr in enumerate(bstrs):
        for j in range(34 - length):
            bstr += bstr[j]
        bstrs[i] = bstr
    nums = [int(bstr, 2) for bstr in bstrs]

    RANGE_MIN2, RANGE_MAX2, NUM_VALUES = 0, 10 ** 10, 10 ** 5
    st = set(nums)
    while len(st) < NUM_VALUES:
        st.add(random.randint(RANGE_MIN2, RANGE_MAX2 + 1))
    nums = list(st)
    OUTPUT_FILE = f'test_data/input_ps{x}.in'
elif TYPE == 3:
    x = int(input())
    NUM_VALUES = 10 ** 5
    length, no_blocks = 21, 300
    st = set()
    for i in range(no_blocks):
        val = random.randint(0, 2 ** 21 - 1)
        for j in range(10 ** 5 // no_blocks):
            val2 = random.randint(0, 2 ** 13 - 1)
            st.add(val * 2 ** 13 + val2)
    
    while len(st) < NUM_VALUES:
        st.add(random.randint(0, 10 ** 10 + 1))

    nums = list(st)
    OUTPUT_FILE = f'test_data/input_pref{x}.in'

with open(OUTPUT_FILE, "w") as f:
    f.write(f"{NUM_VALUES}\n")
    for num in nums:
        f.write(f"{num}\n")

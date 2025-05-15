import sys
import subprocess
import os
import math
from collections import defaultdict as dd

if len(sys.argv) < 3:
    print("Invalid number of arguments")
    sys.exit(1)

program_path = sys.argv[1]
second_arg = sys.argv[2]

if os.path.isdir(second_arg):
    input_files = sorted([os.path.join(second_arg, f) for f in os.listdir(second_arg) if f.endswith(".in")])
elif os.path.isfile(second_arg):
    input_files = [second_arg]
else:
    print("Invalid second argument")
    sys.exit(1)

for input_file in input_files:

    with open(input_file, "r") as f:
        input_data = f.read()

    try:
        result = subprocess.run(
            [program_path],
            input=input_data,
            text=True,
            capture_output=True,
            timeout=2
        )

        output = result.stdout

        lines = output.strip().splitlines()
        passed = True
        n, w, max_hash = len(lines), 0, 0
        hashes = dd(int)

        for line in lines:
            key, hash = map(int, line.split())
            hashes[hash] += 1
            max_hash = max(hash, max_hash)
            w = int(math.log2(max(key, 1))) + 1

        collisions = 0
        for sizes in hashes.values():
            collisions += sizes * (sizes - 1) // 2
        
        if collisions > 0 or max_hash >= 2 ** max(w / 2, int(math.log2(n) + 4)):
            passed = False
        
        if passed:
            print(f"Test {input_file}: Passed")
        else:
            print(f"Test {input_file}: Failed")
            print(f"No. of collisions: {collisions}")

    except subprocess.TimeoutExpired:
        print(f"Test {input_file}: Timeout")
    except Exception as e:
        print(f"Test {input_file}: Error: {e}")
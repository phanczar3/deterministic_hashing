import sys
import subprocess
import os
import math
from collections import defaultdict as dd

class FormatError(Exception):
    """Exception for wrong format code"""
    def __init__(self, message="Wrong format of a file"):
        super().__init__(message)


if len(sys.argv) not in [3, 4]:
    print("Invalid number of arguments")
    sys.exit(1)

program_path = sys.argv[1]
second_arg = sys.argv[2]
collision_bound = int(sys.argv[3]) if len(sys.argv) == 4 else 0

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
        input_lines = input_data.strip().splitlines()
        keys = set()
        msg = "Wrong input file format"    

        for i, line in enumerate(input_lines):
            lst = line.strip().split()
            if len(lst) > 1:
                raise FormatError(msg)
            try:
                x = int(lst[0])
                if i == 0:
                    n = x
                else:
                    if x in keys:
                        raise FormatError(msg)
                    keys.add(x)
            except ValueError:
                raise FormatError(msg)

        if n != len(keys):
            raise FormatError(msg)


        result = subprocess.run(
            [program_path],
            input=input_data,
            text=True,
            capture_output=True,
            timeout=2
        )

        output = result.stdout

        lines = output.strip().splitlines()

        msg = "Wrong output file format"
        keys_output = set()
        for line in lines:
            lst = line.strip().split()
            if len(lst) != 2:
                raise FormatError(msg)
            try:
                key, hash = map(int, line.split())

                if key in keys_output:
                    raise FormatError(msg)
                keys_output.add(key)
            except ValueError:
                raise FormatError(msg)

        if len(keys_output) != n:
            raise FormatError(msg)


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
        
        if collisions > collision_bound:
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
#!/usr/bin/python


import sys

if __name__ == "__main__":
    args = sys.argv[1:]
    if len(args) != 2:
        print("Usage: python3 CompareLists.py <list1> <list2>")
        sys.exit(1)

    list1_path = args[0]
    list2_path = args[1]

    with open(list1_path, 'r') as f:
        list1 = f.readlines()

    with open(list2_path, 'r') as f:
        list2 = f.readlines()

    # Remove whitespace and convert to sets
    set1 = set(line.strip() for line in list1)
    set2 = set(line.strip() for line in list2)

    # Find differences
    only_in_list1 = set1 - set2
    only_in_list2 = set2 - set1

    # Print results
    if only_in_list1:
        print("Items only in list 1:")
        for item in only_in_list1:
            print(item)

    if only_in_list2:
        print("Items only in list 2:")
        for item in only_in_list2:
            print(item)

    if not only_in_list1 and not only_in_list2:
        print("The lists are identical.")

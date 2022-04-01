"""
Predict Secondary Structure With Stochastic Context Free Grammar.
This program will find a secondary structure using a Context Free Grammar.

author: Reis Gadsden
version: 04/01/22
git: https://github.com/reismgadsden/predict_sec_struct_with_cfg

class: CS-5531-101 @ Appalachian State University
instructor: Dr. Mohammad Mohebbi
"""
import sys

memo = dict()


def main():
    sequence = "ACGU"
    print(predict_second_struct(sequence))
    print(memo)


def predict_second_struct(sequence):
    if sequence in memo:
        return memo[sequence]
    if len(sequence) > 1:
        k_args = []
        if is_pair(sequence[0], sequence[-1]):
            k_args.append("(" +  predict_second_struct(sequence[1:-1]) + ")")
        if len(sequence) > 2:
            for i in range(1, len(sequence) - 1):
                if is_pair(sequence[0], sequence[len(sequence) - i]):
                    k_args.append(predict_second_struct(sequence[0:len(sequence)-i]) + predict_second_struct(sequence[len(sequence)-i:len(sequence)]))
        k_args.append("." + predict_second_struct(sequence[1:]))
        best = ""
        for k in k_args:
            if count_base_pairs(k) > count_base_pairs(best):
                best = k
        memo[sequence] = best
        return memo[sequence]
    elif len(sequence) == 1:
        return "."
    elif len(sequence) == 0:
        return ""
    else:
        print("Something is really messed up")
        sys.exit(0)


def count_base_pairs(sequence) -> int:
    count = 0
    for c in sequence:
        if c == "(":
            count += 1
    return count


def is_pair(c1, c2) -> bool:
    look_up = {
        'A': ['U'],
        'C': ['G'],
        'G': ['C', 'U'],
        'T': ['A'],
        'U': ['A', 'G']
    }

    return c2 in look_up[c1]


if __name__ == "__main__":
    main()
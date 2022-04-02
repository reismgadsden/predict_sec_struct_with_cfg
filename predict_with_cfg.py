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
import get_stochastic_values as gsv

memo = dict()
rule_memo = dict()

memo_s = dict()
rule_memo_s = dict()

prob_dict = gsv.compute_stochastic_vals()


def main():
    sequence = "GGCGGUGAAAUGCC"
    standard_pair = predict_secondary_structure(sequence)
    stochastic_pair = predict_stochastic_secondary_structure(sequence)
    print("Sequence: " + sequence)
    print()
    print("Standard Pairing Method: " + standard_pair[0])
    print("Standard Pairing Probability: ", count_base_pairs_stochastically(sequence, standard_pair[0]))
    print("Standard Pairing Rules:")
    print(standard_pair[1])
    print()
    print("Stochastic CFG Method: " + stochastic_pair[0])
    print("Stochastic CFG Probability: ", count_base_pairs_stochastically(sequence, stochastic_pair[0]))
    print("Stochastic CFG Rules:")
    print(stochastic_pair[1])


def predict_secondary_structure(sequence):
    # base case 1: sequence has already been computed
    if sequence in memo and sequence in rule_memo:
        return [memo[sequence], rule_memo[sequence]]

    # base case 2: sequence has length of 0; represents rule S -> λ
    if len(sequence) == 0:
        return ["", "S -> λ"]

    # base case 3: sequence has length of 1; represents rule S -> A | C | G | U
    elif len(sequence) == 1:
        return [".", "S -> " + sequence.lower()]

    # recursive cases
    elif len(sequence) > 1:
        k_args = []

        # case 1: i and j are matching; represents rule S -> cSc (c ∈ {A, C, G, U})
        if is_pair(sequence[0], sequence[-1]):
            val = predict_secondary_structure(sequence[1:-1])
            val[0] = "(" + val[0] + ")"
            val[1] = "S -> " + sequence[0].lower() + "S" + sequence[-1].lower() + "\n" + val[1]
            k_args.append(val)

        # case 2: i binds to some k, use recursion to solve all possible k's;
        # represents rule S -> SS
        if len(sequence) > 2:
            # loop through all possible k values and append their resulting secondary structure
            # to an array for further validation
            for i in range(1, len(sequence) - 1):
                if is_pair(sequence[0], sequence[len(sequence) - i - 1]):
                    val1 = predict_secondary_structure(sequence[0:len(sequence) - i])
                    val2 = predict_secondary_structure(sequence[len(sequence) - i:len(sequence)])

                    val = [val1[0] + val2[0], "S -> S1S2 {\n" + "\tS1:\n\t\t" + val1[1].replace("\n",
                                                                                                "\n\t\t") + "\n" + " \tS2:\n\t\t" +
                           val2[1].replace("\n", "\n\t\t") + "\n}"]
                    k_args.append(val)

        # case 3: i binds to no other base; represents rule S -> cS (c ∈ {A, C, G, U})
        val = predict_secondary_structure(sequence[1:])
        val[0] = "." + val[0]
        val[1] = "S -> " + sequence[0].lower() + "S\n" + val[1]
        k_args.append(val)

        # loop through all possible cases and choose the one that has the most base pairs
        best = ["", ""]
        for k in k_args:
            if count_base_pairs(k[0]) >= count_base_pairs(best[0]):
                best = k
        memo[sequence] = best[0]
        rule_memo[sequence] = best[1]
        return [best[0], best[1]]

    # base case ?: I am a bad programmer
    else:
        print("Something is really messed up")
        sys.exit(0)


def predict_stochastic_secondary_structure(sequence):
    # base case 1: sequence has already been computed
    if sequence in memo_s and sequence in rule_memo_s:
        return [memo_s[sequence], rule_memo_s[sequence]]

    # base case 2: sequence has length of 0; represents rule S -> λ
    if len(sequence) == 0:
        return ["", "S -> λ"]

    # base case 3: sequence has length of 1; represents rule S -> A | C | G | U
    elif len(sequence) == 1:
        return [".", "S -> " + sequence.lower()]

    # recursive cases
    elif len(sequence) > 1:
        k_args = []

        # case 1: i and j are matching; represents rule S -> cSc (c ∈ {A, C, G, U})
        if is_pair(sequence[0], sequence[-1]):
            val = predict_stochastic_secondary_structure(sequence[1:-1])
            val[0] = "(" + val[0] + ")"
            val[1] = "S -> " + sequence[0].lower() + "S" + sequence[-1].lower() + "\n" + val[1]
            k_args.append(val)

        # case 2: i binds to some k, use recursion to solve all possible k's;
        # represents rule S -> SS
        if len(sequence) > 2:
            # loop through all possible k values and append their resulting secondary structure
            # to an array for further validation
            for i in range(1, len(sequence) - 1):
                if is_pair(sequence[0], sequence[len(sequence) - i - 1]):
                    val1 = predict_stochastic_secondary_structure(sequence[0:len(sequence)-i])
                    val2 = predict_stochastic_secondary_structure(sequence[len(sequence)-i:len(sequence)])

                    val = [val1[0] + val2[0], "S -> S1S2 {\n" + "\tS1:\n\t\t" + val1[1].replace("\n", "\n\t\t") + "\n" + " \tS2:\n\t\t" + val2[1].replace("\n", "\n\t\t") + "\n}"]
                    k_args.append(val)

        # case 3: i binds to no other base; represents rule S -> cS (c ∈ {A, C, G, U})
        val = predict_stochastic_secondary_structure(sequence[1:])
        val[0] = "." + val[0]
        val[1] = "S -> " + sequence[0].lower() + "S\n" + val[1]
        k_args.append(val)

        # loop through all possible cases and choose the one that has the most base pairs
        best = ["", ""]
        for k in k_args:
            if count_base_pairs_stochastically(sequence, k[0]) >= count_base_pairs_stochastically(sequence, best[0]):
                best = k
        memo_s[sequence] = best[0]
        rule_memo_s[sequence] = best[1]
        return [best[0], best[1]]

    # base case ?: I am a bad programmer
    else:
        print("Something is really messed up")
        sys.exit(0)


# this method returns a int value representing the total amount of base pairs
# in a given secondary structure
def count_base_pairs(sec_struct) -> int:
    count = 0
    for c in sec_struct:
        if c == "(":  # only need to count one bracket
            count += 1
    return count


def count_base_pairs_stochastically(sequence, secondary_struct):
    if secondary_struct == "":
        return 0
    counter = 0
    stack = []
    prob = 0
    for c in secondary_struct:
        if c == "(":
            stack.append(["(", counter])
        elif c == ")":
            val = stack.pop()
            prob += prob_dict[determine_pair(sequence[val[1]], sequence[counter])]
        counter += 1

    return prob


# this method returns a boolean value representing if two
# bases can bind
def is_pair(c1, c2) -> bool:
    look_up = {
        'A': ['U'],
        'C': ['G'],
        'G': ['C', 'U'],
        'T': ['A'],
        'U': ['A', 'G']
    }

    return c2 in look_up[c1]


def determine_pair(s1, s2):
    if ord(s1) > ord(s2):
        return s2 + s1
    return s1 + s2


# if this file is not being imported, run main method
if __name__ == "__main__":
    main()

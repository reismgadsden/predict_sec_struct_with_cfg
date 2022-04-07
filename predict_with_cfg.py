"""
Predict Secondary Structure With Stochastic Context Free Grammar.
This program will find a secondary structure using a Context Free Grammar.

author: Reis Gadsden
version: 04/01/22
git: https://github.com/reismgadsden/predict_sec_struct_with_cfg

class: CS-5531-101 @ Appalachian State University
instructor: Dr. Mohammad Mohebbi
"""
# necessary import
import sys
import get_stochastic_values as gsv

# memorization dictionaries for standard method
memo = dict()
rule_memo = dict()

# memorization dictionaries for stochastic method
memo_s = dict()
rule_memo_s = dict()

# probability values from get_stochastic values.py
prob_dict = gsv.compute_stochastic_vals()

"""
Main method of the program, builds the output.
Takes no arguments only makes calls to and print
out the output of our two methods.
"""
def main():
    sequence = "UGAAAGGGUUAUUACUACUAAUACAAGACCAGUUCGUGUCGCACAUCCACGCUUUUGACCAGUCAGAUAGCGGUGGAGCCCAACGUAACGAGCAGUGGGU"
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


"""
Prediction of secondary structure via method of
maximizing base pairs.

:param sequence - the sequence that is to be recursively processed.

:returns a list of length two of our outputs, where index 0 is the
         secondary structure, and index 1 is the parse tree.
"""
def predict_secondary_structure(sequence) -> list:
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

        # case 1: i and j are matching; represents rule S -> xSy (x,y ∈ {A, C, G, U})
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

                    val = [val1[0] + val2[0], "S -> S1S2 {\n" + "\tS1:\n\t\t" + val1[1].replace("\n", "\n\t\t") \
                           + "\n" + " \tS2:\n\t\t" + val2[1].replace("\n", "\n\t\t") + "\n}"]
                    k_args.append(val)

        # case 3: i binds to no other base; represents rule S -> xS (x ∈ {A, C, G, U})
        val = predict_secondary_structure(sequence[1:])
        val[0] = "." + val[0]
        val[1] = "S -> " + sequence[0].lower() + "S\n" + val[1]
        k_args.append(val)

        # loop through all possible cases and choose the one that has the most base pairs
        # append our results to memorization tables for future use
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


"""
Prediction of secondary structure via method of
maximizing stochastic probability.

:param sequence - the sequence that is to be recursively processed.

:returns a list of length two of our outputs, where index 0 is the
         secondary structure, and index 1 is the parse tree.
"""
def predict_stochastic_secondary_structure(sequence) -> list:
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
                    val1 = predict_stochastic_secondary_structure(sequence[0:len(sequence) - i])
                    val2 = predict_stochastic_secondary_structure(sequence[len(sequence) - i:len(sequence)])

                    val = [val1[0] + val2[0], "S -> S1S2 {\n" + "\tS1:\n\t\t" + val1[1].replace("\n", "\n\t\t") \
                           + "\n" + " \tS2:\n\t\t" + val2[1].replace("\n", "\n\t\t") + "\n}"]
                    k_args.append(val)

        # case 3: i binds to no other base; represents rule S -> cS (c ∈ {A, C, G, U})
        val = predict_stochastic_secondary_structure(sequence[1:])
        val[0] = "." + val[0]
        val[1] = "S -> " + sequence[0].lower() + "S\n" + val[1]
        k_args.append(val)

        # loop through all possible cases and choose the one that has the highest probability
        # append our results to memorization tables for future use
        best = ["", ""]
        for k in k_args:
            if count_base_pairs_stochastically(sequence, k[0]) \
                    >= count_base_pairs_stochastically(sequence, best[0]):
                best = k
        memo_s[sequence] = best[0]
        rule_memo_s[sequence] = best[1]
        return [best[0], best[1]]

    # base case ?: I am a bad programmer
    else:
        print("Something is really messed up")
        sys.exit(0)

"""
This method returns an int value representing the total amount of base pairs
in a given secondary structure.

:param sec_struct - the secondary structure to count the base pairs of

:returns an integer value representing number of base pairs
"""
def count_base_pairs(sec_struct) -> int:
    # value to keep count of base pairs
    count = 0

    # loop through each character in the secondary structure
    # and iterate the count everytime an open parenthesis appears.
    # we only need to check for opening, because it is impossible
    # for there to be an unmatched open or closed bracket.
    for c in sec_struct:
        if c == "(":
            count += 1

    # return the count of base pairs
    return count


"""
This method returns a float value representing the stochastic probability
of a given secondary structure for a given sequence.

:param sequence - sequence that the secondary structure represents
:param secondary_struct - corresponding secondary struct for the sequence

:returns a float value representing the stochastic probability
"""
def count_base_pairs_stochastically(sequence, secondary_struct) -> float:
    # in the comparison part of the stochastic method we initialize
    # a variable best, best starts off with an empty string for the
    # secondary structure. This statement is to avoid a fatal error.
    if secondary_struct == "":
        return 0

    # counter variable to keep track of the current index
    counter = 0

    # a "stack" to keep track of matching base pairs
    stack = []

    # a float value to keep track of total probability
    prob = 0

    # loop through each character in the secondary structure and push the indexes
    # of open parenthesis onto the stack whenever an open parenthesis is read. When
    # a closed parenthesis is read, we pop the stack, and get the corresponding
    # probabilities of the pair found at the two indexes, and iterate the
    # probability with the returned value.
    for c in secondary_struct:
        if c == "(":
            stack.append(counter)
        elif c == ")":
            val = stack.pop()
            prob += prob_dict[determine_pair(sequence[val], sequence[counter])]
        counter += 1

    # return dat bih
    return prob


"""
This method will tell you if two bases are a pair.

# order of params is not important
:param c1 - the first base
:param c2 - the second base 

:returns boolean representation of whether two bases can make a pair
"""
def is_pair(c1, c2) -> bool:
    # a look-up table for quick access
    look_up = {
        'A': ['U'],
        'C': ['G'],
        'G': ['C', 'U'],
        'T': ['A'],
        'U': ['A', 'G']
    }

    # return whether the two bases can pair
    return c2 in look_up[c1]


"""
This method will take in two bases and arrange them alphabetically.
This is to improve QOL when accessing the stochastic probability dictionary.

# order of params is not important, thats the point of this method, lol.
:param c1 - the first base
:param c2 - the second base

:returns the two bases in alphabetical order
"""
def determine_pair(c1, c2):
    if ord(c1) > ord(c2):
        return c2 + c1
    return c1 + c2


"""
If this file is not being imported, run main method.
"""
if __name__ == "__main__":
    main()

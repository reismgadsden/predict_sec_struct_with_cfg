import os


stochastic_vals = {
    "AU": 0,
    "CG": 0,
    "GU": 0
}


total_sequence_dict = {
    "AU": 0,
    "CG" : 0,
    "GU" : 0,
    "Total": 0
}


def main():
    print(compute_stochastic_vals())


def read_all_files() -> list:
    container = list()
    dir_path = ".\\real_sec_structures"
    dir_content = os.listdir(dir_path)

    for x in dir_content:
        f = open(dir_path + "\\" + x, 'r')
        output = f.readlines()
        container.append([output[1], output[2], output[0][1:5]])
        f.close()

    return container


def find_matching_offset(second_struct):
    stack = []
    i = 0
    for x in second_struct:
        if x == "(":
            stack.append(i)
        elif x == ")":
            stack.pop()

        if not stack:
            break
        i += 1
    return i


def compute_stochastic_vals() -> dict:
    all_seq = read_all_files()
    for a in all_seq:
        sequence_dict = {
            "AU": 0,
            "CG": 0,
            "GU": 0,
            "Total": 0
        }
        if not get_indexes(a[0], a[1], sequence_dict):
            print("Secondary Structure is not Properly Formatted.")

    for x in total_sequence_dict.keys():
        if x != "Total":
            stochastic_vals[x] = total_sequence_dict[x]/ total_sequence_dict["Total"]

    stochastic_vals["A"] = 0
    stochastic_vals["C"] = 0
    stochastic_vals["G"] = 0
    stochastic_vals["U"] = 0
    return stochastic_vals


def get_indexes(sequence, second_struct, single_dict=None) -> bool:
    stack = []
    i = 0
    for x in second_struct:
        if x == "(":
            stack.append(i)
        elif x == ")":
            if single_dict is None:
                count_sequences(stack.pop(), i, sequence)
            else:
                count_sequences(stack.pop(), i, sequence, single_dict)
        i += 1
    return len(stack) == 0


def count_sequences(open_index, close_index, sequence, single_dict=None) -> None:
    open_seq_char = sequence[open_index].upper()
    close_seq_char = sequence[close_index].upper()

    if ord(open_seq_char) > ord(close_seq_char):
        temp = close_seq_char
        close_seq_char = open_seq_char
        open_seq_char = temp
    if close_seq_char == 'T':
        close_seq_char = 'U'
    total_sequence_dict[open_seq_char + close_seq_char] += 1
    total_sequence_dict["Total"] += 1
    if single_dict is not None:
        single_dict[open_seq_char + close_seq_char] += 1
        single_dict["Total"] += 1

if __name__ == "__main__":
    main()
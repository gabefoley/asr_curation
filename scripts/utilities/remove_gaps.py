import sys
import sequence

if __name__ == "__main__":

    print("\nRemoving gaps")

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    print(input_file)

    seqs = sequence.readFastaFile(input_file)

    for seq in seqs:
        seq.sequence = [x for x in seq.sequence if x != "-"]

    sequence.writeFastaFile(output_file, seqs)

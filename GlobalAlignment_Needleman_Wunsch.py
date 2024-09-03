
def extract_sequence_fasta(filename: str) -> str:
    """
    This function extracts the sequence from the provided FASTA file.
    :param filename (str): name of the FASTA file

    :return str: extracted sequence
    """
    with open(filename, 'r') as f:
        first_line = f.readline() # Check the first line to understand if it's a FASTA file
        if not first_line.startswith('>'):
            print('\x1b[1;31;49m' + "[ERROR] <{}> is not a FASTA file.".format(filename) + '\x1b[0m')
            sys.exit()
        else:
            sequence = ''
            for line in f:
                sequence += line.strip('\n')
            return sequence


def check_sequence (sequence: str, sequence_number: int, type: str, alphabet: set, q: bool) -> None:
    """
    This function verifies that all characters in the input string are within
    the set of allowed characters provided as input.
    The function stops the execution of the script when an unexpected character is found.

    :param sequence (str): sequence to be checked
    :param sequence_number (int): number to indicate if the input sequence is the first or the second to align (possible values: 1 or 2)
    :param type (str): type of the sequence (possible values: 'd','r' or 'p')
    :param alphabet (set): set containing all the allowed characters
    :param q (bool): boolean variable to set quiet mode (True = quiet mode activated)

    :return: None
    """
    # To show the progress bar (if quiet mode is disabled).
    for index, letter in tqdm(enumerate(sequence), total = len(sequence),
                              desc = '\x1b[1;32;49m' + '[CHECK {}] Checking sequence correctness...'.format(sequence_number)
                                     + '\x1b[0m' , ncols = 100, disable = q):
        if letter not in alphabet:
            print('\x1b[1;31;49m' + "\n[ERROR] {} sequence is invalid.".format('First' if sequence_number==1 else 'Second') + '\x1b[0m')
            print('\x1b[1;31;49m' + "[ERROR] '{}' at position {} (zero-based index) is an invalid character.".format(letter, index) + '\x1b[0m')
            print('\x1b[1;31;49m' + "[ERROR] Please enter a valid {} sequence.".format("DNA" if type=='d' else "RNA" if type=='r' else "protein") + '\x1b[0m')
            sys.exit()


def compute_diagonal_indices(rows: 'numpy.array', cols: 'numpy.array', nrows: int) -> tuple['numpy.array','numpy.array']:
    """
        This function computes column and row indices of the elements of the diagonal to compute,
        given in input the column and row indices of the elements of the previous diagonal.

        :param rows (numpy.array): row indices of the elements of the previous diagonal
        :param cols (numpy.array): column indices of the elements of the previous diagonal
        :param nrows (int): number of rows of the matrix (to know its boundaries)

        :return: tuple['numpy.array','numpy.array']: tuple containing as first element the new row indices and as second
                                                     element the new column indices
        """
    if np.min(cols) > 0:
        cols -= 1
        if np.max(rows) < nrows - 1:
            cols = np.append(cols, np.max(cols) + 1)
            rows = np.append(rows, np.max(rows) + 1)
    else:
        cols = cols[:-1:1]
        rows = rows[1::1]
    return rows, cols

def arguments_generator(seq1: str, seq2: str, curr: 'numpy.array', abv: 'numpy.array', abv_abv: 'numpy.array',
                        GAP: int, MISMATCH: int, MATCH: int, row_ind: 'numpy.array',
                        col_ind: 'numpy.array', score_matrix: 'blosum._blosum.BLOSUM' or None) -> list[list, ...]:
    """
    This function generates all the arguments required for each child process, which is responsible of the computation
    of one cell of the diagonal to be currently computed.

    :param seq1 (str): first sequence to be aligned
    :param seq2 (str): second sequence to be aligned
    :param curr (numpy.array): current diagonal
    :param abv (numpy.array): diagonal above the current one
    :param abv_abv (numpy.array): diagonal above the 'abv' diagonal
    :param GAP (int): GAP penalty
    :param MISMATCH (int): MISMATCH penalty
    :param MATCH (int): MATCH score
    :param row_ind (numpy.array): array containing row indices of the elements of 'curr' diagonal
    :param col_ind (numpy.array): array containing column indices of the elements of 'curr' diagonal
    :param score_matrix (None | blosum._blosum.BLOSUM): blosum score matrix

    :return: list[list,...]: list containing lists, each of which corresponds to the arguments necessary to compute one cell of the current diagonal
    """
    all_arguments = []
    if len(curr) < len(abv):  # Case 1 (the current diagonal is shorter than the above diagonal, so we are in the lower triangle of the matrix).
        for i in range(len(curr)):
            cell_above = int(abv[i])
            cell_lateral = int(abv[i+1])

            if len(curr) == len(abv_abv) or len(abv) == len(abv_abv):
                cell_diagonal = int(abv_abv[i])
            else:
                cell_diagonal = int(abv_abv[i + 1])
            characters_aligned = seq1[col_ind[i]] + seq2[row_ind[i]]

            if score_matrix is not None and not '-' in characters_aligned:
                blosum_score = score_matrix[characters_aligned[0]][characters_aligned[1]]
            else:
                blosum_score = None

            all_arguments.append([cell_above, cell_lateral, cell_diagonal, GAP, MISMATCH,
                                     MATCH, characters_aligned, blosum_score])
        return all_arguments

    elif len(curr) > len(abv):  # Case 2 (the current diagonal is longer than the above diagonal, so we are in the upper triangle of the matrix).
        for i in range(len(curr)):

            if i == 0: # The first cell won't have cells above and diagonal to it.
                cell_above = None
                cell_lateral = int(abv[max(0, i-1)])
                cell_diagonal = None

            elif i == len(curr) - 1: # The last cell won't have cells lateral and diagonal to it.
                cell_above = int(abv[max(0, i - 1)])
                cell_lateral = None
                cell_diagonal = None

            else:
                cell_above = int(abv[i-1])
                cell_lateral = int(abv[i])
                cell_diagonal = int(abv_abv[i-1])

            characters_aligned = seq1[col_ind[i]] + seq2[row_ind[i]]
            if score_matrix is not None and cell_diagonal is not None and not '-' in characters_aligned:
                blosum_score = score_matrix[characters_aligned[0]][characters_aligned[1]]
            else:
                blosum_score = None

            all_arguments.append([cell_above, cell_lateral, cell_diagonal, GAP, MISMATCH,
                                      MATCH, characters_aligned, blosum_score])
        return all_arguments

    else: # Case 3 (the current diagonal is long as the above diagonal).
        for i in range(len(curr)):
            if i == 0:
                cell_above = None
                cell_diagonal = None

            else:
                cell_above = int(abv[i-1])
                cell_diagonal = int(abv_abv[i-1])

            cell_lateral = int(abv[i])
            characters_aligned = seq1[col_ind[i]] + seq2[row_ind[i]]
            if score_matrix is not None and cell_diagonal is not None and not '-' in characters_aligned:
                blosum_score = score_matrix[characters_aligned[0]][characters_aligned[1]]
            else:
                blosum_score = None

            all_arguments.append([cell_above, cell_lateral, cell_diagonal, GAP, MISMATCH,
                                  MATCH, characters_aligned, blosum_score])
        return all_arguments


def compute_score_and_move(cell_above: int | None, cell_lateral: int | None, cell_diagonal: int | None,
                    GAP: int, MISMATCH: int, MATCH: int, characters_aligned: str,
                           blosum_score: int | None) -> tuple[int, int]:
    """
    This function computes both the partial alignment score and the move done for a single cell, given in input
    the scores needed to evaluate the partial alignment and the partial alignment values of all adjacent cells.

    :param cell_above (int | None): score of the cell above the cell to be computed by the call of this function
                                    (can be None when the current cell belongs to the first row)
    :param cell_lateral (int | None): score of the cell right to the cell to be computed by the call of this function
                                      (can be None when the current cell belongs to the last column)
    :param cell_diagonal (int | None): score of the cell diagonally up and to the right from the cell to be computed by the call of this function
                                       (can be None when the current cell belongs to either the first row or the last column)
    :param GAP (int): GAP penalty value
    :param MISMATCH (int): MISMATCH value
    :param MATCH (int): MATCH value
    :param characters_aligned (str): the two characters alignable by this cell
    :param blosum_score (int | None): the BLOSUM score of the alignment of the two characters in 'characters_aligned'
                                      (can be None when one of the two characters in 'characters_aligned' is the gap symbol '-' or if
                                      the user didn't want to use it or if the two sequences to be aligned are not protein sequences)

    :return tuple([int], [int]): tuple with the value and the move associated to this cell, i.e. tuple(value, move)
    """
    if cell_lateral is None:
        value = cell_above + GAP
        move = 1                         # Move = 1: vertical move.
        return value, move
    elif cell_above is None:
        value = cell_lateral + GAP
        move = 0                         # Move = 0: horizontal move.
        return value, move
    vertical = cell_above + GAP
    horizontal = cell_lateral + GAP

    if blosum_score is None:
        diagonal = cell_diagonal + (MATCH if characters_aligned[0]==characters_aligned[1] else MISMATCH)
    else:
        diagonal = cell_diagonal + blosum_score

    value = max(vertical, horizontal, diagonal)
    if value==vertical and value==horizontal and value==diagonal:
        move = 6             # Move = 6: vertical move + horizontal move + diagonal move.
        return value, move
    if value==vertical and value==horizontal:
        move = 4             # Move = 4: vertical move + horizontal move.
        return value, move
    if value==horizontal and value==diagonal:
        move = 3             # Move = 3: horizontal move + diagonal move.
        return value, move
    if value==vertical and value==diagonal:
        move = 5             # Move = 5: vertical move + diagonal move.
        return value, move
    if value==diagonal:
        move = 2             # Move = 2: diagonal move.
        return value, move
    if value==vertical:
        move = 1
        return value, move

    move = 0
    return value, move

def recursive_traceback(df: 'pd.DataFrame', row_idx: int, col_idx: int,
                        heuristic_token: bool, seq = "") -> str:
    """
    This function performs the traceback starting from the lower-rightmost cell
    to the upper-leftmost cell.
    The traceback is recursively implemented to identify all the alignments
    giving the same score, which is also the best possible score.
    The function returns a string containing all the informations needed to print
    the possibly several alignments: the informations for each alignment are divided by ':'.
    For example, if there are two equivalent alignments, the output string will be 'AATTCA:A-TT-C',
    which corresponds to:
        Alignment 1
     SEQ1: ATC
     SEQ2: ATA
        Alignment 2
     SEQ1: AT-
     SEQ2: -TC

    :param df (pandas.core.frame.DataFrame): Pandas dataframe containing the 'moves' to perform the traceback
    :param row_idx (int): Current row index
    :param col_idx (int): Current column index
    :param seq (str): Dynamically growing string

    :return str: String shaped like 'alignment1:alignment2:...:alignmentN' in order to
                 have all the possible alignments with the same (best) score
    """
    move = int(df.iloc[row_idx, col_idx])   # Current move to evaluate.

    if heuristic_token == True:   # In this case, whenever a 'complex' move is encountered (move >= 3), a simpler move is chosen at random so to not introduce bias in the choice.
        if move == 6:
            move = random.choice([3,4,5])  # Moves 3,4 and 5 are less complex than move 6 because they comprise two movements instead of three.
        elif move == 5:
            move = random.choice([1,2])    # Move 1 and 2 are less complex than move 5 because they comprise one movement instead of two.
        elif move == 4:
            move = random.choice([0,1])    # Same holds for move 0.
        elif move == 3:
            move = random.choice([0,2])

    lett1 = df.columns[col_idx]  # Letter of sequence 1 to be considered.
    lett2 = df.index[row_idx]  # Letter of sequence 2 to be considered.

    if move == -1:  # END OF TRACEBACK (-1 is the value in the upper-leftmost cell).
        return seq
    elif move == 0: # horizontal move
        bit_seq = lett1 + "-"
        seq = bit_seq + seq
        return recursive_traceback(df, row_idx, col_idx - 1, heuristic_token, seq = seq)
    elif move == 1: # vertical move
        bit_seq = "-" + lett2
        seq = bit_seq + seq
        return recursive_traceback(df, row_idx - 1, col_idx, heuristic_token, seq = seq)
    elif move == 2: # diagonal move
        bit_seq = lett1 + lett2
        seq = bit_seq + seq
        return recursive_traceback(df, row_idx - 1, col_idx -1, heuristic_token, seq = seq)
    elif move == 3: # diagonal move + horizontal move, so this is a fork that gives two alignments with the same (best) score.
        bit_seq_diag = lett1 + lett2
        seq_diag = bit_seq_diag + seq
        bit_seq_horiz = lett1 + '-'
        seq_horiz = bit_seq_horiz + seq
        return (recursive_traceback(df, row_idx - 1, col_idx - 1, heuristic_token, seq = seq_diag)
                + ":" + recursive_traceback(df, row_idx, col_idx - 1, heuristic_token, seq = seq_horiz))
    elif move == 4: # horizontal move + vertical move, so this is a fork that gives two alignments with the same (best) score.
        bit_seq_horiz = lett1 + '-'
        seq_horiz = bit_seq_horiz + seq
        bit_seq_vert = '-' + lett2
        seq_vert = bit_seq_vert + seq
        return (recursive_traceback(df, row_idx, col_idx - 1, heuristic_token, seq = seq_horiz)
                + ":" + recursive_traceback(df, row_idx - 1, col_idx, heuristic_token, seq = seq_vert))
    elif move == 5: # diagonal move + vertical move, so this is a fork that gives two alignments with the same (best) score.
        bit_seq_diag = lett1 + lett2
        seq_diag = bit_seq_diag + seq
        bit_seq_vert = '-' + lett2
        seq_vert = bit_seq_vert + seq
        return (recursive_traceback(df, row_idx - 1, col_idx - 1, heuristic_token, seq = seq_diag)
                + ":" + recursive_traceback(df, row_idx - 1, col_idx, heuristic_token, seq = seq_vert))
    elif move == 6: # diagonal move + horizontal move + vertical move, so this is a trivium that gives three alignments with the same (best) score.
        bit_seq_diag = lett1 + lett2
        seq_diag = bit_seq_diag + seq
        bit_seq_horiz = lett1 + '-'
        seq_horiz = bit_seq_horiz + seq
        bit_seq_vert = '-' + lett2
        seq_vert = bit_seq_vert + seq
        return ((recursive_traceback(df, row_idx - 1, col_idx - 1, heuristic_token, seq = seq_diag)
                + ":" + recursive_traceback(df, row_idx - 1, col_idx, heuristic_token, seq = seq_vert))
                + ":" + recursive_traceback(df, row_idx, col_idx - 1, heuristic_token, seq = seq_horiz))


if __name__ == '__main__':
    # Libraries are imported in the __main__ otherwise all child-processes will need to import
    # them every time (in Windows), increasing tremendously the execution time.
    import pandas as pd, numpy as np, sys, multiprocessing, argparse, random, blosum as bl
    from tqdm import tqdm

    parser = argparse.ArgumentParser(prog = 'GlobalAlignment_Needleman_Wunsch.py',
                                     description = 'This program computes in parallel the optimal global alignment'
                                                   ' of two sequences (DNA, RNA or protein sequences).',
                                     epilog = 'Marco Cominelli, 2024')

    group = parser.add_mutually_exclusive_group(required = True)
    group.add_argument('-s', '--sequence', type = str, nargs = 2, required = False,
                       help = 'the two sequences to be aligned')
    group.add_argument('-f', '--fasta', type = str, nargs = 2, required = False,
                       help = 'the two FASTA files containing the sequences to be aligned')
    parser.add_argument('-t', '--type', type = str, choices = ['d', 'r', 'p'], default = 'd',
                       help = "type of sequences to align: 'd' for DNA sequences, 'r' for RNA sequences and "
                              "'p' for protein sequences. Default is 'd' ")
    parser.add_argument('-g', '--gap', type = int, default = '-4',
                        help = 'negative GAP penalty. Default is -4')
    parser.add_argument('-mm', '--mismatch', type = int, default = '-5',
                        help = 'negative MISMATCH penalty. Default is -5')
    parser.add_argument('-m', '--match', type = int, default = '5',
                        help = 'positive MATCH score. Default is 5')
    parser.add_argument('-b', '--blosum', type = int, default = None,
                        help = 'use the specified BLOSUM matrix for MATCHES/MISMATCHES '
                               '(compatible only with protein sequences)')
    parser.add_argument('-o','--output', type = str, default = None,
                        help = "save the alignment(s) in the specified output file")
    parser.add_argument('-q', '--quiet', action = 'store_true',
                        help = "don't display the output")
    parser.add_argument('-c', '--cores', type = int, default = 3,
                        help = 'number of cores to use. Default is 3 (if available)')
    parser.add_argument('-v', '--verbose', type = int, default = 1, choices = [1, 2],
                        help = 'increase verbosity. Default is 1')
    parser.add_argument('-a', '--approximation', action = 'store_true',
                        help = "don't show all possible alignments but a reduced number of them")

    args = parser.parse_args()  # It contains all the arguments and their respective values.

    if args.fasta is None:    # Case in which the user entered the sequences from command line.
        sequence1, sequence2 = args.sequence[:]
        sequence1, sequence2 = sequence1.upper(), sequence2.upper()
    else:     # Case in which the user provided two FASTA files containing the sequences to be aligned.
        fasta1, fasta2 = args.fasta
        sequence1 = extract_sequence_fasta(fasta1).upper()
        sequence2 = extract_sequence_fasta(fasta2).upper()

    sequence_type = args.type
    GAP = args.gap
    MISMATCH = args.mismatch
    MATCH = args.match

    if GAP > 0:
        print('\x1b[1;31;49m' + "[ERROR] Gap penalty must be negative." + '\x1b[0m')
        sys.exit()

    if MISMATCH > 0:
        print('\x1b[1;31;49m' + "[ERROR] Mismatch penalty must be negative." + '\x1b[0m')
        sys.exit()

    if MATCH < 0:
        print('\x1b[1;31;49m' + "[ERROR] Match score must be positive." + '\x1b[0m')
        sys.exit()

    blosum_number = args.blosum
    output_filename = args.output
    quiet = args.quiet
    verbose = args.verbose
    heuristic_token = args.approximation
    if args.cores <= multiprocessing.cpu_count():
        num_processes = args.cores
    else:
        print('\x1b[1;93;49m' + "[WARNING] The number of cores asked is higher than the number of available cores."
              " The number of cores will be set to the maximum possible: {}.".format(multiprocessing.cpu_count())
              + '\x1b[0m')
        num_processes = multiprocessing.cpu_count()

    # Check: when the user want to use a BLOSUM matrix as scoring matrix,
    # the sequences to be aligned must be protein sequences.
    if blosum_number is not None and sequence_type != 'p':
        print('\x1b[1;31;49m' + "[ERROR] The usage of any BLOSUM matrix as scoring matrix is "
                                "available only for protein sequences."
                                " Use '-t p' or '--type p' to enable alignment of protein sequences."  + '\x1b[0m')
        sys.exit()

    blosum_matrix = None
    if blosum_number is not None:
        try:
            blosum_matrix = bl.BLOSUM(blosum_number)
        except BaseException:
            print('\x1b[1;31;49m' + "[ERROR] Unknown BLOSUM number", blosum_number,"\b. Choose number ϵ "
                                                                                   "{45,50,62,80,90}." + '\x1b[0m')
            sys.exit()

    DNA_ALPHABET = {'A', 'T', 'C', 'G'}  # Set of allowed characters in DNA sequences.

    RNA_ALPHABET = {'A', 'U', 'C', 'G'}  # Set of allowed characters in RNA sequences.

    PROTEIN_ALPHABET = {'A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I',  # Set of allowed characters in protein sequences.
                        'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'}

    if sequence_type == 'd':
        alphabet = DNA_ALPHABET
    elif sequence_type == 'r':
        alphabet = RNA_ALPHABET
    else:
        alphabet = PROTEIN_ALPHABET

    # Check for the correctness of the first sequence.
    check_sequence(sequence1, 1, sequence_type, alphabet, quiet)

    # Check for the correctness of the second sequence.
    check_sequence(sequence2, 2, sequence_type, alphabet, quiet)

    # The longest sequence is placed horizontally in the first line of the matrix, while the shortest sequence is placed
    # vertically in the first column. Specifically, in this implementation, the sequence placed in the first line
    # will be reversed (see the documentation for the explanation).

    original_sequence1, original_sequence2 = sequence1, sequence2     # Used at the end when printing the alignment.

    if len(sequence2) <= len(sequence1):
        sequence1 = sequence1[::-1] + "-"
        sequence2 = "-" + sequence2
    else:
        seq2_copy = sequence2
        sequence2 = "-" + sequence1
        sequence1 = seq2_copy[::-1] + "-"

    # This 3D array will contain the partial alignment scores in the first layer and
    # the values corresponding to the moves done in the second layer.
    array_scores_and_moves = np.zeros(( 2, len(sequence2), len(sequence1) ))

    # Offset corresponding to the first diagonal to compute.
    offset = len(sequence1) - 2

    # Offset corresponding to the last diagonal to compute.
    last_offset = -len(sequence2) + 1

    # Those variables are used later to compute the indices of the diagonal elements, at each iteration.
    number_rows = len(sequence2)
    row_indices = np.array(0)
    col_indices = np.array(len(sequence1) - 1)

    total_iterations = len(sequence1) - 2 + len(sequence2) - 1  # To show the progress bar (if quiet mode is disabled)
    with tqdm(total = total_iterations, desc =  '\x1b[1;32;49m' + '[SCORING] Computing scores and respective moves...' + '\x1b[0m', disable = quiet) as pbar:
        i = 0
        # 'While' loop for looping on diagonals.
        while offset >= last_offset:

            row_indices, col_indices = compute_diagonal_indices(row_indices, col_indices, number_rows)

            current = np.diag(array_scores_and_moves[0], k = offset)           # Diagonal to be computed.
            above = np.diag(array_scores_and_moves[0], k = offset + 1)         # Diagonal above the 'current' diagonal.
            above_above = np.diag(array_scores_and_moves[0], k = offset + 2)   # Diagonal above the 'above' diagonal.

            effective_num_processes = min(len(current), num_processes)  # If the cells to be computed in parallel are less than the cores asked in input,
                                                                        # there is no meaning in starting more processes than just the ones needed.
            p = multiprocessing.Pool(effective_num_processes)

            # 'arguments' is a list of list, in which each list contains the input parameters for one child process,
            # responsible of the computation of one cell of the current diagonal.
            arguments = arguments_generator(sequence1, sequence2, current, above, above_above, GAP,
                                            MISMATCH, MATCH, row_indices, col_indices, blosum_matrix)

            output = np.array(p.starmap(compute_score_and_move, arguments))  # 'output' is a numpy array --> [[value1, move1], [value2, move2],...,[valueN, moveN]].

            array_scores_and_moves[0][row_indices, col_indices] = output[:,0]
            array_scores_and_moves[1][row_indices, col_indices] = output[:,1]

            offset -= 1    # To loop on the next diagonal at the next iteration.

            # To update the progress bar.
            pbar.update(1)
            i += 1

    # Pandas dataframe with the scores, but rotated in the usual way (i.e. with the sequence above written in the original orientation).
    df_scores = pd.DataFrame(np.fliplr(array_scores_and_moves[0]),index=[i for i in sequence2], columns = [i for i in sequence1[::-1]])
    # Pandas dataframe with the moves, but rotated in the usual way (i.e. with the sequence above written in the original orientation).
    df_moves = pd.DataFrame(np.fliplr(array_scores_and_moves[1]), index=[i for i in sequence2], columns=[i for i in sequence1[::-1]])

    if quiet is False:
        if verbose == 2:
            print("\nMatrix of scores:")
            print(df_scores, "\n")

            print("Matrix of moves for the traceback:")
            print("The numbers have the following meanings:")
            print("0) Horizontal move")
            print("1) Vertical move")
            print("2) Diagonal move")
            print("3) Horizontal move + diagonal move")
            print("4) Vertical move + horizontal move")
            print("5) Vertical move + diagonal move")
            print("6) Vertical move + diagonal move + horizontal move")

            print(df_moves, "\n")

        print('\x1b[1;32;49m' + "[TRACEBACK] Starting traceback..." + '\x1b[0m' )

    df_moves.loc["-","-"] = -1    # The first cell is set to '-1' to stop the recursive function when it reaches this cell.

    # If two sequences are too long there is the risk of RecursionError or a waiting time too long,
    # so in this case an heuristic is used. See the documentation for the full explanation.

    # Even if the sequences are not too long, the risk of RecursionError still exists when they are very different
    # (RecursionError can be caused also by particular combinations of the input parameters GAP, MATCH and MISMATCH).

    try:
        alignments = recursive_traceback(df_moves, len(sequence2)-1, len(sequence1)-1, heuristic_token)
    except RecursionError:
        print("All possible alignments were too many to compute {}, so only the alignment which"
              "maximizes the letters aligned will be shown:".format('(even if heuristics was applied)' if heuristic_token else '')  )
        df_moves[df_moves == 6] = 2      # To 'maximize the letters aligned' means to choose the diagonal movement (move = 2) whenever is possible.
        df_moves[df_moves == 5] = 2
        df_moves[df_moves == 3] = 2
        df_moves[df_moves == 4] = random.choice([0,1])
        alignments = recursive_traceback(df_moves, len(sequence2)-1, len(sequence1)-1, heuristic_token)

    alignments_list = alignments.split(":")

    final_score = df_scores.iloc[-1,-1]
    string_separator = "\n\n_____________________________________________________________\n"

    if not quiet or output_filename is not None:
        for idx, alignment in enumerate(alignments_list):
            gaps_counter, identity_counter = 0, 0
            string_alignment_nr = 'Alignment nr. {}'.format(idx+1)
            if len(original_sequence1) >= len(original_sequence2):
                alignment_sequence1 = alignment[0::2]
                alignment_sequence2 = alignment[1::2]
            else:
                alignment_sequence1 = alignment[1::2]
                alignment_sequence2 = alignment[0::2]
            string_symbols = ''

            for idx, letter in enumerate(alignment_sequence1):
                if letter != '-' and alignment_sequence2[idx] != '-':
                    if letter == alignment_sequence2[idx]:
                        identity_counter += 1
                        string_symbols += '|'
                    else:
                        string_symbols += ':'
                else:
                    gaps_counter += 1
                    string_symbols += ' '

            string_alignment = ''
            n_whitespaces = 9 + len(str(len(alignment_sequence1)))
            for idx in range(0, len(alignment_sequence1), 100):
                string_alignment += ('\nSEQ1  ' + str(idx+1) + ' '*(n_whitespaces-6-len(str(idx+1))) + '{}   {}'.format(alignment_sequence1[idx:idx+100], min(idx+100, len(alignment_sequence1)))
                                     + '\n' + ' '*n_whitespaces + string_symbols[idx:idx+100]
                                   + '\nSEQ2  ' + str(idx+1) + ' '*(n_whitespaces-6-len(str(idx+1))) + '{}   {}'.format(alignment_sequence2[idx:idx+100], min(idx+100, len(alignment_sequence2))) +'\n' )

            string_score = '\nScore: {}'.format(final_score)
            identity_percentage = round( (identity_counter / len(alignment_sequence1))*100, 1)
            identity_string = '\nIdentities: {}/{} ({}%)'.format(identity_counter,len(alignment_sequence1),identity_percentage)
            gaps_percentage = round((gaps_counter / len(alignment_sequence1)) * 100, 1)
            gaps_string = '\nGaps: {}/{} ({}%)'.format(gaps_counter,len(alignment_sequence1),gaps_percentage)


            if not quiet:
                print(string_alignment_nr + string_alignment + string_score + identity_string + gaps_string + string_separator)

            if output_filename is not None:
                with open('./'+output_filename, 'a') as f:
                    f.write(string_alignment_nr + string_alignment + string_score + identity_string + gaps_string + string_separator + '\n')

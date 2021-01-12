from fetch_data import retrieve_ENST_seq, retrieve_ENST_overlap


class SeqParser:
    """Class containing useful methods for parsing a transcript sequence"""

    start_codon = ['ATG']
    stop_codon = ['TAG', 'TAA', 'TGA']

    def get_n_terminal_homologous_arms(self, ENST_ID: str, arm_len: int = 50):
        """Get region of cDNA that is translated

        Rule:
        1. Contains a start codon (ATG)
        2. Contains a stop codon (TAG, TAA, TGA)

        Note: Stop codon is included as part of CDS (https://m.ensembl.org/info/website/glossary.html)

        Refresher on terminology: https://www.biostars.org/p/65162/
        exon: sequence present in mature RNA (5'UTR + CDS + 3'UTR)
        cds: sequence present in mature RNA and codes for protein
        """
        seq = retrieve_ENST_seq(ENST_ID=ENST_ID)
        overlap = retrieve_ENST_overlap(ENST_ID=ENST_ID)
        cds = []
        exons = []
        if overlap:
            for i in overlap:
                if i.get('Parent') == ENST_ID:
                    if i.get('feature_type') == 'exon':
                        exons.append(i)
                    elif i.get('feature_type') == 'cds':
                        cds.append(i)
        # Set flag to True if reverse strand as the farthest position on the index
        # is the beginning of the sequence
        reverse = True if seq.strand == 'reverse' else False
        cds = sorted(cds, key=lambda x: x.get('start'), reverse=reverse)[0]
        # For development
        # print(reverse)
        # print(seq.left_seq_pos)
        # print(seq.right_seq_pos)
        # print(cds.get('start'))
        # print(cds.get('end'))
        # if reverse:
        #     end_index = abs(cds.get('start') - seq.right_seq_pos) - 1
        #     start_index = abs(cds.get('end') - seq.right_seq_pos)
        # else:
        #     start_index = abs(cds.get('start') - seq.left_seq_pos)
        #     end_index = abs(cds.get('end') - seq.left_seq_pos - 1)
        # print(start_index, end_index)
        # print(seq.seq[start_index:end_index])
        # print('  ', seq.seq[start_index+3:end_index])

        if reverse:
            # right most side is where sequence starts (reverse strand)
            print(cds.get('end'))
            left_arm = cds.get('end') - arm_len - 2, cds.get('end') - 3
            right_arm = cds.get('end') - 2, cds.get('end') + arm_len - 3
            return f"{seq.chromosome}:{left_arm[0]}-{left_arm[1]}", f"{seq.chromosome}:{right_arm[0]}-{right_arm[1]}"
        else:
            # left most side is where sequence starts (forward strand)
            left_arm = cds.get('start') - arm_len, cds.get('start') - 1
            right_arm = cds.get('start'), cds.get('start') + arm_len - 1
            return f"{seq.chromosome}:{left_arm[0]}-{left_arm[1]}", f"{seq.chromosome}:{right_arm[0]}-{right_arm[1]}"


if __name__ == "__main__":
    # Reverse strand
    print(SeqParser().get_n_terminal_homologous_arms(ENST_ID='ENST00000288602'))
    # Forward strand
    print(SeqParser().get_n_terminal_homologous_arms(ENST_ID='ENST00000413518'))
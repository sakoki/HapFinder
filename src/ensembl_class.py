from typing import Dict
import logging

logging.basicConfig(level=logging.DEBUG,
                    format="%(asctime)s %(levelname)s %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S")
logger = logging.getLogger(__name__)


class EnsemblSeq:
    """Class for ensembl sequence data"""

    nucleotide_pairs = {
        'A': 'T',
        'T': 'A',
        'G': 'C',
        'C': 'G',
        'a': 't',
        't': 'a',
        'g': 'c',
        'c': 'g'
    }

    def __init__(self,
                 seq: str = None,
                 enst_id: str = None,
                 molecule: str = None,
                 version: int = None,
                 description: str = None):
        self.seq = seq
        self.enst_id = enst_id
        self.molecule = molecule
        self.version = version
        self.description = description

        # Extract the following from description
        self.ref = None
        self.chromosome = None
        self.left_seq_pos = None
        self.right_seq_pos = None
        self.strand = None

        # Calculated fields
        self.seq_length = None
        self.start_index = None
        self.end_index = None
        self.ref_seq = None

        # Call methods for calculating fields
        self.__parse_desc()
        self.__get_seq_length()
        self.__get_index()
        self.__get_reference_seq()

    def __repr__(self):
        return f"{self.enst_id} - {self.description}"

    def parse_resp(self, resp: Dict):
        """Parse required fields from ENSEMBL api json response"""
        try:
            return EnsemblSeq(seq=resp.get('seq'),
                              enst_id=resp.get('id'),
                              molecule=resp.get('molecule'),
                              version=resp.get('version'),
                              description=resp.get('desc'))
        except Exception as e:
            logger.error(f"Failed to create EnsemblSeq class: {e}")
            return

    def __parse_desc(self):
        """Parse description field

        Example: 'chromosome:GRCh38:7:140734486:140924732:-1'
        """
        description = self.description
        if description:
            try:
                description = description.split(':')
                self.ref = description[1]
                self.chromosome = description[2]
                self.left_seq_pos = int(description[3])
                self.right_seq_pos = int(description[4])
                self.strand = 'forward' if int(description[5]) == 1 else 'reverse'
            except Exception as e:
                logger.error(f'Failed to parse desc field: {e}')

    def __get_seq_length(self):
        """Get sequence length"""
        left_pos = self.left_seq_pos
        right_pos = self.right_seq_pos
        if left_pos and right_pos:
            try:
                # starting 0 index
                self.seq_length = right_pos - left_pos + 1
            except Exception as e:
                logger.error(f'Failed calculating sequence length: {e}')

    def __get_index(self):
        """Get start and end index of sequence (starting at 0)"""
        left_pos = self.left_seq_pos
        right_pos = self.right_seq_pos
        if left_pos and right_pos:
            try:
                self.start_index = left_pos - left_pos
                self.end_index = right_pos - left_pos
            except Exception as e:
                logger.error(f'Failed calculating start and end index: {e}')

    def __get_reference_seq(self):
        """If strand is reverse, reference strand is the opposite

        Keep 5' --> 3' order
        """
        seq = self.seq
        strand = self.strand
        nucleotide_pairs = self.nucleotide_pairs
        ref_seq = ''
        if seq and strand == 'reverse':
            try:
                for i in range(len(seq)):
                    complement = nucleotide_pairs.get(seq[-1 * (i + 1)])
                    if isinstance(complement, str):
                        ref_seq += complement
                # Validity check - lengths must match
                if len(seq) == len(ref_seq):
                    self.ref_seq = ref_seq
                else:
                    raise ValueError("Reference seq must have the same length as the input seq")
            except Exception as e:
                logger.error(f'Failed building reference strand: {e}')

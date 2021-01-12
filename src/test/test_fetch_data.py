import unittest
import json
import os
from fetch_data import retrieve_ENST_seq, retrieve_ENST_overlap
from ensembl_class import EnsemblSeq


class TestClassConstruction(unittest.TestCase):

    # Initialize test vars
    seq = 'ATGCGC'
    ref_seq = 'GCGCAT'
    enst_id = 'ENST00000TESTING'
    molecule = 'dna'
    version = 11
    description = 'chromosome:GRCh38:7:100123456:100123461:-1'
    representation = 'ENST00000TESTING - chromosome:GRCh38:7:100123456:100123461:-1'

    # Extracted from description
    ref = 'GRCh38'
    chromosome = '7'
    left_seq_pos = 100123456
    right_seq_pos = 100123461
    strand = 'reverse'
    seq_length = 6
    start_index = 0
    end_index = 5

    # Load modified sample json response from ENSEMBL sequence API
    # Contains custom field for testing purposes
    sample_json = os.path.join(os.path.dirname(__file__), 'BRAF_sample.json')
    with open(sample_json, 'r') as f:
        sample_json = json.load(f)

    def test_1_manual_construction(self):
        """Construction by passing required fields manually"""
        seq = EnsemblSeq(
            seq=self.seq,
            enst_id=self.enst_id,
            molecule=self.molecule,
            version=self.version,
            description=self.description
        )
        self.assertEqual(seq.__repr__(), self.representation)
        self.assertEqual(seq.seq, self.seq)
        self.assertEqual(seq.enst_id, self.enst_id)
        self.assertEqual(seq.molecule, self.molecule)
        self.assertEqual(seq.version, self.version)
        self.assertEqual(seq.ref, self.ref)
        self.assertEqual(seq.chromosome, self.chromosome)
        self.assertEqual(seq.left_seq_pos, self.left_seq_pos)
        self.assertEqual(seq.right_seq_pos, self.right_seq_pos)
        self.assertEqual(seq.strand, self.strand)
        self.assertEqual(seq.seq_length, self.seq_length)
        self.assertEqual(seq.start_index, self.start_index)
        self.assertEqual(seq.end_index, self.end_index)
        self.assertEqual(seq.ref_seq, self.ref_seq)

    def test_2_parse_resp_construction(self):
        """Construction by parsing json response"""
        seq = EnsemblSeq().parse_resp(resp=self.sample_json)
        self.assertEqual(seq.__repr__(), self.sample_json.get('repr'))
        self.assertEqual(seq.seq, self.sample_json.get('seq'))
        self.assertEqual(seq.enst_id,  self.sample_json.get('id'))
        self.assertEqual(seq.molecule, self.sample_json.get('molecule'))
        self.assertEqual(seq.version, self.sample_json.get('version'))
        self.assertEqual(seq.ref, self.sample_json.get('ref'))
        self.assertEqual(seq.chromosome, self.sample_json.get('chromosome'))
        self.assertEqual(seq.left_seq_pos, self.sample_json.get('left_seq_pos'))
        self.assertEqual(seq.right_seq_pos, self.sample_json.get('right_seq_pos'))
        self.assertEqual(seq.strand, self.sample_json.get('strand'))
        self.assertEqual(seq.seq_length, self.sample_json.get('seq_length'))
        self.assertEqual(seq.start_index, self.sample_json.get('start_index'))
        self.assertEqual(seq.end_index, self.sample_json.get('end_index'))


class TestEnsemblAPICall(unittest.TestCase):

    # Test with BRAF gene
    enst_id = 'ENST00000288602'

    def test_1_call_api_and_format_resp(self):
        """Call Ensembl sequence API"""
        res = retrieve_ENST_seq(ENST_ID=self.enst_id)
        self.assertIsInstance(res, EnsemblSeq)

    def test_3_all_overlap_api(self):
        """Call Ensembl overlap API"""
        res = retrieve_ENST_overlap(ENST_ID=self.enst_id)
        cds = []
        exons = []
        for i in res:
            if i.get('Parent') == self.enst_id:
                if i.get('feature_type') == 'exon':
                    exons.append(i)
                elif i.get('feature_type') == 'cds':
                    cds.append(i)
        self.assertEqual(len(exons), 19)


if __name__ == "__main__":
    unittest.main(verbosity=2)
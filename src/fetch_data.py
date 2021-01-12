import requests
from ensembl_class import EnsemblSeq
import logging

logging.basicConfig(level=logging.DEBUG,
                    format="%(asctime)s %(levelname)s %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S")
logger = logging.getLogger(__name__)


class EnsemblAPIClient:
    """API client for querying ensembl database"""

    def __init__(self):
        self.seq_url = 'https://rest.ensembl.org/sequence/id/{}'
        self.overlap_url = 'https://rest.ensembl.org/overlap/id/{}'

    def get_sequence(self, ENST_ID: str, rtype: str = 'genomic'):
        """Get the requested ENST transcript sequence 5' -> 3

        API documentation: https://rest.ensembl.org/documentation/info/sequence_id
        '"""
        logger.debug(f'Requesting {rtype} sequence for ENST ID: {ENST_ID}')
        query = self.seq_url.format(ENST_ID)
        resp = requests.get(
            url=query,
            params={
                'type': rtype,
                'content-type': 'application/json'
            }
        )
        try:
            resp.raise_for_status()
        except Exception as e:
            logger.error(f'Ensembl sequence api request failed: {e}')
            return None
        if resp.status_code == 200:
            logger.info('Sequence api request success')
            return EnsemblSeq().parse_resp(resp=resp.json())

    def get_overlapping_sequences(self, ENST_ID: str):
        logger.debug(f'Requesting overlapping sequences for ENST ID: {ENST_ID}')
        query = self.overlap_url.format(ENST_ID)
        resp = requests.get(
            url=query,
            params={
                'feature': ['cds', 'exon'],
                'content-type': 'application/json'
            }
        )
        try:
            resp.raise_for_status()
        except Exception as e:
            logger.error(f'Ensembl overlap api request failed: {e}')
            return None
        if resp.status_code == 200:
            logger.info('Overlap api request success')
            return resp.json()


ensemblAPI = EnsemblAPIClient()


def retrieve_ENST_seq(ENST_ID: str, rtype: str = 'genomic'):
    return ensemblAPI.get_sequence(ENST_ID=ENST_ID, rtype=rtype)


def retrieve_ENST_overlap(ENST_ID: str):
    return ensemblAPI.get_overlapping_sequences(ENST_ID=ENST_ID)


if __name__ == "__main__":
    # Test sample request and construction of ensembl response object
    # ensemblAPI = EnsemblAPIClient()
    res = retrieve_ENST_seq(ENST_ID='ENST00000288602', rtype='genomic')
    print(res)
    print(res.molecule)
    print(res.version)
    print(res.ref)
    print(res.chromosome)
    print(res.left_seq_pos)
    print(res.right_seq_pos)
    print(res.strand)
    print('Sequence length: ', len(res.seq))
    print(res.seq[:500])
    res = retrieve_ENST_overlap(ENST_ID='ENST00000288602')
    print(len(res))
    t = []
    exons = []
    cds = []
    for i in res:
        if i.get('Parent') == 'ENST00000288602':
            if i.get('feature_type') == 'exon':
                exons.append(i)
            elif i.get('feature_type') == 'cds':
                cds.append(i)
    print(len(exons))
    print(len(cds))
    print(exons)
    print(cds)
    # print(res)

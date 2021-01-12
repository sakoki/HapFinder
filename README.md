# HapFinder

The following script generates homologous arms for making a DSB at the N-terminus of a transcript of interest.
Written in vanilla python version 3.8.6

File Overview
```angular2html
.
├── __init__.py
├── ensembl_class.py
├── fetch_data.py
├── seq_parser.py
└── test
   ├── __init__.py
   ├── __pycache__
   ├── BRAF_sample.json
   └── test_fetch_data.py
```

### Testing 
```angular2html
python3 -m unittest discover -v
```

### Sample Usage
```python
from seq_parser import SeqParser

SeqParser().get_n_terminal_homologous_arms(ENST_ID='ENST00000288602', arm_len=50)
# Output: ('7:140924651-140924700', '7:140924701-140924750')
```

#!/usr/bin/env python3

import os
import sys
import tempfile

sys.path.append('../../../')
from bin.trait_mapping.create_efo_table import *


class TestCreateEFOTableComponents:

    def setup_method(self):
        self.input_file = tempfile.NamedTemporaryFile(delete=False)
        self.output_file = tempfile.NamedTemporaryFile(delete=False)
        with open(self.input_file.name, 'w') as infile:
            infile.write('http://purl.obolibrary.org/obo/MONDO_0009796\n'
                         'http://www.orpha.net/ORDO/Orphanet_414\n')

    def teardown_method(self):
        os.remove(self.input_file.name)
        os.remove(self.output_file.name)

    def test_get_ols_details(self):
        """Check that the output of the script is adequate in terms of format; also check two example traits."""
        create_efo_table(self.input_file.name, self.output_file.name)
        with open(self.output_file.name) as outfile:
            out_lines = outfile.read().splitlines()
        assert len(out_lines) == 2, 'Should get two output lines for two input lines'
        trait_1, trait_2 = [out_line.split('\t') for out_line in out_lines]
        for trait in [trait_1, trait_2]:
            assert len(trait) == 22, 'Each output line must contain correct number of columns'

        # Trait 1
        assert trait_1[0] == 'ornithine aminotransferase deficiency'  # disease
        assert 'inherited retinal dystrophy' in trait_1[1]  # child of
        assert trait_1[7] == 'MeSH:C538071 || MeSH:D015799'  # MeSH cross-refs
        assert trait_1[8] == 'NCIT:C84744'  # NCIT cross-ref
        assert trait_1[12] == 'OMIM:258870'  # OMIM cross-ref
        assert trait_1[13] == 'DOID:1415'  # DOID cross-ref
        assert trait_1[15] == 'UMLS:C0018425'  # ULMS cross-ref
        assert trait_1[18] == 'Orphanet:414'  # Orphanet cross-ref
        assert trait_1[21] == 'MONDO:0009796'  # MONDO ref

        # Trait 2
        assert trait_2[0] == 'Gyrate atrophy of choroid and retina'  # disease
        assert trait_2[1] == 'Disease'  # child of
        assert 'is a very rare, inherited retinal dystrophy' in trait_2[2]  # description
        assert 'Ornithine aminotransferase deficiency' in trait_2[3]  # synonyms
        assert trait_2[7] == 'MeSH:C537132 || MeSH:C538071 || MeSH:D015799'  # MeSH cross-ref
        assert trait_2[8] == 'NCIT:C84744'  # NCIT cross-ref
        assert trait_2[12] == 'OMIM:258870'  # OMIM cross-ref
        assert trait_2[13] == 'DOID:1415'  # DOID cross-ref
        assert trait_2[15] == 'UMLS:C0018425 || UMLS:C0599035'  # ULMS cross-refs
        assert trait_2[18] == 'Orphanet:414'  # Orphanet cross-ref
        assert trait_2[21] == 'MONDO:0009796'  # MONDO ref

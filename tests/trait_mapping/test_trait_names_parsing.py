import os.path
import unittest

import eva_cttv_pipeline.trait_mapping.trait_names_parsing as trait_names_parsing


class TestGetTraitNames(unittest.TestCase):

    def test_trait_names_parsing(self):
        # Test file contains two records: one with a Pathogenic variant another with a Benign one. Only trait names
        # from the Pathogenic one should get into the
        test_filename = os.path.join(os.path.dirname(__file__), 'resources/variant_summary.tsv.gz')
        traits = trait_names_parsing.parse_trait_names(test_filename)
        self.assertEqual(
            sorted(traits),
            sorted(['breast-ovarian cancer, familial 1', 'hereditary breast and ovarian cancer syndrome',
             'hereditary cancer-predisposing syndrome', 'not provided'])
        )

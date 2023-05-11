import os

from cmat.clinvar_xml_io import ClinVarDataset


resources_dir = os.path.join(os.path.dirname(__file__), 'resources')


def test_dataset_write():
    input_file = os.path.join(resources_dir, 'test_clinvar_dataset.xml.gz')
    output_file = os.path.join(resources_dir, 'test_output.xml.gz')

    input_dataset = ClinVarDataset(input_file)
    input_dataset.write(output_file)
    output_dataset = ClinVarDataset(output_file)
    for r_in, r_out in zip(input_dataset, output_dataset):
        assert r_in.accession == r_out.accession

    if os.path.exists(output_file):
        os.remove(output_file)

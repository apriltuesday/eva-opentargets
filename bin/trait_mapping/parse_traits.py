import argparse

import cmat.trait_mapping.main as main


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Parse traits from ClinVar XML")
    parser.add_argument("-i", dest="input_filepath", required=True,
                        help="ClinVar XML dump file. One record per line.")
    parser.add_argument("-o", dest="output_traits_filepath", required=True,
                        help="path to output file for all traits for downstream processing")
    parser.add_argument("-u", dest="output_for_platform", required=False,
                        help="path to output file for all traits, for use with curation platform")
    args = parser.parse_args()
    main.parse_traits(args.input_filepath, args.output_traits_filepath, args.output_for_platform)

from eva_cttv_pipeline.clinvar_xml_io.clinvar_xml_io.hgvs_variant import HgvsVariant, SequenceType, VariantType


def test_lrg_sequence():
    h = HgvsVariant('LRG_274t1:c.(940+1_941-1)_(1186+1_1187-1)del')
    assert h.reference_sequence == 'LRG_274t1'
    assert h.sequence_type == SequenceType.CODING
    # Nothing else in this expression is supported yet
    assert h.variant_type == h.start == h.stop == h.precise_span() == None


def test_precise_deletion():
    h = HgvsVariant('NC_000011.8:g.67967551_67974756del')
    assert h.reference_sequence == 'NC_000011.8'
    assert h.sequence_type == SequenceType.GENOMIC
    assert h.variant_type == VariantType.DELETION
    assert h.start == 67967551
    assert h.stop == 67974756
    assert h.precise_span() == 67974756 - 67967551 + 1


def test_coordinate_pivots():
    h = HgvsVariant('NM_001256054.2(C9orf72):c.-45+63_-45+80GGGGCC(2_25)')
    assert h.reference_sequence == 'NM_001256054.2'
    assert h.sequence_type == SequenceType.CODING
    assert h.start == 63
    assert h.stop == 80
    assert h.repeat_sequence == 'GGGGCC'
    assert h.precise_span() == 80 - 63 + 1


def test_unknown_spans_with_known_endpoints():
    # No known span if start > stop (technically not well-formed HGVS)
    h = HgvsVariant('NC_000011.8:g.42_13del')
    assert h.start == 42
    assert h.stop == 13
    assert h.precise_span() is None

    # Similarly for coordinate pivots (which is well-formed HGVS, but not yet supported)
    h = HgvsVariant('NM_000548.3(TSC2):c.5068+27_5069-47dup34')
    assert h.start == 27
    assert h.stop == -47
    assert h.precise_span() is None

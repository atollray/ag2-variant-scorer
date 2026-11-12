"""Build and verify the synthetic DEL and tandem DUP used in the AG2 pilot."""

from alphagenome.data import genome

from reference import IndexedReference


# The affected interval is chr22:36202000–36202999 (1-based inclusive).
AFFECTED_REGION = genome.Interval("chr22", 36201999, 36202999)


def prepare_sv_examples(reference: IndexedReference) -> dict[str, genome.Variant]:
    region = AFFECTED_REGION
    flank = 32
    context = genome.Interval(region.chromosome, region.start - flank, region.end + flank)
    bases = reference.extract(context)
    before = bases[:flank]
    segment = bases[flank:-flank]
    after = bases[-flank:]
    if set(before[-1] + segment) - set("ACGT"):
        raise ValueError("The synthetic SV region or its anchor contains an ambiguous base.")

    # Anchor DEL immediately before the segment and insert the extra DUP copy after it.
    originals = {
        "DEL_1000bp": genome.Variant(
            chromosome=region.chromosome,
            position=region.start,
            reference_bases=before[-1] + segment,
            alternate_bases=before[-1],
        ),
        "TANDEM_DUP_1000bp": genome.Variant(
            chromosome=region.chromosome,
            position=region.end,
            reference_bases=segment[-1],
            alternate_bases=segment[-1] + segment,
        ),
    }
    expected_alt = {
        "DEL_1000bp": before + after,
        "TANDEM_DUP_1000bp": before + segment + segment + after,
    }

    cases = {}
    for case_id, original in originals.items():
        reference.verify_variant(original)
        variant = genome.normalize_variant(original, reference)
        reference.verify_variant(variant)
        start = variant.start - context.start
        end = variant.end - context.start
        if not 0 <= start < end <= len(bases):
            raise ValueError(f"Normalized {case_id} lies outside its verification interval.")
        reconstructed = bases[:start] + variant.alternate_bases + bases[end:]
        if reconstructed != expected_alt[case_id]:
            raise ValueError(f"Normalized {case_id} changes the intended ALT sequence.")
        cases[case_id] = variant
    return cases

"""Check the pilot SNV against an indexed reference without changing the index."""

from pathlib import Path
import warnings

from alphagenome.data import genome
from pyfaidx import Fasta, FastaIndexingError, FetchError


def verify_snv_reference(fasta: Path, variant: genome.Variant) -> None:
    """Check the REF base using the variant's zero-based reference interval."""
    if not variant.is_snv:
        raise ValueError("The pilot requires a single-nucleotide variant.")

    interval = variant.reference_interval
    # Adapted from the AG1 Reference reader and SNV branch of verified_normalize.
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            with Fasta(
                str(fasta),
                indexname=f"{fasta}.fai",
                as_raw=True,
                strict_bounds=True,
                sequence_always_upper=True,
                mutable=False,
                build_index=False,
                rebuild=False,
                duplicate_action="stop",
            ) as reference:
                observed = reference[interval.chromosome][interval.start:interval.end]
    except (FastaIndexingError, FetchError, KeyError, OSError, ValueError, Warning) as exc:
        raise ValueError("Could not read the REF base from the FASTA and its existing .fai.") from exc

    if len(observed) != 1 or observed not in "ACGT":
        raise ValueError("The reference base is ambiguous or missing.")
    if observed != variant.reference_bases:
        raise ValueError(f"REF mismatch for {variant}: FASTA contains {observed}.")

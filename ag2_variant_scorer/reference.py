"""Read an indexed reference without changing the FASTA or index."""

from pathlib import Path
import warnings

from alphagenome.data import genome
from pyfaidx import Fasta, FastaIndexingError, FetchError


class IndexedReference:
    def __init__(self, fasta: Path):
        self.fasta = fasta

    def extract(self, interval: genome.Interval) -> str:
        """Read a nonempty 0-based half-open interval for the SDK normalizer."""
        if not 0 <= interval.start < interval.end:
            raise ValueError(f"Invalid reference interval: {interval}")

        # Adapted from the AG1 Reference reader.
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("error")
                with Fasta(
                    str(self.fasta),
                    indexname=f"{self.fasta}.fai",
                    as_raw=True,
                    strict_bounds=True,
                    sequence_always_upper=True,
                    mutable=False,
                    build_index=False,
                    rebuild=False,
                    duplicate_action="stop",
                ) as reference:
                    sequence = reference[interval.chromosome][interval.start:interval.end]
        except (FastaIndexingError, FetchError, KeyError, OSError, ValueError, Warning) as exc:
            raise ValueError(f"Could not read {interval} from the FASTA and its existing .fai.") from exc

        if len(sequence) != interval.end - interval.start or set(sequence) - set("ACGTN"):
            raise ValueError(f"Invalid reference sequence or index for {interval}.")
        return sequence

    def verify_variant(self, variant: genome.Variant) -> None:
        """Check the entire REF allele before scoring."""
        observed = self.extract(variant.reference_interval)
        if "N" in observed:
            raise ValueError(f"Ambiguous REF sequence at {variant.reference_interval}.")
        if observed != variant.reference_bases:
            raise ValueError(f"REF mismatch at {variant.reference_interval}.")


def verify_snv_reference(fasta: Path, variant: genome.Variant) -> None:
    """Check the pilot SNV's REF base."""
    if not variant.is_snv:
        raise ValueError("The pilot requires a single-nucleotide variant.")
    IndexedReference(fasta).verify_variant(variant)

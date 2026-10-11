"""Pin the parse_strandedness.py semantics: the verbatim port of the legacy
check_strandedness awk mapping (fractions margin 0.4; secondstrand -> '1',
firststrand -> '2', unstranded -> '0')."""
import pytest

from parse_strandedness import parse_strand


def report(second_fraction, first_fraction):
    """An infer_experiment.py-shaped report (PE layout: values on the label lines)."""
    return (
        f'Fraction of reads explained by "1++,1--": {second_fraction}\n'
        f'Fraction of reads explained by "1+-,-1+": {first_fraction}\n'
    )


def test_secondstrand_maps_to_1():
    assert parse_strand(report(0.9, 0.05)) == "1"


def test_firststrand_maps_to_2():
    assert parse_strand(report(0.05, 0.9)) == "2"


def test_unstranded_maps_to_0():
    assert parse_strand(report(0.5, 0.5)) == "0"
    assert parse_strand(report(0.2, 0.2)) == "0"


def test_margin_is_exact():
    # 0.4 difference exactly -> inside the margin -> unstranded
    assert parse_strand(report(0.7, 0.3)) == "0"
    # just over the margin -> stranded call
    assert parse_strand(report(0.71, 0.3)) == "1"


def test_only_last_two_fractions_count():
    # First two values are equal (would give '0' if wrongly counted); the
    # last two decide, and their margin gives firststrand.
    text = report(0.1, 0.1) + 'Fraction of reads explained by "1++,1--": 0.95\n'
    assert parse_strand(text) == "2"


def test_malformed_report_raises():
    with pytest.raises(ValueError):
        parse_strand("no fractions here")

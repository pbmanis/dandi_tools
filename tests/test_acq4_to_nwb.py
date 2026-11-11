import pytest

from src.acq4_to_nwb import parse_age, protocol_type, select, sex_code, weight_with_units


@pytest.mark.parametrize(
    "agestr, expected",
    [
        ("P30D", ("P30D", 30.0, [])),
        ("30D", ("P30D", 30.0, [])),
        ("P95", ("P95D", 95.0, [])),
        ("P6W", ("P6W", 42.0, [])),
        ("P150+D", ("P150D/", None, ["lower bound"])),
        ("P150 +", ("P150D/", None, ["lower bound"])),
        ("P150 ish", ("P150D", None, ["approximate"])),
        ("P23?", ("P23D", None, ["approximate"])),
        ("", (None, None, [])),
        (None, (None, None, [])),
    ],
)
def test_parse_age(agestr, expected):
    assert parse_age(agestr) == expected


def test_parse_age_rejects_unknown_format():
    with pytest.raises(ValueError):
        parse_age("P30his")


@pytest.mark.parametrize(
    "sex, expected",
    [("M", "M"), ("male", "M"), ("F", "F"), ("Female", "F"), ("", "U"), (None, "U"), ("?", "U"), ("other", "O")],
)
def test_sex_code(sex, expected):
    assert sex_code(sex) == expected


@pytest.mark.parametrize(
    "input_weight, expected_output",
    [
        ("20.0", "20.0 g"),  # string, but no unit: should assume grams
        ("20 g", "20.0 g"),  # string with unit, just reformat
        ("20g", "20.0 g"),  # string with unit, but no space: should still work
        ("20 grams", "20.0 g"),  # string with full unit name: make standard format
        ("0.02 kg", "20.0 g"),  # string with different unit: should convert to grams
        (20, "20.0 g"),  # int without unit: should assume grams
        (20.0, "20.0 g"),  # float without unit: should assume grams
        (20.25, "20.25 g"),  # numbers and strings are formatted the same way
        ("20.25", "20.25 g"),
        ("", None),  # empty string: no weight
        (None, None),  # None: no weight
        (0, None),  # zero: not measured
        ("0 g", None),
        ("unknown", None),  # cannot be interpreted
    ],
)
def test_weight_conversion(input_weight, expected_output):
    assert weight_with_units(input_weight) == expected_output


@pytest.mark.parametrize(
    "protocol, expected_mode",
    [
        ("CCIV_long_HK_000", "CC"),
        ("Ic_LED_stim_20Hz_001", "CC"),
        ("Map_NewBlueLaser_IC_10Hz_000", "CC"),
        ("VCIV_simple_000", "VC"),
        ("Vc_LED_stim_10Hz_000", "VC"),
        ("VC_LED_Increasing_000", "VC"),
        ("Map_NewBlueLaser_VC_increase_1ms_HK_001", "VC"),
        ("Map_Camera_000", None),
        ("CC_Spont_ZeroI_000", None),
    ],
)
def test_protocol_type(protocol, expected_mode):
    ptype = protocol_type(protocol)
    assert (ptype[0] if ptype is not None else None) == expected_mode


def test_select_matches_protocol_type():
    protocols = ["CCIV_long_HK_000", "VCIV_simple_000", "Map_Camera_000", "CC_Spont_ZeroI_000"]
    assert select(protocols) == ["CCIV_long_HK_000", "VCIV_simple_000"]

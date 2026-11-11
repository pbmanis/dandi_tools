import numpy as np
import pytest

from src.display_nwb import parse_location


@pytest.mark.parametrize(
    "location, expected",
    [
        ("x=102555 um, y=7769 um (scanner coordinates)", [102555.0, 7769.0]),
        ("x=0 um, y=0 um (widefield, centered on the field of view)", [0.0, 0.0]),
        ("x=-12 um, y=5 um (scanner coordinates)", [-12.0, 5.0]),
        ("(np.int32(102554), np.int32(7769))", [102554.0, 7769.0]),  # older files
        ("[0, 0]", [0.0, 0.0]),  # older files
    ],
)
def test_parse_location(location, expected):
    np.testing.assert_array_equal(parse_location(location), expected)

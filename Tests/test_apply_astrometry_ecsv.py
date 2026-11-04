""" Applying a platepar to the measurements of an ECSV file.

Run with: python -m pytest Tests/test_apply_astrometry_ecsv.py
"""

import os

import pytest

np = pytest.importorskip("numpy")
astropy_table = pytest.importorskip("astropy.table")

from RMS.Astrometry.ApplyAstrometry import xyToRaDecPP
from RMS.Astrometry.ApplyAstrometryECSV import applyAstrometryECSV
from RMS.Astrometry.Conversions import trueRaDec2ApparentAltAz
from RMS.Formats.Platepar import Platepar
from RMS.Misc import getRmsRootDir


TEMPLATE = os.path.join(getRmsRootDir(), 'share', 'platepar_templates', 'template_generic_720p_4mm.cal')

ECSV = """# %ECSV 0.9
# ---
# datatype:
# - {name: datetime, datatype: string}
# - {name: x_image, unit: pix, datatype: float64}
# - {name: y_image, unit: pix, datatype: float64}
# - {name: integrated_pixel_value, datatype: int64}
# delimiter: ','
# schema: astropy-2.0
datetime,x_image,y_image,integrated_pixel_value
2026-09-05T06:00:00.000000,600.0,300.0,1000
2026-09-05T06:00:00.040000,620.0,310.0,1000
2026-09-05T06:00:00.080000,640.0,320.0,1000
"""


def test_azimuth_and_altitude_of_the_measurements(tmp_path):
    """ The azimuth and altitude written are those of the measured positions, converted from J2000. """

    platepar = Platepar()
    platepar.read(TEMPLATE, use_flat=None)
    platepar.star_list = []

    (tmp_path / "meteor.ecsv").write_text(ECSV)

    applyAstrometryECSV(str(tmp_path), "meteor.ecsv", None, platepar=platepar)

    table = astropy_table.QTable.read(str(tmp_path / "meteor.ecsv"), format='ascii.ecsv')

    time_data = [[2026, 9, 5, 6, 0, 0, milliseconds] for milliseconds in (0.0, 40.0, 80.0)]
    jd_data, ra_data, dec_data, _ = xyToRaDecPP(time_data, [600.0, 620.0, 640.0], [300.0, 310.0, 320.0],
                                                [1000, 1000, 1000], platepar, measurement=True)
    azim_data, alt_data = trueRaDec2ApparentAltAz(ra_data, dec_data, jd_data, platepar.lat, platepar.lon,
                                                  refraction=False)

    assert np.allclose(np.array(table['azimuth'], dtype=np.float64), azim_data, rtol=0, atol=1e-9)
    assert np.allclose(np.array(table['altitude'], dtype=np.float64), alt_data, rtol=0, atol=1e-9)

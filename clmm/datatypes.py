"""
Define the custom data type
"""

from astropy import table
from collections import namedtuple

"""
GCData: A namedtuple tying values with units to the metadata of where the values came from and thus how the values are used

Parameters
----------
creator: string, type of object (i.e. model, summarizer, inferrer) that made this data
specs: specifications of how the data was created
values: astropy table with column names and units

Notes
-----

"""
GCData = namedtuple('GCData', ['creator', 'specs', 'values'])

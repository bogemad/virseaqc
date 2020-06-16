from setuptools import setup
import os

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "virseaqc",
    version = "0.2",
    author = "Daniel Bogema",
    author_email = "daniel.bogema@dpi.nsw.gov.au",
    description = ("A combined pipeline for VIRal genome SEarching, Assembly, and Quality Control."),
    license = "GPL-3.0",
    keywords = "genomics assembly searching quality control",
    url = "https://github.com/bogemad/virseaqc",
    py_modules=['virseaqc'],
    scripts=['bin/virseaqc', 'bin/virseaqc_filter_host_reads'],
    long_description=read('README.md'),
)

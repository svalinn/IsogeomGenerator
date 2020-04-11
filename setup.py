from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="isosurface-generator",
    version="0.0.1",
    author="Kalin R. Kiesling",
    author_email="krkiesling@gmail.com",
    description="A package for generating meshed isosurface geometries",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/kkiesling/ww",
    packages=setuptools.find_packages(),
    python_requires='2.7',
)

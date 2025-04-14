from setuptools import setup, find_packages

VERSION = '0.0.2'

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='spectool',
    version=VERSION,
    description='FWHM determination tool',
    author="Joao Bento",
    author_email="jpsbento@gmail.com",
    packages=find_packages('src'),
    package_dir={'': 'src'},
    long_description=long_description,
    long_description_content_type="text/markdown",
    entry_points={},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        "matplotlib",
        "scipy"
    ]
)

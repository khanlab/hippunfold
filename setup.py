import setuptools
import json

with open("README.rst", "r") as fh:
    long_description = fh.read()

with open('hippunfold/pipeline_description.json', 'r') as fh:
    pipeline = json.load(fh)
    name = pipeline['GeneratedBy'][0]['Name']
    description = pipeline['Name']
    version = pipeline['GeneratedBy'][0]['Version']
    url = pipeline['GeneratedBy'][0]['CodeURL']
    author = pipeline['GeneratedBy'][0]['Author']
    author_email = pipeline['GeneratedBy'][0]['AuthorEmail']
 
setuptools.setup(
    name=name,
    version=version,
    author=author,
    author_email=author_email,
    description=description,
    long_description=long_description,
    long_description_content_type="text/x-rst",
    url=url,
    packages=setuptools.find_packages(),
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    entry_points={'console_scripts': [
        'hippunfold=hippunfold.run:main'
    ]},
    install_requires=[
        "snakebids>=0.2.0",
        "snakemake>=5.28.0",
        "nnunet>=1.6.6",
        "pandas",
        "nibabel",
        "numpy",
        "scipy"
    ],
    python_requires='>=3.7'
)

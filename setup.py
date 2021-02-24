
with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="HCRProbeDesign", # Replace with your own username
    version="0.0.1",
    author="Loyal A. Goff",
    author_email="author@example.com",
    description="Probe Design tool for Hybridization Chain Reaction",
    long_description=long_description,
    long_description_content_type="text/markdown",
    license="MIT",
    url="https://github.com/pypa/sampleproject",
    packages=setuptools.find_packages(),
    install_requires=[
        'primer3-py',
        'biopython',
        'beautifulsoup4',
    ]
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
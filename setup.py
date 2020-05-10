import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="byase",
    version="1.0.2",
    author="Lili Dong",
    author_email="ncjllld@hit.edu.cn",
    description="A library that uses Bayesian inference to identify gene-level and isoform-level ASE",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ncjllld/byase",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    entry_points={
        'console_scripts': [
            'byase=byase.run:main'],
    },
)

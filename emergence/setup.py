from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open("requirements.txt", "r", encoding="utf-8") as fh:
    requirements = fh.read().splitlines()

setup(
    name="quantum-emergence",
    version="1.0.0",
    author="Quantum Matter Research Group",
    author_email="contact@quantum-matter.example",
    description="A package for simulating and analyzing quantum systems with emergent time",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/quantum-matter/emergence",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Mathematics",
    ],
    python_requires=">=3.8",
    install_requires=requirements,
    entry_points={
        "console_scripts": [
            "quantum-emergence=emergence.main:main",
        ],
    },
)

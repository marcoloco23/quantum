from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open("requirements.txt", "r", encoding="utf-8") as fh:
    requirements = fh.read().splitlines()
    # Filter out comments and platform-specific dependencies
    requirements = [
        req for req in requirements if not req.startswith("#") and ";" not in req
    ]

setup(
    name="quantum-matter",
    version="1.0.0",
    author="Marc Sperzel",
    author_email="marcsperzel1@gmail.com",
    description="A comprehensive quantum simulation framework for emergent time, Ising models, and SYK models",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/marcoloco23/quantum",
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
        "Intended Audience :: Science/Research",
    ],
    python_requires=">=3.8",
    install_requires=requirements,
    entry_points={
        "console_scripts": [
            "quantum-emergence=emergence.main:main",
            "quantum-app=app:main",
        ],
    },
)

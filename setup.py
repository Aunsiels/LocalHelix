from setuptools import setup, find_packages
import pathlib

here = pathlib.Path(__file__).parent.resolve()

# Get the long description from the README file
long_description = (here / "README.md").read_text(encoding="utf-8")

setup(
    name="LocalHelix",
    version="0.1.0",
    description="This program analyzes your DNA to find interesting insights. Do not use for medical advice and"
                "always consult your doctor.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Aunsiels/LocalHelix",
    author="Julien Romero",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research ",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3 :: Only",
    ],
    keywords="dna, analysis, snp, clinvar, snpedia, gwas",
    package_dir={"": "localhelix"},
    packages=find_packages(where="localhelix"),
    python_requires=">=3.7, <4",
    install_requires=["requests", "beautifulsoup4", "pandas", "pickledb", "tqdm", "wikitextparser"],
    entry_points={  # Optional
        "console_scripts": [
            "sample=analyzer:main",
        ],
    },
    project_urls={  # Optional
        "Bug Reports": "https://github.com/Aunsiels/LocalHelix/issues",
        "Funding": "https://www.paypal.com/donate/?hosted_button_id=GR3D64Y7S7TU2",
        "Source": "https://github.com/Aunsiels/LocalHelix",
    },
)

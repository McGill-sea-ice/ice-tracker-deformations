import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="2021_SeaIceDeformation-bdu002",
    version="0.0.1",
    author="Beatrice Duval",
    author_email="beatrice.duval@mail.mcgill.ca",
    description="A package for the computation of arctic sea-ice deformations from icetracker data (Sentinel-1 and RCM).",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://gitlab.science.gc.ca/bdu002/2021_SeaIceDeformation",
    project_urls={
        "Bug Tracker": "https://gitlab.science.gc.ca/bdu002/2021_SeaIceDeformation/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
    ],
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    python_requires=">=3.6",
)
import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="cryoorigami-turalaksel", # Replace with your own username
    version="1.0.0",
    author="Tural Aksel",
    author_email="turalaksel@gmail.com",
    description="Python toolset for cryo-EM data analysis",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/douglaslab/pyOrigamiEM",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.5',
)

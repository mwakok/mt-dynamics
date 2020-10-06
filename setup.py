from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="mtdynamics", # Replace with your own username
    version="0.1.0",
    author="Florian Huber and Maurits Kok",
    author_email="f.huber@esciencecenter.nl",
    description="MT dynamics simulation",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/florian-huber/mt-dynamics",
    packages=find_packages(exclude=['*tests*']),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=[
        "matplotlib",
        "numpy",
        "pandas",
        "parameters",
        "scipy",
    ],
    setup_requires=['pytest-runner'],
    tests_require=['pytest'],
)

from setuptools import setup, find_packages

setup(
    name="vqf",
    version="0.1.0",
    description="This repository contains implementation of the algorithm presented in the article \"Variational Quantum Factoring\", by Eric R. Anschuetz, Jonathan P. Olson, Alán Aspuru-Guzik, Yudong Cao.",
    author="Michal Stechly",
    license="MIT",
    packages=find_packages(),
    python_requires=">=3.10",
    install_requires=open("requirements.txt").readlines(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.10",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)

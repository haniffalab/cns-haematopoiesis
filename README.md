# Standard Project

A template repo for the standard Haniffa Lab project

## About

We want to learn about the development of different bones (skull/femur/hipbone/shoulder/spine) as sites of haematopoesis. 

How will we study this? 

As a first step, look at the distribution of the progenitors binned by gestation and by organ. Compare with the bone marrow, liver, yolk sac cell composition by gestation. 

Hypothesis:   haematopoeisis in the skull begins before it begins at the other bone sites.  The skull is an independent site of haematopoesis. 

In the lymphoid organs (liver, BM), the macrophages and mast cells populate the haematopoetic organ, and interact with the microenvironment (stromal cells) to prepare the organ to start the production of progenitors. So, we’ll expect the infiltration and hence abundance of mast cells + macrophages at an earlier time-point, and an increase in the progenitors at later time points. Will we observe this in the skull? What’s our negative control? Knee as we don’t expect haematopoeisis there. 

How will the transcriptional profile/abundance of the immune cells+ its microenvironment change across time? What does this tell about the change in the function of these cells across gestation. Can answer these questions by looking at the change in the ratio of the cell types. 


Primary lymphoid organs:  bone marrow, liver and  thymus
Secondary lymphoid organs: lymph nodes, the spleen, mucosa-associated lymphoid tissue (MALT) 

### Project Team

Dr L. Ipsum, Newcastle University ([lorem.ipsum@newcastle.ac.uk](mailto:lorem.ipsum@newcastle.ac.uk))  
Professor D. Sit Amet, XY University ([d.sit.amet@newcastle.ac.uk](mailto:d.sit.amet@example.com))

### Contact

C. Adipiscing
Newcastle University  
([consectetur.adpiscing@newcastle.ac.uk](mailto:consectetur.adpiscing@newcastle.ac.uk))

## Built With

This section is intended to list the frameworks and tools you're using to develop this software. Please link to the home page or documentatation in each case.

[Framework 1](https://something.com)  
[Framework 2](https://something.com)  
[Framework 3](https://something.com)

## Getting Started

### Prerequisites

Any tools or versions of languages needed to run code. For example specific Python or Node versions. Minimum hardware requirements also go here.

### Installation

How to build or install the applcation.

### Running Locally

How to run the application on your local system.

### Running Tests

How to run tests on your local system.

## Deployment

### Local

Deploying to a production style setup but on the local system. Examples of this would include `venv`, `anaconda`, `Docker` or `minikube`.

### Production

Deploying to the production system. Examples of this would include cloud, HPC or virtual machine.

## Usage

Any links to production environment, video demos and screenshots.

## Roadmap

- [x] Initial Research
- [ ] Minimum viable product <-- You are Here
- [ ] Alpha Release
- [ ] Feature-Complete Release

## Contributing

Please install `nbconvert` with `pip` or `conda` and run the following command in your cloned repo to enable automatic filtering of Jupyter notebook output when creating commits.

```sh
echo '[filter "strip-notebook-output"]\nclean = "jupyter nbconvert --ClearOutputPreprocessor.enabled=True --ClearMetadataPreprocessor.enabled=True --to=notebook --stdin --stdout --log-level=ERROR"' >> .git/config
```

If a notebook's output must be preserved in the repo then name the notebook with the extensions `.output.ipynb` so it doesn't go through the filter.

### Main Branch

Protected and can only be pushed to via pull requests. Should be considered stable and a representation of production code.

### Dev Branch

Should be considered fragile, code should compile and run but features may be prone to errors.

### Feature Branches

A branch per feature being worked on.

https://nvie.com/posts/a-successful-git-branching-model/

## License

## Citation

Please cite the associated papers for this work if you use this code:

```
@article{xxx2023paper,
  title={Title},
  author={Author},
  journal={arXiv},
  year={2023}
}
```

## Acknowledgements

This work was funded by a grant from the UK Research Councils, EPSRC grant ref. EP/L012345/1, “Example project title, please update”.

**Copyright &copy; 2021, Sebastian Andreß**\
All rights reserved. Please find the license [here](https://github.com/sebastianandress/Slicer-SurfaceFragmentsRegistration/blob/master/LICENSE.md).

Please cite the corresponding paper when using this algorithm for publications:

    @article{SimilaritySubgroups-Andress,
        author      = {Andreß, Sebastian and Achilles, Felix and Bischoff, Jonathan and Calvalcanti Kußmaul, Adrian and Böcker, Wolfgang and Weidert},
        title       = {A Method for Finding High Accuracy Surface Zones on 3D Printed Bone Models},
        journal     = {Computers in Biology and Medicine},
        publisher   = {Elsevier},
        year        = {2021},
        doi         = {https://doi.org/10.1016/j.compbiomed.2021.104590}
    }


![Header](/Resources/header.png)

# Surface Fragments Registration for 3D Slicer

## Introduction
This tool is intended to find different surface error types. Unless traditional methods, it tries to find correlating subgroups of vertices, that show high accuracy (within this subgroup) and uses those for registration.

It tries to close the gap between rigid and deformable registration and validation methods. As described in the paper mentioned above, this method is tested for four different types of model errors.

The main purpose of this module is to find high accuracy surface zones (above a certain area size), that show low deviation (within a certain threshold) when taken for themselves, regardless of whether they deviate significantly from the entire model.

![Screenshot](/Resources/screenshot1.png)

The algorithm can be divided into four sections:
A optional pre-registration, B automatic search of an initialization subgroup and initial registration with it, C further cleaning of the subgroup, and D final registration and deviation calculation for this registration. It is an iterative process, till most vertices of the source model were assigned to a subgroup.

Parameters mentioned in the diagram can also be set in the module.

![Flowchart](/Resources/flowchart.png)

## Example

### Boxes - Minimal Working Example
Two boxes are compared with each other. Both are identical, but one is broken in the middle into two pieces. The algorithm recognizes this fact and registers both halves of the source box with the target box, thus creating a "Similarity Subgroup" for each half. For each subgroup considered as such, the deviation is zero.

![ExampleOutput](/Resources/exampleBoxes.png)

### Hemipelvis - Realistic Example
Two segmentations of a pelvis (e.g. one before and one after 3D printing) are compared. One shows a fracture with a displaced fragment (Fig A, red, source), the other is a correct match to the CT (Fig B, green, target).

The algorithm automatically detects the fragment (Fig. B & D, yellow), registers it to the target and calculates the surface deviation of the source to the target per registration (Fig. C & E).

![ExampleOutput](/Resources/exampleHemipelvis.png)


## How to install
The Extension will be available in the [Extension Manager](http://slicer.kitware.com/midas3/slicerappstore/extension/view?extensionId=330842) of 3D Slicer soon.
To install it manually, please follow the description on the official [3D Slicer page](https://www.slicer.org/wiki/Documentation/Nightly/Developers/FAQ/Extensions). 

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
This extension tries to close the gap between rigid and deformable registration. The method is only applicable to surface models. It is an extension of the [ICP-Registration](https://de.wikipedia.org/wiki/Iterative_Closest_Point_Algorithm), where in a complementary step contiguous high-accuracy zones ("fragments") are searched, separated, and registered individually. This happens iteratively, so that as result in each case for the surface fragment suitable registrations are found.

Background of this method is that it can come with the production of models (e.g. in the 3D printing) to errors, with which an entire region is consistently shifted from the remaining model (e.g. if a part of the model breaks off/warps). Both areas in themselves are correct and consistent, only the combination of both cannot be used. This method aims to automatically identify both areas and create a registration with surface deviation in respectively, in order to make the model still usable, depending on the purpose of the application.

As can be seen in the above paper, the use case is specifically aimed at FDM 3D printing of bone models. In this case, it often happens that a fracture fragment shifts in comparison to the rest of the model, e.g. due to improper detachment from the print bed. Since this is not relevant intraoperatively after reduction, a method was required which allows for automatic, fragment-wise validation.

![Screenshot](/Resources/screenshot1.png)

The algorithm can be divided into four sections:
A optional pre-registration, B automatic search of an initialization subgroup and initial registration with it, C further cleaning of the subgroup, and D final registration and deviation calculation for this registration. It is an iterative process, till most vertices of the source model were assigned to a subgroup.

Parameters mentioned in the diagram can also be set in the module.

![Flowchart](/Resources/flowchart.png)

## Parameters
- Inputs
    - `Source Model` (*S*): Surface fragments will be calculated as well as distance scalars will be applied to this model.
    - `Target Model` (*T*): The ground truth model, here the Source Model will be compared to.
- Settings
    - `Initialization Candidate Radius` (*ζ*): The algorithm uses an initialization region, which is randomly chosen. Its size is determined by this radius. The size should be smaller than the expected fragment size. Too large regions will lead to less fragments, too small regions may lead to multiple inadequate fragments.
    - `Minimal Fragment Area` (*ε*): No fragment will be smaller than this area. The algorithm will continue till no continuous surface part of the model larger than this value is assigned to a fragment. To small values may result in an infinite loop (terminated by `Maximal Iterations` value).
    - `Cutoff Deviation` (*δ*): All fragments, that are generated, will show a maximum continuous deviation of this value. A smaller value typically results in more fragments. Too small values may result in multiple inadequate fragments.
- Advanced
    - `Initialization LM-Registration` with Source and Target Landmark Fiducials: It is recommended to perform a pre-registration, optionally using the modules own landmark registration method. Please select the Source and Target Fiducials Node. Both need to contain >3 fiducials each in corresponding order.
    - `Initialization Iterations` (*n*): An initialization region is randomly chosen and selected by its smallest overall deviation. The number of the random zones per iteration can be selected. Too small values may result in bad initial registration, large values may lead to a longer computation time.
    - `Opening Width` (*χ*): To remove thin bridging vertex lines that remain after thresholding due to the three dimensional overlapping, an opening method similar to the morphological operator in image processing is used. The original operator works by first eroding and then dilating a mask image. This value variates the number of vertex rows eroded and dilated in this process. Using small values may result in bridging vertices between fragments, too large values may smooth the outlines of the fragments.
    - `Maximal Iterations`: For some models it might happen that the algorithm gets caught in a infinite loop, depending on the used settings. This is a technical selector terminating the iterative loop. It is not possible for the algorithm to find more fragments as number of iterations selected here.
    - `Create Transformations`: Create a transformation node for each fragment registration.
    - `Mark Deviations`: Generate a scalar showing the surface deviation of each source models vertex for each registration.
    - `Mark Fragments`: Generate a scalar showing highlighting each vertex of a respective fragment.
- Outputs
    - `Select Fragment`: After running this method, the generated transformations and surface deviation scalars can be applied here in a fast forward manner. It may cause errors if nodes are renamed or deleted manually.


## Example

Two segmentations of a pelvis (e.g. one before and one after 3D printing) are compared. One shows a fracture with a displaced fragment (screenshot above, red, Source Model), the other is a correct match to the CT (screenshot above, green, Target Model).

The algorithm automatically detects the fragment, registers it to the target and calculates the surface deviation of the source to the target per registration. Deviation shown by heatmap cold to hot.

![ExampleOutput](/Resources/exampleOutput.gif)


## How to install
The Extension will be available in the [Extension Manager](http://slicer.kitware.com/midas3/slicerappstore/extension/view?extensionId=330842) of 3D Slicer soon.
To install it manually, please follow the description on the official [3D Slicer page](https://www.slicer.org/wiki/Documentation/Nightly/Developers/FAQ/Extensions). 

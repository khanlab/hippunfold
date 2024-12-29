# WIP - workflow for updating atlases to new v2 specs

This example is run from the existing multihist atlas folder.

Still needs updates for cifti structure labels etc (will add this in once we have an updated wb_command/wb_view in our deps)

We currently have templates and atlases in hippunfold (selected by --template and --atlas), where the former selects the volume-based template used for volumetric registration, and the latter selects the unfolded-space subfields (and metrics for unfolded registration). I am proposing the latter become purely surface-based, so instead of providing labels and metrics in the volumetric unfolded space, we use giftis for everything, and this now becomes surface-templates instead of atlases. 
